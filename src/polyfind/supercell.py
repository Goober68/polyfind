"""Cells of many *independent* chains: staggered construction, energy, dipole.

:class:`polyfind.pack.CrystalPacker` places one or two chains per cell and reaches larger
cells by tiling (:func:`polyfind.topology.build_cell`).  A tiling is an exact replication,
so it adds atoms without adding packing degrees of freedom: a chain cannot sit at a
different axial height from its own image, because it *is* its own image.  For a copolymer
that is not a cosmetic limitation.  An ``11 VDF : 1 VDCN`` chain carries one bulky nitrile
per twelve monomers, and a ``2x2x1`` tiling puts **every** chain's nitrile at the same
axial height, turning eight isolated substituents into a continuous plane of them that the
whole structure has to clear.

This module adds the missing degree of freedom and nothing else:

* :func:`staggered_cell` builds the tiled cell and then slides each chain along its own
  axis by a whole number of monomers.  A slide is a rigid translation of a neutral chain,
  so it changes no bond length, no bond angle and no torsion, and (see
  :func:`cell_dipole`) it leaves the cell dipole *exactly* where it was -- a staggered
  antipolar cell is still exactly antipolar, which is why the stagger and the polarity can
  be chosen independently.
* :class:`SupercellEnergy` evaluates the lattice energy of any such cell with the packer's
  own pair potential, so a staggered cell and the aligned cell it came from are scored by
  the same function.  It is a direct sum over every periodic image within the cutoff rather
  than a re-derivation: :meth:`SupercellEnergy.check_against_packer` asserts that it
  reproduces :meth:`polyfind.pack.CrystalPacker.energy` per monomer on the packer's own
  two-chain cell, which it must, since the energy per monomer of a pair potential is
  invariant under tiling.

What it deliberately does **not** do is search.  It is an evaluator; the caller chooses the
stagger pattern and drives the optimiser (see ``examples/copolymer_staggered.py``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import backend as bk
from .ewald import exclusion_correction
from .forcefield import COULOMB
from .pack import E_PER_A2_TO_C_PER_M2, CrystalPacker, PeriodicChain, _erfc_pos
from .topology import Cell, build_cell

__all__ = [
    "STAGGER_PATTERNS",
    "SupercellEnergy",
    "cell_density",
    "cell_dipole",
    "cell_polarization",
    "cell_to_cif",
    "column_layout",
    "stagger_pattern",
    "staggered_cell",
]


# --------------------------------------------------------------- the lateral layout
def column_layout(nx: int = 2, ny: int = 2, k_chains: int = 2) -> list[tuple[int, int, int]]:
    """``(ix, iy, k)`` of each chain, in the order :func:`polyfind.topology.build_cell` emits.

    ``k`` indexes the chain within the packer's own cell: ``k = 0`` sits at the cell origin
    and ``k = 1`` at ``(a + b) / 2``, i.e. half a lattice vector along each of the two
    transverse axes.  So in units of ``b / 2`` a chain's transverse height is
    ``2 iy + k``, which is the layer index the ladder pattern below ramps along.
    """
    return [(ix, iy, k) for ix in range(nx) for iy in range(ny) for k in range(k_chains)]


def stagger_pattern(name: str, nx: int = 2, ny: int = 2, k_chains: int = 2,
                    monomers: int = 12) -> np.ndarray:
    """Per-chain axial offset in whole monomers, for one of :data:`STAGGER_PATTERNS`.

    The offsets are integers because the copolymer's backbone is a PVDF backbone: every
    monomer is geometrically the same, so sliding a chain by a whole monomer maps its
    carbon/hydrogen/fluorine skeleton onto itself and moves only the substituent that makes
    the monomer a comonomer.  A fractional slide would move the skeleton too and is a
    different (and strictly worse packed) object; ``examples/copolymer_staggered.py``
    measures the skeleton's periodicity rather than assuming it.
    """
    if name not in STAGGER_PATTERNS:
        raise ValueError(f"unknown stagger pattern {name!r} (known: {sorted(STAGGER_PATTERNS)})")
    f = STAGGER_PATTERNS[name]
    return np.array([f(ix, iy, k, monomers) % monomers
                     for ix, iy, k in column_layout(nx, ny, k_chains)], dtype=int)


# Each entry maps a chain's lateral position to its axial offset in monomers.  ``m`` is the
# number of monomers in the repeat, so every pattern is stated as a fraction of the repeat
# and none of them hard-codes twelve.
STAGGER_PATTERNS = {
    # every chain at the same height: the 2x2x1 tiling, kept here so the aligned reference
    # goes through exactly the same code path as the staggered cells
    "aligned": lambda ix, iy, k, m: 0,
    # a uniform ladder: one quarter of the repeat per b/2 layer.  Affine in the transverse
    # height, so it is a c-glide rather than an arbitrary decoration, and it is the
    # arrangement that maximises the *minimum* axial offset over every close column pair
    # (see examples/copolymer_staggered.py, which measures that rather than claiming it).
    "ladder": lambda ix, iy, k, m: (2 * iy + k) * (m // 4),
    # the ladder plus one monomer per step along a, so all eight chains sit at eight
    # different heights: the most axial slots occupied, at a smaller nearest-neighbour offset
    "spread": lambda ix, iy, k, m: (2 * iy + k) * (m // 4) + ix,
    # the two sublattices of the centred cell half a repeat apart, leaving both the 4.80 A
    # b-axis and the 10.5 A a-axis neighbours aligned.  Not an independent structure: the
    # packer's own ``dz`` already offsets chain 2 against chain 1, so this is the aligned cell
    # at a different point of its own registry sweep -- which makes it a check that the search
    # is reproducible rather than a candidate.
    "antiphase": lambda ix, iy, k, m: k * (m // 2),
    # half a repeat per step along a, so the two sheets stacked along the *long* axis are
    # offset and both near shells stay aligned.  The other half of the decomposition by
    # neighbour shell.
    "sheet": lambda ix, iy, k, m: ix * (m // 2),
}


def staggered_cell(chain: PeriodicChain, params, stagger, nx: int = 2, ny: int = 2,
                   packer: CrystalPacker | None = None) -> Cell:
    """:func:`polyfind.topology.build_cell`, then slide chain ``cid`` by ``stagger[cid]`` monomers.

    ``stagger`` is one integer per chain of the tiled cell, in monomers along ``+c``.  The
    slide is applied after the tiling, so the cell's lattice, composition, atom order and
    chain blocking are byte-for-byte those of the aligned cell; only the ``z`` of whole
    chains moves.  ``stagger = 0`` everywhere reproduces :func:`build_cell` exactly.

    Coordinates are **not** wrapped, for the same reason ``build_cell`` does not wrap them:
    each chain's crystallographic repeat is kept whole so that the cell dipole is a property
    of the cell rather than of where the origin was put.
    """
    cell = build_cell(chain, params, nx=nx, ny=ny, packer=packer)
    m = np.asarray(stagger, dtype=float).ravel()
    if m.shape != (cell.n_chains,):
        raise ValueError(f"stagger has {m.shape[0]} entries but the cell has {cell.n_chains} chains")
    dz = m * (chain.c / chain.n_monomers)
    shift = np.zeros((cell.n_atoms, 3))
    shift[:, 2] = dz[cell.chain_of]
    cell.coords = cell.coords + shift
    return cell


# ------------------------------------------------------------------ cell properties
def cell_dipole(cell: Cell, charges) -> np.ndarray:
    """``sum_i q_i r_i`` (e.A) over the cell as built, with ``charges`` one repeat's worth.

    Every repeat unit of these polymers is neutral, so this is independent of the origin and
    of which image of each chain is used -- and therefore **independent of the stagger**: a
    rigid translation ``d`` of a neutral chain changes its moment by ``(sum q) d = 0``.  A
    staggered cell's polarity is exactly its aligned parent's.  ``examples`` measures it
    anyway; the helper that builds antipolar cells has shipped two defects in this area.
    """
    q = np.tile(np.asarray(charges, dtype=float).ravel(), cell.n_chains)
    if q.shape[0] != cell.n_atoms:
        raise ValueError(f"{q.shape[0]} charges for {cell.n_atoms} atoms")
    tot = float(q.sum())
    if abs(tot) > CrystalPacker.NEUTRAL_TOL:
        raise ValueError(f"the cell carries {tot:+.3e} e: sum_i q_i r_i is not a dipole moment")
    return np.einsum("n,nc->c", q, cell.coords)


def cell_polarization(cell: Cell, charges) -> np.ndarray:
    """``mu / V`` in C/m^2; the rigid-ion polarization of the point-charge model."""
    return cell_dipole(cell, charges) / cell.volume * E_PER_A2_TO_C_PER_M2


def cell_density(cell: Cell, mass_per_chain: float) -> float:
    """g/cm^3 from the cell volume and ``mass_per_chain`` (one repeat's formula mass)."""
    return cell.n_chains * mass_per_chain / (cell.volume * 0.602214)


def cell_to_cif(cell: Cell, title: str = "polyfind") -> str:
    """P1 CIF of a :class:`~polyfind.topology.Cell`, in fractional coordinates.

    :func:`polyfind.pack.to_cif` for a cell that did not come from a two-chain
    :class:`~polyfind.pack.PackResult`.  Fractional coordinates are taken modulo one, which
    wraps the overhanging pendant atoms a CIF reader would otherwise place outside the box;
    the extended-XYZ file keeps each chain whole, so use that one for anything that reads the
    cell dipole off the coordinates.
    """
    lat = np.asarray(cell.lattice, dtype=float)
    lengths = np.linalg.norm(lat, axis=1)
    ang = [float(np.degrees(np.arccos(np.clip(
        lat[(i + 1) % 3] @ lat[(i + 2) % 3] / (lengths[(i + 1) % 3] * lengths[(i + 2) % 3]),
        -1.0, 1.0)))) for i in range(3)]
    frac = np.linalg.solve(lat.T, cell.coords.T).T % 1.0
    lines = [f"data_{title}",
             f"_cell_length_a {lengths[0]:.4f}",
             f"_cell_length_b {lengths[1]:.4f}",
             f"_cell_length_c {lengths[2]:.4f}",
             f"_cell_angle_alpha {ang[0]:.3f}",
             f"_cell_angle_beta {ang[1]:.3f}",
             f"_cell_angle_gamma {ang[2]:.3f}",
             "_symmetry_space_group_name_H-M 'P 1'",
             "loop_",
             "_atom_site_label",
             "_atom_site_type_symbol",
             "_atom_site_fract_x",
             "_atom_site_fract_y",
             "_atom_site_fract_z"]
    for i, (e, f) in enumerate(zip(cell.elements, frac)):
        lines.append(f"{e}{i + 1} {e} {f[0]:.5f} {f[1]:.5f} {f[2]:.5f}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------- energy
@dataclass
class SupercellTerms:
    """The pieces of :meth:`SupercellEnergy.energy`, so a caller can report them apart."""

    pair: float  # the pair kernel over every image, bonded scales applied
    torsion: float  # the Fourier torsion term, one repeat's worth per chain
    ewald: float  # the reciprocal-space sum, zero for a dsf packer
    exclusion: float  # what the bonded exclusions take back out of that sum
    n_chains: int
    n_monomers: int

    @property
    def total(self) -> float:
        return self.pair + self.torsion + self.ewald + self.exclusion

    @property
    def per_monomer(self) -> float:
        return self.total / (self.n_chains * self.n_monomers)


class SupercellEnergy:
    """:class:`polyfind.pack.CrystalPacker`'s potential, on a cell of arbitrarily many chains.

    Construct it from the packer whose potential you want; every setting that packer carries
    -- the Lennard-Jones cutoff convention, the damping, ``eps_r``, the charges, and whether
    the electrostatics are the truncated damped-shifted-force sum or an Ewald sum -- is taken
    from its own pair tables rather than rebuilt, so the two agree by construction and not by
    coincidence.  :meth:`check_against_packer` verifies that on the packer's own cell.

    Not supported, and refused rather than ignored: a ``valence`` term or a ``charge_flux``
    (both of which make the energy a function of the chain's internal geometry, which this
    rigid-chain sum does not vary), and an applied ``field``.
    """

    def __init__(self, packer: CrystalPacker, chain: PeriodicChain | None = None):
        if packer.valence is not None:
            raise ValueError("a valence term is a per-chain constant for a rigid chain, but it is "
                             "not this module's to guess at: pass a packer without one")
        if packer.charge_flux is not None:
            raise ValueError("charge flux makes the charges a function of the chain geometry; this "
                             "rigid multi-chain sum does not carry it")
        if getattr(packer, "_field_on", False):
            raise ValueError("an applied field is not carried here; clear it with set_field(None)")
        self.pk = packer
        self.chain = chain if chain is not None else packer.chain
        self.rc = float(packer.rc)
        self.rc2 = float(packer.rc2)
        self.alpha = float(packer.alpha)
        self.n = int(packer.n)
        to = lambda Z: np.asarray(bk.to_numpy(Z), dtype=float)  # noqa: E731
        self._A, self._B = to(packer._A_nn), to(packer._B_nn)
        self._qq, self._qqf = to(packer._qq_nn), to(packer._qqf_nn)
        self._const = to(packer._const_nn)
        self._scales = dict(packer._scale_nn_base)  # z-image offset -> (n, n) bonded scale
        self._q_repeat = np.asarray(packer._q_cell[: self.n], dtype=float)
        self._ewald = packer._ewald
        self._cache: dict[int, tuple] = {}

    # --- tables -------------------------------------------------------------------
    def _tiled(self, n_chains: int):
        t = self._cache.get(n_chains)
        if t is None:
            tile = lambda Z: np.tile(Z, (n_chains, n_chains))  # noqa: E731
            t = self._cache[n_chains] = (tile(self._A), tile(self._B), tile(self._qq),
                                         tile(self._qqf), tile(self._const))
        return t

    def _v(self, r2, A, B, qq, qqf, const):
        """The pair potential on a flat selection of squared distances, all inside the cutoff.

        Character for character :meth:`polyfind.pack.CrystalPacker._pair_energy`'s expression,
        and fed from that packer's own tables, so the only difference between the two sums is
        the order the terms are added in: measured at 4e-10 kcal/mol per monomer on a 592-atom
        cell (:meth:`check_against_packer`), which is accumulated rounding over some 10^5 pairs
        rather than a difference in the potential.
        """
        rr = np.sqrt(r2)
        ir6 = 1.0 / (r2 * r2 * r2)
        return ir6 * (A * ir6 - B) + qq * (_erfc_pos(self.alpha * rr) / rr) + qqf * rr + const

    def images(self, cell: Cell, extra: int = 0) -> np.ndarray:
        """The integer images a minimum-image-reindexed pair sum has to cover.

        Pair separations are taken in the minimum-image fractional frame first
        (:meth:`pair_energy`), exactly as :meth:`polyfind.ewald.Ewald.terms` does, so the
        images needed are ``ceil(rc / h + 1/2)`` per axis rather than however many the
        unwrapped coordinates happen to span.  On a 592-atom cell that is 45 images instead
        of 343, and it is a tighter sum over the *same* set of physical pairs, not a smaller
        one: :meth:`energy_is_converged` adds ``extra`` shells and checks nothing moves.
        """
        lat = cell.lattice
        V = abs(np.linalg.det(lat))
        widths = V / np.linalg.norm(np.cross(np.roll(lat, -1, axis=0), np.roll(lat, -2, axis=0)), axis=1)
        r = np.ceil(self.rc / widths + 0.5).astype(int) + int(extra)
        return np.array([(ia, ib, ic) for ia in range(-r[0], r[0] + 1)
                         for ib in range(-r[1], r[1] + 1) for ic in range(-r[2], r[2] + 1)], dtype=float)

    # --- the sum -------------------------------------------------------------------
    def column_correction(self, cell: Cell) -> float:
        """What the bonded exclusions take out of an unscaled sum, over the whole cell.

        ``sum_t sum_{k,p,q} (s_k[p,q] - 1) v(|r_tp - r_tq - k c|)``: the only pairs a bond path
        can join are pairs in the same chain at the same transverse image, so this is the whole
        difference between the scaled sum :meth:`pair_energy` wants and the unscaled sum it
        actually evaluates.  For a chain placed antiparallel the repeat-closing bond runs
        towards ``-c`` and the stack is read at ``-k``.

        **Evaluated on the cell's own atoms, chain by chain, rather than once on the ideal
        chain.** Every chain is the same rigid object up to an isometry, so in exact arithmetic
        one evaluation times the chain count would do -- but the terms being corrected are
        1-2 and 1-3 pairs at 1.1 to 2.5 A, where ``r^-12`` is of order ``10^5`` kcal/mol, and
        the correction is of order ``10^2``.  That is a difference of large numbers, and it only
        cancels if the coordinates the two halves see are the *same* coordinates.  They are not
        for a cell read back from a file: extended XYZ carries eight decimals, and on a pair
        whose potential has a gradient of ``10^6`` kcal/mol/A a 5e-9 A difference is worth
        millikelvin-scale energy -- measurably, about 2e-3 kcal/mol per monomer on a 592-atom
        cell.  Taking both halves from the cell makes the cancellation exact whatever
        coordinates the cell carries, at the cost of a 74x74 sum per chain.
        """
        lat2 = np.asarray(cell.lattice[2], dtype=float)
        K = int(np.ceil(self.rc / float(np.linalg.norm(lat2)))) + 1
        n = cell.n_per_chain
        tot = 0.0
        for t in range(cell.n_chains):
            X = np.asarray(cell.coords[t * n : (t + 1) * n], dtype=float)
            D = X[:, None, :] - X[None, :, :]
            for k in range(-K, K + 1):
                S = self._scales.get(-k if cell.reversed_of[t] else k)
                if S is None:
                    continue
                d = D - k * lat2
                r2 = d[:, :, 0] ** 2 + d[:, :, 1] ** 2 + d[:, :, 2] ** 2
                sel = (r2 < self.rc2) & (r2 > 1e-8) & (S != 1.0)
                if not sel.any():
                    continue
                idx = np.nonzero(sel)
                v = self._v(r2[idx], self._A[idx], self._B[idx], self._qq[idx], self._qqf[idx],
                            self._const[idx])
                tot += float(((S[idx] - 1.0) * v).sum())
        return tot

    def pair_energy(self, cell: Cell, images=None) -> float:
        """``0.5 sum_{i,j,T} s_ij(T) v(|r_i - r_j - T|)`` over every image within the cutoff.

        Evaluated as an *unscaled* sum over every pair plus one correction
        (:meth:`column_correction`), because the bonded scales touch only pairs inside one
        chain at one transverse image.  That removes the need to know each pair's true image
        index, which is what lets the separations be taken in the minimum-image frame.
        """
        A, B, qq, qqf, const = self._tiled(cell.n_chains)
        X = np.asarray(cell.coords, dtype=float)
        lat = cell.lattice
        hinv = np.linalg.inv(lat)
        S = X @ hinv
        dS = S[:, None, :] - S[None, :, :]
        base = (dS - np.round(dS)) @ lat  # (N, N, 3) minimum-image separations
        ijk = self.images(cell) if images is None else np.asarray(images, dtype=float)
        # Most of those images cannot contain a pair at all, and finding out costs three
        # reductions rather than a 592 x 592 distance array: the separations live inside the
        # box [lo, hi], so the closest any pair in image ``shift`` can come is the distance
        # from the origin to the shifted box.  Conservative (it ignores the correlation
        # between components), exact arithmetic, and on an orthorhombic 592-atom cell it
        # leaves 3 images of 45 -- the rest are genuinely empty, which
        # :meth:`energy_is_converged` confirms by keeping them.
        lo = base.reshape(-1, 3).min(axis=0)
        hi = base.reshape(-1, 3).max(axis=0)
        tot = 0.0
        for shift in ijk @ lat:
            near = np.clip(0.0, lo + shift, hi + shift)
            if float(near @ near) >= self.rc2:
                continue
            dx = base[:, :, 0] + shift[0]
            dy = base[:, :, 1] + shift[1]
            dz = base[:, :, 2] + shift[2]
            r2 = dx * dx + dy * dy + dz * dz
            sel = (r2 < self.rc2) & (r2 > 1e-8)
            if not sel.any():
                continue
            idx = np.nonzero(sel)
            tot += float(self._v(r2[idx], A[idx], B[idx], qq[idx], qqf[idx], const[idx]).sum())
        return 0.5 * (tot + self.column_correction(cell))

    def terms(self, cell: Cell, images=None) -> SupercellTerms:
        """Every piece of the energy of ``cell``, in kcal/mol per cell."""
        pair = self.pair_energy(cell, images)
        tors = self.pk.torsion_energy(self.chain.dihedrals) * cell.n_chains
        e_ew = e_ex = 0.0
        if self._ewald is not None:
            q = np.tile(self._q_repeat, cell.n_chains)
            e_ew = float(self._ewald.terms(cell.coords, q, cell.lattice).total)
            # Per chain, on the cell's own atoms and the cell's own axial period, for the
            # reason :meth:`column_correction` gives: the correction is a small difference of
            # large 1/r terms and only cancels against coordinates it actually shares.
            cz = float(cell.lattice[2, 2])
            K = int(np.ceil(self.rc / cz)) + 1
            stack = self.pk._scale_column_np(K)
            n = cell.n_per_chain
            e_ex = sum(float(exclusion_correction(
                cell.coords[t * n : (t + 1) * n], self._q_repeat, cz,
                np.flip(stack, axis=0) if cell.reversed_of[t] else stack,
                pref=COULOMB / self.pk.eps_r)[0]) for t in range(cell.n_chains))
        return SupercellTerms(pair=pair, torsion=tors, ewald=e_ew, exclusion=e_ex,
                              n_chains=cell.n_chains, n_monomers=self.chain.n_monomers)

    def energy(self, cell: Cell) -> float:
        """Lattice energy of ``cell`` in kcal/mol per cell."""
        return self.terms(cell).total

    def energy_per_monomer(self, cell: Cell) -> float:
        return self.terms(cell).per_monomer

    # --- self-checks ----------------------------------------------------------------
    def energy_is_converged(self, cell: Cell, extra: int = 1) -> float:
        """``|E(images) - E(images + extra shells)|``: zero when the image bound is enough."""
        return abs(self.pair_energy(cell, self.images(cell, extra))
                   - self.pair_energy(cell, self.images(cell)))

    def check_against_packer(self, params, nx: int = 1, ny: int = 1) -> float:
        """``|this per monomer - packer per monomer|`` on the packer's own cell.

        The energy per monomer of a pair potential is invariant under tiling, so for any
        ``nx, ny`` this must be zero to rounding.  It is the whole basis for comparing a
        staggered cell with the aligned one: the same function scores both.
        """
        cell = build_cell(self.chain, params, nx=nx, ny=ny, packer=self.pk)
        mine = self.energy_per_monomer(cell)
        theirs = float(self.pk.energy(np.asarray(params, dtype=float)[None])[0]) / (
            self.pk.n_chains * self.chain.n_monomers)
        return abs(mine - theirs)
