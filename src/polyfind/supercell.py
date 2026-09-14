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
* :class:`SupercellEnergy` delegates complete energy, Cartesian torsion, valence,
  charge flux, Ewald/induced polarization, field and all derivatives to the packer.
  Its cells may carry independently perturbed nuclei and a full triclinic lattice;
  :meth:`SupercellEnergy.check_against_packer` checks replication of the same model.
  Rigid-charge :func:`cell_dipole` does not include flux or induced moments; use
  :meth:`SupercellEnergy.dipole` for the actual complete-model moment.

What it deliberately does **not** do is search.  It is an evaluator; the caller chooses the
stagger pattern and drives the optimiser (see ``examples/copolymer_staggered.py``).
"""
from __future__ import annotations

import numpy as np

from .pack import E_PER_A2_TO_C_PER_M2, CellEnergyTerms, CrystalPacker, PeriodicChain
from .periodic_geometry import ChainLayout, PlacedCell
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
    layout = ChainLayout.from_cell(cell)
    repeat_charges = np.asarray(charges)
    if (repeat_charges.shape != (cell.n_per_chain,) or repeat_charges.dtype.kind not in "iuf"
            or not np.isfinite(repeat_charges).all()):
        raise ValueError("charges must be a finite real vector for one declared chain repeat")
    q = np.empty(cell.n_atoms,dtype=float)
    q[layout.atoms] = repeat_charges
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
class SupercellEnergy:
    """Complete packer model on independently placed homogeneous chain repeats.

    The packer owns topology, charges, all Hamiltonian terms and derivatives.
    Cell metadata declares source/local atom order and each chain's reversal;
    no canonical setting-angle reconstruction or contiguous-block assumption.
    """

    def __init__(self, packer: CrystalPacker, chain: PeriodicChain | None = None):
        if chain is not None and chain is not packer.chain:
            raise ValueError("chain topology belongs to the supplied packer; use its chain")
        self.pk = packer

    @property
    def chain(self):
        return self.pk.chain

    def evaluate(self, cell: Cell, *, images=None):
        """Complete energy/decomposition/derivatives, in declared source order."""
        return self.pk.evaluate_chain_cell(cell.coords,cell.lattice,ChainLayout.from_cell(cell),
                                           images=images)

    def images(self, cell: Cell, extra: int = 0) -> np.ndarray:
        """The geometry owner's minimum-fractional image representatives."""
        return np.asarray(list(PlacedCell(cell.coords,cell.lattice).image_indices(self.pk.rc,extra)),
                          dtype=np.int64)

    def column_correction(self, cell: Cell) -> float:
        """UNHALVED bonded pair diagnostic, already included in pair energy.

        Computed by the complete owner on these exact nuclei, charges and
        effective image indices. Energy itself applies scales directly; it
        never subtracts this large diagnostic from an unscaled pair sum.
        """
        return self.evaluate(cell).terms.pair_correction

    def pair_energy(self, cell: Cell, images=None) -> float:
        return self.evaluate(cell,images=images).terms.pair

    def terms(self, cell: Cell, images=None) -> CellEnergyTerms:
        """All complete-model terms, including Cartesian torsion/valence/field."""
        return self.evaluate(cell,images=images).terms

    def energy(self, cell: Cell) -> float:
        return self.evaluate(cell).terms.total

    def energy_per_monomer(self, cell: Cell) -> float:
        return self.evaluate(cell).terms.per_monomer

    def energy_and_grad(self, cell: Cell):
        """(E, source-order atom gradient, all nine lattice derivatives)."""
        result = self.evaluate(cell)
        return result.terms.total,result.grad_coords.copy(),result.grad_lattice.copy()

    def dipole(self, cell: Cell):
        """(permanent, induced) cell dipoles on these actual nuclei, in e.A."""
        return self.pk.chain_dipole(cell.coords,cell.lattice,ChainLayout.from_cell(cell))

    def energy_is_converged(self, cell: Cell, extra: int = 1) -> float:
        """Absolute pair-energy difference after adding cutoff image shells.

        This verifies real-space cutoff coverage only, not Ewald k-space,
        force convergence, a relaxed minimum, stability or physical validity.
        """
        return abs(self.pair_energy(cell,self.images(cell,extra))
                   - self.pair_energy(cell,self.images(cell)))

    def check_against_packer(self, params, nx: int = 1, ny: int = 1) -> float:
        """Per-monomer replication check of the SAME complete model."""
        cell = build_cell(self.chain,params,nx=nx,ny=ny,packer=self.pk)
        mine = self.energy_per_monomer(cell)
        theirs = float(self.pk.energy(np.asarray(params,dtype=float)[None])[0]) / (
            self.pk.n_chains*self.chain.n_monomers)
        return abs(mine-theirs)
