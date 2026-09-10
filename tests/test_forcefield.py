import warnings

import numpy as np
import pytest

from polyfind.chain import build_chain, build_chain_batch
from polyfind.forcefield import (
    Frame,
    SimpleFF,
    _reversal_angle_residual,
    _reversal_image3,
    _reversal_images,
    erfc_approx,
    fit_ris,
)
from polyfind.polymers import PE, PVDF, THREE_STATE, get_polymer
from scipy.special import erfc


def test_erfc_approx():
    x = np.linspace(-3, 3, 61)
    assert np.allclose(erfc_approx(x), erfc(x), atol=2e-7)


def test_energy_batch_matches_single():
    ff = SimpleFF()
    rng = np.random.default_rng(0)
    structs = [build_chain(PVDF, rng.uniform(-180, 180, 8)) for _ in range(6)]
    eb = ff.energy_batch(structs)
    for s, e in zip(structs, eb):
        assert ff.energy(s) == pytest.approx(e, rel=1e-8)


def test_pe_trans_below_gauche_and_pentane_effect():
    ff = SimpleFF()
    rep = fit_ris(PE, ff, step=15.0, n_monomers=8, scan="dense")
    m = rep.model
    T, GP, GM = 0, 1, 2
    assert abs(m.first_order[0, T]) < 1e-9
    assert 0.2 < m.first_order[0, GP] < 2.0  # gauche costs a fraction of a kcal/mol
    assert m.first_order[0, GP] == pytest.approx(m.first_order[0, GM])  # mirror symmetry
    assert m.second_order[0, GP, GM] > 1.0  # pentane effect: G+G- strongly disfavoured
    assert abs(m.second_order[0, T, T]) < 1e-9
    assert rep.n_evaluations == 1 + 24 + 24 * 24


def test_pvdf_fit_has_two_pair_types_and_is_mirror_symmetric():
    ff = SimpleFF()
    rep = fit_ris(PVDF, ff, step=20.0, n_monomers=5)
    m = rep.model
    assert m.B == 2 and m.second_order.shape == (2, 3, 3)
    mir = np.array(THREE_STATE.mirror)
    assert np.allclose(m.second_order, m.second_order[:, mir][:, :, mir])
    # the two pair matrices differ (CH2- vs CF2-centred)
    assert not np.allclose(m.second_order[0], m.second_order[1])


# ------------------------------------------------ reflection-with-reversal, orders 1-3
def _raw_fit(name, third_order=True, step=20.0):
    """An unsymmetrised fit, so the relations can be measured rather than assumed."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the chiral polymers warn once on build_chain
        return fit_ris(get_polymer(name), SimpleFF(), step=step, n_monomers=3,
                       third_order=third_order, symmetrize=False)


@pytest.mark.parametrize("name", ["pvdf", "pvdc", "cfe", "cdfe"])
def test_third_order_terms_obey_reflection_with_reversal(name):
    """The third-order image of E(phi_1..phi_N) == E(-phi_N..-phi_1), measured.

    Reversal reverses the order of a triple as well as mirroring its states, and shifts
    the bond type by one more than the pair term does, so

        e3[b, x, y, z] == e3[(c - 2 - b) % B, m(z), m(y), m(x)].

    PVDF and PVDC are achiral controls: any true symmetry of the model has to hold there
    too, and does.  The near misses -- mirroring without reversing the triple, reversing
    without mirroring, and the right form at the wrong bond-type shift -- are all wrong
    for the chiral pair, which is what makes this a measurement and not a relabelling.
    """
    p = get_polymer(name)
    rep = _raw_fit(name)
    e1, e2, e3 = rep.model.first_order, rep.model.second_order, rep.model.third_order
    B, m = p.bonds_per_repeat, np.array(THREE_STATE.mirror)
    *_, c = _reversal_images(e1, e2, THREE_STATE.mirror, B)
    assert c == 1  # for B = 2 the shift is decided by the first- and second-order fit

    _, resid = _reversal_image3(e3, THREE_STATE.mirror, B, c)
    assert resid < 0.5, resid  # grid noise: 0.00 (pvdf, cfe, cdfe), 0.02 (pvdc)

    plain = float(np.abs(e3 - e3[:, m][:, :, m][:, :, :, m]).max())
    wrong_shift = _reversal_image3(e3, THREE_STATE.mirror, B, c + 1)[1]
    no_reverse = float(np.abs(e3 - e3[[(c - 2 - b) % B for b in range(B)]][:, m][:, :, m][:, :, :, m]).max())
    no_mirror = float(np.abs(e3 - np.transpose(e3[[(c - 2 - b) % B for b in range(B)]], (0, 3, 2, 1))).max())
    if p.is_chiral:
        assert plain > 5.0 and wrong_shift > 5.0 and no_reverse > 5.0 and no_mirror > 5.0
    else:
        assert plain < 0.5  # an achiral chain has the plain mirror as well


@pytest.mark.parametrize("name", ["cfe", "cdfe"])
def test_a_chiral_fit_symmetrises_all_three_orders(name):
    """The default fit averages a chiral model over the operation it just validated."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = fit_ris(get_polymer(name), SimpleFF(), step=20.0, n_monomers=3, third_order=True).model
    B, m = model.B, np.array(model.states.mirror)
    *_, c = _reversal_images(model.first_order, model.second_order, model.states.mirror, B)
    assert _reversal_image3(model.third_order, model.states.mirror, B, c)[1] < 1e-9
    # and the real G+/G- asymmetry is still there, in the third order too
    assert np.abs(model.third_order - model.third_order[:, m][:, :, m][:, :, :, m]).max() > 5.0


@pytest.mark.parametrize("name", ["pvdf", "pvdc", "cfe", "cdfe"])
def test_state_angles_obey_reflection_with_reversal(name):
    """``adapt_angles``' mirror pairing is valid for a chiral fit too, and is checked.

    The 1-D scan of a bond of type ``b`` is the scan of type ``(c - b) % B`` read at the
    negated angle, so ``arg1[b, s] == -arg1[(c - b) % B, m(s)]``.  Averaging that over the
    bond types (a bijection) is what makes the adapted angles of a mirror pair equal and
    opposite -- for a chiral chain as much as an achiral one, even though the *per bond
    type* minima are wildly asymmetric there (80 deg apart for CFE and CDFE).
    """
    p = get_polymer(name)
    rep = _raw_fit(name, third_order=False)
    B = p.bonds_per_repeat
    *_, c = _reversal_images(rep.model.first_order, rep.model.second_order, THREE_STATE.mirror, B)
    assert _reversal_angle_residual(rep.argmin1, THREE_STATE, B, c) == 0.0  # on a dense grid, exactly

    gp, gm = THREE_STATE.names[1], THREE_STATE.names[2]
    per_type = max(abs(rep.argmin1[b][gp] + rep.argmin1[b][gm]) for b in range(B))
    assert (per_type > 10.0) if p.is_chiral else (per_type == 0.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = fit_ris(get_polymer(name), SimpleFF(), step=20.0, n_monomers=3).model
    assert model.states.angles[1] == -model.states.angles[2]  # aggregated: exactly paired


def test_energy_coords_matches_energy_batch_list_and_single():
    ff = SimpleFF()
    rng = np.random.default_rng(3)
    dih = rng.uniform(-180, 180, size=(5, 8))
    template, coords = build_chain_batch(PVDF, dih)
    e_tuple = ff.energy_batch((template, coords))
    structs = [build_chain(PVDF, dih[m]) for m in range(5)]
    e_list = ff.energy_batch(structs)
    np.testing.assert_allclose(e_tuple, e_list, atol=1e-8)
    for m, e in enumerate(e_tuple):
        assert ff.energy(structs[m]) == pytest.approx(e, rel=1e-8)


def test_fit_ris_unknown_scan_raises():
    with pytest.raises(ValueError):
        fit_ris(PE, SimpleFF(), scan="bogus")


def test_dense_scan_batched_matches_unbatched_reference():
    """The batched dense scan (build_chain_batch + energy_batch) must reproduce the
    original per-point build_chain/calc.energy_batch(list) path exactly."""
    ff = SimpleFF()
    rep = fit_ris(PVDF, ff, step=20.0, n_monomers=5, scan="dense")
    grid = rep.grid
    N = max(5 * 2, 2 * 2 + 6)
    base = np.full(N, 180.0)
    ref = ff.energy(build_chain(PVDF, base))
    j0 = (N // 2 - 1) - ((N // 2 - 1) - 0) % 2
    structs = []
    for phi in grid:
        d = base.copy()
        d[j0] = phi
        structs.append(build_chain(PVDF, d))
    expected = ff.energy_batch(structs) - ref
    np.testing.assert_allclose(rep.scan1[0], expected, atol=1e-9)


@pytest.mark.parametrize("polymer,n_monomers", [(PE, 8), (PVDF, 6)])
def test_adaptive_matches_dense_10deg_fit(polymer, n_monomers):
    """The adaptive (coarse-to-fine) fit must reproduce the dense (default) 10-degree
    scan's first- and second-order energies: never meaningfully worse, often better since
    it is not confined to a grid, with the same mirror symmetry, the same basin-angle
    signs, and several times fewer calculator evaluations.

    One matrix entry is excluded from the ``second_order`` comparison here: PVDF's
    CF2-centred G+/G- pair (bond type 1), whose basin has no interior minimum -- see
    ``test_adaptive_pvdf_g_plus_g_minus_has_no_interior_minimum``, which checks it
    explicitly instead. Everywhere else -- including PE's own G+/G- pair (whose basin
    *does* have an interior minimum, just close to the basin edge) and PVDF's CH2-centred
    G+/G+ pair (a real, fairly sharp bowl around (86, 86) deg that the dense grid samples
    only at (90, 80), verified by a 2-degree brute-force scan) -- adaptive matches dense
    within [-0.6, +0.1] kcal/mol: it is never meaningfully worse, and where it is better
    by more than refinement noise, that reflects genuine sub-grid curvature dense's grid
    is too coarse to resolve, not an error.
    """
    dense = fit_ris(polymer, SimpleFF(), step=10.0, n_monomers=n_monomers, third_order=True, scan="dense")
    adap = fit_ris(polymer, SimpleFF(), step=10.0, n_monomers=n_monomers, third_order=True, scan="adaptive")

    d1, a1 = dense.model.first_order, adap.model.first_order
    d2, a2 = dense.model.second_order, adap.model.second_order

    mask = np.ones_like(d2, dtype=bool)
    if polymer is PVDF:
        gp, gm = THREE_STATE.index("G+"), THREE_STATE.index("G-")
        mask[1, gp, gm] = mask[1, gm, gp] = False

    assert np.all(a1 <= d1 + 0.1)
    assert np.all(a1 >= d1 - 0.5)
    assert np.all(a2[mask] <= d2[mask] + 0.1)
    assert np.all(a2[mask] >= d2[mask] - 0.6)

    mir = np.array(polymer.states.mirror)
    assert np.allclose(a1, a1[:, mir], atol=1e-8)
    assert np.allclose(a2, a2[:, mir][:, :, mir], atol=1e-8)

    for b in range(polymer.bonds_per_repeat):
        assert np.sign(adap.argmin1[b]["G+"]) == np.sign(dense.argmin1[b]["G+"])
        assert np.sign(adap.argmin1[b]["G-"]) == np.sign(dense.argmin1[b]["G-"])

    assert dense.n_evaluations >= 2.5 * adap.n_evaluations


def test_adaptive_pvdf_g_plus_g_minus_has_no_interior_minimum():
    """PVDF's CF2-centred (bond type 1) G+/G- pair is the one place adaptive and dense
    disagree by more than grid/refinement noise, and it is expected, not a bug: this
    basin has no interior critical point near the G+/T boundary. Scanning phi from 100 to
    140 deg (crossing straight through the nominal boundary at 120) at each point's own
    best psi gives a smooth, monotonically *decreasing* profile (E ~= 7.8, 5.1, 4.88 (at
    phi=120), 4.7, 4.04 (at phi=140), ...; no kink, no minimum). So the true infimum of a
    correctly bounded search of the open G+/G- basin (matching ``_basins``' tie-break --
    see ``_basin_bounds``) sits at its own edge, at ~4.85 kcal/mol; excluding the single
    boundary point cannot raise the infimum of a continuous, monotonic function
    approaching it. The dense step=10 deg scan reports a much higher ~8.17 kcal/mol only
    because it is too coarse to sample anywhere near that edge (its neighbouring grid
    points are 110 deg, E ~= 10.2, and 120 deg itself, excluded by the same tie-break).
    This is asserted explicitly, as expected behaviour, rather than folded into a widened
    tolerance in the general accuracy test above.
    """
    dense = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=6, scan="dense")
    adap = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=6, scan="adaptive")
    gp, gm = THREE_STATE.index("G+"), THREE_STATE.index("G-")

    d = dense.model.second_order[1, gp, gm]
    a = adap.model.second_order[1, gp, gm]
    assert d == pytest.approx(8.17, abs=0.2)
    assert a == pytest.approx(4.85, abs=0.2)
    assert a - d < -3.0  # a large, expected improvement -- not refinement noise

    # neither coordinate sits exactly on the G+/T (or G-/T) tie boundary that _basins
    # would reassign to T, but adaptive does ride right up against one of them
    phi, psi = adap.argmin2[1][("G+", "G-")]
    assert phi != 120.0 and psi != -120.0
    assert min(abs(phi - 120.0), abs(psi + 120.0)) < 1e-3


# ------------------------------------------------------- the settable parameters
def _probe_dihedrals(n, k):
    """Deterministic conformers for the golden-energy test; ``k = 0`` is all-trans."""
    if k == 0:
        return np.full(n, 180.0)
    return -180.0 + ((np.arange(n) * (71 + 40 * k) + 29 * k) % 360).astype(float)


# Energies (total, LJ, Coulomb) of a default-constructed SimpleFF, in hex float so the
# comparison is bit-for-bit rather than "close".  Recorded from the code as it stood
# before the parameters became settable; every default in SimpleFF is still that
# potential, and this is what says so.
GOLDEN_ENERGIES = {
    ("pe", 0): ("0x1.1e3fc00e8da83p+2", "0x1.6da9b630c824dp-1", "0x1.e1151290e9472p+1"),
    ("pe", 1): ("0x1.a3b106a39b6b4p+11", "0x1.a2304c690a679p+11", "0x1.3f6baf81062a0p+1"),
    ("pe", 2): ("0x1.339abec54a8d0p+9", "0x1.2b05af4b379bap+9", "0x1.27a0f2f1adc5ap+1"),
    ("pvdf", 0): ("-0x1.3fde07fe80080p+3", "0x1.d9cf7d9295a83p+2", "-0x1.1662e363e56e1p+4"),
    ("pvdf", 1): ("0x1.ed8901cec023ep+15", "0x1.ed9ea3604bb38p+15", "-0x1.457c33b4de016p+4"),
    ("pvdf", 2): ("0x1.f424484648882p+10", "0x1.f74621458c1f0p+10", "-0x1.b62410b50e2b7p+4"),
    ("pvdc", 0): ("0x1.6cd76b45a239ep+7", "0x1.790398e958758p+7", "-0x1.8585b476c774fp+2"),
    ("pvdc", 1): ("0x1.0b7800064280cp+20", "0x1.0b780f7b4dc39p+20", "-0x1.dcce6531704c4p+2"),
    ("pvdc", 2): ("0x1.36594713d009cp+14", "0x1.365185143630cp+14", "-0x1.4a35315572436p+3"),
}
_PROBES = {"pe": (PE, 8), "pvdf": (PVDF, 8), "pvdc": (get_polymer("pvdc"), 6)}


@pytest.mark.parametrize("name,k", sorted(GOLDEN_ENERGIES))
def test_default_simpleff_is_bit_for_bit_the_illustrative_potential(name, k):
    """A default ``SimpleFF()`` must reproduce the pre-parameterisation energies exactly.

    Not ``approx``: the point of the defaults is that selecting a fitted preset is the
    *only* way the potential changes, so any drift at all -- a reordered sum, a promoted
    dtype -- is a regression, and the last bit is where such a change shows up first.
    """
    poly, n = _PROBES[name]
    s = build_chain(poly, _probe_dihedrals(n, k))
    e, lj, coulomb = (float.fromhex(v) for v in GOLDEN_ENERGIES[(name, k)])
    ff = SimpleFF()
    assert ff.energy(s) == e
    comp = ff.components(s)
    assert comp["lj"] == lj
    assert comp["coulomb"] == coulomb
    assert SimpleFF.from_preset("illustrative").energy(s) == e


def test_torsion_by_bond_defaults_to_the_shared_triple():
    """Per-bond-type coefficients, with the same triple everywhere, change nothing."""
    s = build_chain(PVDF, _probe_dihedrals(8, 2))
    shared = SimpleFF()
    per_bond = SimpleFF(torsion_by_bond=(shared.torsion, shared.torsion))
    assert per_bond.energy(s) == shared.energy(s)
    assert per_bond.torsion_coefficients(PVDF, 8).shape == (8, 3)
    # bond type of dihedral j is j % bonds_per_repeat, the convention fit_ris scans with
    V = SimpleFF(torsion_by_bond=((1.0, 0.0, 0.0), (0.0, 0.0, 2.0))).torsion_coefficients(PVDF, 4)
    assert np.array_equal(V, np.array([[1.0, 0, 0], [0, 0, 2.0], [1.0, 0, 0], [0, 0, 2.0]]))


def test_torsion_by_bond_separates_the_two_bond_types():
    """Each dihedral is driven by its own bond type's triple, summed by hand."""
    dih = _probe_dihedrals(8, 1)
    s = build_chain(PVDF, dih)
    a, b = (2.0, 0.3, 1.0), (-1.0, 0.7, 3.0)
    nb = SimpleFF(torsion=(0.0, 0.0, 0.0)).energy(s)  # the nonbonded part alone
    mixed = SimpleFF(torsion_by_bond=(a, b)).energy(s)

    def V(t, phi):
        return 0.5 * (t[0] * (1 + np.cos(phi)) + t[1] * (1 - np.cos(2 * phi)) + t[2] * (1 + np.cos(3 * phi)))

    phi = np.deg2rad(dih)
    expect = nb + sum(V(a if j % 2 == 0 else b, phi[j]) for j in range(len(dih)))
    assert mixed == pytest.approx(expect, abs=1e-9)
    assert mixed != SimpleFF(torsion=a).energy(s) and mixed != SimpleFF(torsion=b).energy(s)
    with pytest.raises(ValueError, match="must be 2 triples"):
        SimpleFF(torsion_by_bond=((1.0, 0.0, 0.0),)).energy(s)


def test_charge_scale_and_eps_r_are_degenerate_in_the_energy():
    """Coulomb goes as charge_scale^2 / eps_r -- the degeneracy the fit has to live with."""
    s = build_chain(PVDF, _probe_dihedrals(8, 2))
    base = SimpleFF().components(s)["coulomb"]
    scaled = SimpleFF(charge_scale=0.8).components(s)
    assert scaled["coulomb"] == pytest.approx(0.64 * base, rel=1e-12)
    assert SimpleFF(charge_scale=0.8).energy(s) == pytest.approx(SimpleFF(eps_r=1 / 0.64).energy(s), rel=1e-12)
    assert SimpleFF(eps_r=2.0).components(s)["coulomb"] == pytest.approx(0.5 * base, rel=1e-12)


def test_lj_overrides_are_per_instance_and_never_touch_the_global_table():
    from polyfind.polymers import UFF_LJ

    before = dict(UFF_LJ)
    s = build_chain(PVDF, _probe_dihedrals(8, 0))
    plain = SimpleFF()
    fat = SimpleFF(lj={"F": (3.9, 0.050)})
    assert fat.energy(s) != plain.energy(s)
    assert fat.lj_table()["F"] == (3.9, 0.050)
    assert fat.lj_table()["C"] == UFF_LJ["C"]
    assert UFF_LJ == before  # the module-level table is untouched
    # the table is UFF plus the main-group fallbacks a Frame of some other chemistry needs;
    # every element polymers.py defines keeps polymers.py's value
    assert {k: v for k, v in plain.lj_table().items() if k in UFF_LJ} == UFF_LJ
    assert plain.energy(s) == float.fromhex(GOLDEN_ENERGIES[("pvdf", 0)][0])  # still exact


def test_topology_cache_is_keyed_on_the_parameters():
    """Two calculators differing only in a parameter must not share a cached topology."""
    s = build_chain(PE, _probe_dihedrals(8, 1))
    a, b = SimpleFF(), SimpleFF(scale14=0.25, eps_r=3.0, lj={"H": (3.1, 0.044)})
    ea, eb = a.energy(s), b.energy(s)
    assert ea != eb
    assert a.energy(s) == ea and b.energy(s) == eb  # cached values are still each its own
    assert len(a._cache) == 1 and len(b._cache) == 1


def test_unknown_preset_names_are_rejected():
    with pytest.raises(KeyError, match="unknown SimpleFF preset"):
        SimpleFF.from_preset("no-such-preset")


# ------------------------------------------------- arbitrary geometries (Frame)
def _frame_from_chain(polymer, dihedrals):
    """A :class:`Frame` built from a chain, so the tests need no external data file."""
    s = build_chain(polymer, dihedrals)
    return Frame.from_geometry(s.elements, s.coords, system=polymer.name, source="synthetic"), s


# Torsion sets to build test frames from: the polymorph conformations plus a couple of
# mixed ones, jittered off their ideal angles.  Deliberately not uniform draws on the
# circle -- with *frozen* bond angles an arbitrary torsion set can fold the chain onto
# itself until two atoms four bonds apart are inside covalent bonding distance (the
# pentane G+G- clash does it on its own), at which point any distance-based bond
# perception, ours or anyone else's, invents a ring.  That is an artifact of the rigid
# geometry; the reference data holds relaxed conformers, which these imitate.
_CONFORMERS = [
    [178.0, -176.0, 179.0, 174.0, -178.0, 177.0, 180.0],      # all-trans (beta)
    [172.0, 62.0, 176.0, -58.0, 168.0, 65.0, -175.0],         # TGTG' (alpha)
    [175.0, 171.0, 178.0, 58.0, -174.0, 176.0, 172.0],        # T3GT3 (gamma-like)
    [-172.0, 178.0, 64.0, 170.0, 176.0, -61.0, 174.0],
    [166.0, -63.0, -168.0, -59.0, 173.0, 178.0, -177.0],
]


def test_inferred_bonds_and_elements_match_build_chain_for_pvdf():
    """The Frame path must see exactly the molecule build_chain builds.

    This is the check that licenses fitting to reference geometries at all: if the
    inferred topology of a five-monomer PVDF oligomer were not bond-for-bond the one the
    package builds for itself, the fit would be scoring a different molecule.
    """
    for dih in _CONFORMERS:
        frame, s = _frame_from_chain(PVDF, dih)
        assert sorted(frame.elements) == sorted(s.elements)
        assert len(frame.bonds) == len(s.bonds)
        assert len(frame.backbone) == len(s.backbone) == 10
        assert [frame.elements[i] for i in frame.backbone] == ["C"] * 10
        # reordering into build_chain's order reproduces the element list and the bond
        # graph exactly -- to_structure raises if either differs
        rebuilt = frame.to_structure(PVDF)
        assert rebuilt.elements == list(s.elements)
        assert (sorted(tuple(sorted(b)) for b in rebuilt.bonds)
                == sorted(tuple(sorted(b)) for b in s.bonds))


def test_frame_energy_equals_structure_energy():
    """The same potential and the same number, whichever path built the topology."""
    ff = SimpleFF(charge_increments=(("C", "H", -0.10), ("C", "F", 0.20)))
    for dih in _CONFORMERS:
        frame, _ = _frame_from_chain(PVDF, dih)
        assert ff.energy_frame(frame) == pytest.approx(ff.energy(frame.to_structure(PVDF)), abs=1e-9)


def test_bond_charge_increments_reproduce_the_polymer_tables():
    """The illustrative charges *are* a bond-charge-increment model; these are its numbers."""
    from polyfind.forcefield import bci_charges

    cases = ((PVDF, {("C", "H"): -0.10, ("C", "F"): 0.20}),
             (PE, {("C", "H"): -0.06}),
             (get_polymer("pvdc"), {("C", "H"): -0.10, ("C", "Cl"): 0.10}))
    for polymer, inc in cases:
        s = build_chain(polymer, np.full(7, 180.0))
        q = bci_charges(s.elements, s.bonds, inc)
        assert abs(q.sum()) < 1e-12  # neutral by construction, for any increments
        # Away from the ends the two agree exactly.  build_chain gives its two cap
        # hydrogens charge zero and leaves the carbons they hang off unbalanced, which a
        # neutral increment model cannot and should not reproduce.
        caps = {s.n_atoms - 1, s.n_atoms - 2}
        touching = caps | {a for a, b in s.bonds if b in caps} | {b for a, b in s.bonds if a in caps}
        interior = [i for i in range(s.n_atoms) if i not in touching]
        assert np.allclose(q[interior], s.charges[interior]), polymer.name


def test_rotor_changes_one_torsion_and_no_bond_or_angle():
    """The projection the force fit uses is the only displacement with this property."""
    from polyfind.forcefield import _side_mask, dihedral_angles, rotate_about_bond

    frame, _ = _frame_from_chain(PVDF, _CONFORMERS[1])
    a0 = frame.angles()
    lengths = lambda c: np.array([np.linalg.norm(c[i] - c[j]) for i, j in frame.bonds])
    l0 = lengths(frame.coords)
    wrap = lambda x: (x + 180.0) % 360.0 - 180.0
    for k, (_, b, c, _) in enumerate(np.asarray(frame.torsions, dtype=int)):
        mask = _side_mask(frame.n_atoms, frame.adjacency, int(b), int(c))
        h = 1e-4
        moved = rotate_about_bond(frame.coords, mask, int(b), int(c), h)
        d = wrap(dihedral_angles(moved, frame.torsions) - a0) / np.degrees(h)
        want = np.zeros(frame.n_dihedrals)
        want[k] = 1.0
        assert np.allclose(d, want, atol=1e-6)
        assert np.allclose(lengths(moved), l0, atol=1e-12)


def test_frame_torques_are_the_projection_of_the_cartesian_gradient():
    """``frame_torques`` must equal ``sum_i F_i . u_ki`` for the model's own forces."""
    ff = SimpleFF(charge_increments=(("C", "H", -0.10), ("C", "F", 0.20)))
    frame, _ = _frame_from_chain(PVDF, _CONFORMERS[3])
    h = 1e-5
    F = np.zeros((frame.n_atoms, 3))
    for i in range(frame.n_atoms):
        for a in range(3):
            cp = frame.coords.copy()
            cp[i, a] += h
            cm = frame.coords.copy()
            cm[i, a] -= h
            F[i, a] = -(ff.energy_frame(frame.with_coords(cp))
                        - ff.energy_frame(frame.with_coords(cm))) / (2 * h)
    assert np.allclose(np.einsum("kia,ia->k", frame.rotors(), F), ff.frame_torques(frame), atol=1e-6)


def test_a_bond_stretch_force_pair_projects_to_exactly_zero():
    """Why forces are fitted through the torsional projection and not in Cartesians.

    A potential with no bonded terms cannot produce a bond-stretch force, and a reference
    method evaluated at a geometry relaxed by some *other* method is dominated by exactly
    that (measured on the reference set: 93% of the sum of squared force components).  The
    projection annihilates it identically, so what survives is the part this potential is
    answerable for.
    """
    frame, _ = _frame_from_chain(PVDF, _CONFORMERS[2])
    rng = np.random.default_rng(7)
    F = np.zeros((frame.n_atoms, 3))
    for i, j in frame.bonds:  # an arbitrary force along every bond, equal and opposite
        u = frame.coords[j] - frame.coords[i]
        u = u / np.linalg.norm(u) * rng.normal()
        F[i] -= u
        F[j] += u
    assert np.abs(np.einsum("kia,ia->k", frame.rotors(), F)).max() < 1e-10


def test_read_frames_needs_a_path_or_the_environment_variable(monkeypatch):
    from polyfind.forcefield import TRAINSET_ENV, read_frames

    monkeypatch.delenv(TRAINSET_ENV, raising=False)
    with pytest.raises(FileNotFoundError, match=TRAINSET_ENV):
        read_frames()
