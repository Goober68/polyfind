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
    """An unsymmetrised fit, so the relations can be measured rather than assumed.

    Rigid angles throughout this section: it measures the fit's symmetry algebra (reversal
    images, bond-type shifts), which does not depend on how the backbone angles are
    treated, and the relaxed scan PVDC, CFE and CDFE now default to costs a hundred times
    more per fit.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the chiral polymers warn once on build_chain
        return fit_ris(get_polymer(name), SimpleFF(), step=step, n_monomers=3,
                       third_order=third_order, symmetrize=False, angles="rigid")


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
        model = fit_ris(get_polymer(name), SimpleFF(), step=20.0, n_monomers=3, third_order=True, angles="rigid").model
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
        model = fit_ris(get_polymer(name), SimpleFF(), step=20.0, n_monomers=3, angles="rigid").model
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


def test_adaptive_pvdf_g_plus_g_minus_infimum_sits_at_the_basin_edge():
    """PVDF's CF2-centred (bond type 1) G+/G- pair is the one place adaptive and dense
    disagree by more than grid/refinement noise, and it is expected, not a bug: the
    infimum of that basin sits at the basin's own *edge*, the G+/T boundary at 120 deg,
    where no grid of step 10 deg can sample it.

    Re-measured (the numbers below, and the two pinned values, moved with the PVDF C-C
    correction from 1.54 to 1.528 A). Profile of the pair energy over the G+ basin, at
    each phi's own best psi in the G- basin, on a 2 deg grid: it falls from 9.87 at
    phi = 106 monotonically to 4.97 at phi = 118, the last grid point inside the basin,
    and the bounded continuous search reaches **4.75 at phi = 120 - 1e-6, psi = -75**.
    There is a second, shallow interior minimum near phi = 72 at 4.97, which is *above*
    the edge infimum and so does not change the answer -- an earlier version of this
    docstring said the basin had no interior critical point at all, and that overstated
    it. The dense step = 10 deg scan can only report **8.77**, because its G+ samples
    stop at 110 deg (120 itself is handed to T by ``_basins``' tie-break; see
    ``_basin_bounds``). This is asserted explicitly, as expected behaviour, rather than
    folded into a widened tolerance in the general accuracy test above.
    """
    dense = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=6, scan="dense")
    adap = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=6, scan="adaptive")
    gp, gm = THREE_STATE.index("G+"), THREE_STATE.index("G-")

    d = dense.model.second_order[1, gp, gm]
    a = adap.model.second_order[1, gp, gm]
    assert d == pytest.approx(8.77, abs=0.2)
    assert a == pytest.approx(4.75, abs=0.2)
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
#
# The PVDF and PVDC rows were RE-RECORDED when the batched geometry corrections landed
# (PVDF C-C 1.54 -> 1.528 A, PVDC backbone angles 114/114 -> 123/114; see
# docs/REFERENCES.md).  Those change the *structures* the potential is evaluated on, not
# the potential, which is exactly why the PE rows below are untouched: PE's geometry did
# not move, so its three energies still match the pre-parameterisation recording bit for
# bit, and they are what still guards the potential itself against drift.
GOLDEN_ENERGIES = {
    ("pe", 0): ("0x1.1e3fc00e8da83p+2", "0x1.6da9b630c824dp-1", "0x1.e1151290e9472p+1"),
    ("pe", 1): ("0x1.a3b106a39b6b4p+11", "0x1.a2304c690a679p+11", "0x1.3f6baf81062a0p+1"),
    ("pe", 2): ("0x1.339abec54a8d0p+9", "0x1.2b05af4b379bap+9", "0x1.27a0f2f1adc5ap+1"),
    ("pvdf", 0): ("-0x1.0adb43a3bf2e9p+3", "0x1.1e9633e714343p+3", "-0x1.14b8bbc569b16p+4"),
    ("pvdf", 1): ("0x1.2759360e35fccp+16", "0x1.2763f688308d3p+16", "-0x1.4477470167d44p+4"),
    ("pvdf", 2): ("0x1.21906e59e0153p+11", "0x1.232546efdb088p+11", "-0x1.b81a1be1a2196p+4"),
    ("pvdc", 0): ("0x1.351d41daeb829p+5", "0x1.67f392c00d8cbp+5", "-0x1.96b287291050ep+2"),
    ("pvdc", 1): ("0x1.b8d2ad49ad013p+15", "0x1.b8d5965338f4ep+15", "-0x1.fc1b69a4a43e5p+2"),
    ("pvdc", 2): ("0x1.d64cf03c1fdddp+9", "0x1.d53d2eb589facp+9", "-0x1.4454cc7ec19a2p+3"),
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


# ------------------------------------------- valence terms and off-site charges
BOND_TERMS = (("C-C", 620.0, 1.526), ("C-H", 680.0, 1.090), ("C-F", 740.0, 1.380), ("*", 500.0, 1.6))
ANGLE_TERMS = (("C-C-C", 120.0, 112.7), ("C-C-H", 100.0, 109.5), ("C-C-F", 120.0, 109.0),
               ("H-C-H", 70.0, 107.8), ("F-C-F", 160.0, 107.0), ("F-C-H", 100.0, 107.0),
               ("*", 110.0, 109.5))
PVDF_INC = (("C", "H", -0.10), ("C", "F", 0.20))


def test_valence_terms_are_off_by_default():
    """The rigid potential is still the default, and nothing about it moved."""
    ff = SimpleFF()
    assert ff.bond_terms is None and ff.angle_terms is None and ff.charge_offsets is None
    assert ff.bond_table() == {} and ff.angle_table() == {} and ff.offset_table() == {}
    s = build_chain(PVDF, np.full(8, 180.0))
    assert ff.components(s)["bond"] == 0.0 and ff.components(s)["angle"] == 0.0
    # the golden test above already pins the energy itself, bit for bit


def test_valence_energy_is_an_additive_constant_on_a_rigid_chain():
    """Why the packing side cannot be disturbed by this, measured rather than argued.

    Every bond length and every bond angle of a ``build_chain`` chain is frozen, so a
    valence term is one number per topology however the torsions are set.  A lattice
    energy *difference*, an RIS energy (measured from all-trans) and a polarization
    therefore cannot move -- which is what lets the valence terms be added without
    re-measuring a single packing result.
    """
    plain = SimpleFF()
    valence = SimpleFF(bond_terms=BOND_TERMS, angle_terms=ANGLE_TERMS)
    dihedrals = _CONFORMERS + [np.full(7, 180.0), np.array([60.0, -60.0] * 3 + [180.0])]
    pairs = [(valence.energy(build_chain(PVDF, d)), plain.energy(build_chain(PVDF, d)))
             for d in dihedrals]
    offsets = [v - p for v, p in pairs]
    scale = max(abs(p) for _, p in pairs)  # the G+G- clash conformer is ~1e5 kcal/mol
    assert max(offsets) - min(offsets) < 1e-11 * scale  # i.e. rounding, not a torsion dependence
    comp = valence.components(build_chain(PVDF, _CONFORMERS[1]))
    assert comp["bond"] > 0 and comp["angle"] > 0
    assert comp["lj"] == plain.components(build_chain(PVDF, _CONFORMERS[1]))["lj"]


def test_valence_terms_are_the_harmonics_they_claim_to_be():
    """Hand-summed over the bond graph, against the potential's own number."""
    from polyfind.chain import angle as chain_angle, distance
    from polyfind.forcefield import valence_topology

    ff = SimpleFF(bond_terms=BOND_TERMS, angle_terms=ANGLE_TERMS)
    s = build_chain(PVDF, _CONFORMERS[2])
    bl, al = valence_topology(s.elements, s.bonds)
    bt, at = ff.bond_table(), ff.angle_table()
    e_b = sum(0.5 * bt.get(n, bt["*"])[0] * (distance(s.coords, i, j) - bt.get(n, bt["*"])[1]) ** 2
              for i, j, n in bl)
    e_a = sum(0.5 * at.get(n, at["*"])[0]
              * np.radians(chain_angle(s.coords, i, j, k) - at.get(n, at["*"])[1]) ** 2
              for i, j, k, n in al)
    comp = ff.components(s)
    assert comp["bond"] == pytest.approx(e_b, abs=1e-9)
    assert comp["angle"] == pytest.approx(e_a, abs=1e-9)


def test_bond_and_angle_types_come_from_the_graph():
    """A nitrile's C-N is not a single C-N, and a linear centre gets no bend term."""
    from polyfind.forcefield import valence_topology

    s = build_chain(get_polymer("vdcn"), np.full(5, 180.0))
    bl, al = valence_topology(s.elements, s.bonds)
    assert "C=N" in {n for _, _, n in bl}  # the nitrile, typed apart from a single C-N
    assert "C-N" not in {n for _, _, n in bl}
    # the nitrile carbon has two neighbours, so its 180-degree angle carries no term
    nitrile_c = [i for i, e in enumerate(s.elements)
                 if e == "C" and sorted(s.elements[j] for a, b in s.bonds
                                        for j in (b,) if a == i) == ["N"]]
    centres = {j for _, j, _, _ in al}
    assert nitrile_c and not (set(nitrile_c) & centres)
    pvdf = build_chain(PVDF, np.full(7, 180.0))
    assert {n for _, _, n in valence_topology(pvdf.elements, pvdf.bonds)[0]} == {"C-C", "C-H", "C-F"}
    # F-C-H is there because build_chain caps the last CF2 with a hydrogen
    assert {n for _, _, _, n in valence_topology(pvdf.elements, pvdf.bonds)[1]} == {
        "C-C-C", "C-C-H", "C-C-F", "H-C-H", "F-C-F", "F-C-H"}


def test_an_off_site_charge_moves_the_coulomb_term_and_nothing_else():
    """The point of the offset: fluorine's charge and its van der Waals centre come apart."""
    s = build_chain(PVDF, _CONFORMERS[1])
    plain = SimpleFF(charge_increments=PVDF_INC)
    zero = SimpleFF(charge_increments=PVDF_INC, charge_offsets=(("F", 0.0),))
    moved = SimpleFF(charge_increments=PVDF_INC, charge_offsets=(("F", 0.18),))
    assert zero.energy(s) == plain.energy(s)  # d = 0 is exactly the atom-centred model
    assert moved.components(s)["lj"] == plain.components(s)["lj"]  # LJ stays on the nucleus
    assert abs(moved.components(s)["coulomb"] - plain.components(s)["coulomb"]) > 0.5


def test_only_terminal_atoms_may_carry_an_off_site_charge():
    from polyfind.forcefield import offset_sites

    s = build_chain(PVDF, np.full(7, 180.0))
    sites = offset_sites(s.elements, s.bonds, {"F": 0.2})
    assert len(sites) == 10 and all(s.elements[a] == "F" for a, _, _ in sites)
    with pytest.raises(ValueError, match="only terminal atoms"):
        offset_sites(s.elements, s.bonds, {"C": 0.2})


def test_the_lattice_projection_keeps_the_charge_and_the_dipole():
    """How an off-site model reaches ``pack``, which has one charge per atom and no sites."""
    from polyfind.forcefield import bci_charges, charge_site_coords, offset_sites, project_offset_charges

    s = build_chain(PVDF, _CONFORMERS[3])
    q = bci_charges(s.elements, s.bonds, {("C", "H"): -0.10, ("C", "F"): 0.20})
    sites = offset_sites(s.elements, s.bonds, {"F": 0.18})
    q_proj = project_offset_charges(s.coords, q, sites)
    mu_off = (q[:, None] * charge_site_coords(s.coords, sites)).sum(0)
    mu_proj = (q_proj[:, None] * s.coords).sum(0)
    assert q_proj.sum() == pytest.approx(q.sum(), abs=1e-12)
    assert np.allclose(mu_proj, mu_off, atol=1e-12)  # same dipole, exactly
    assert not np.allclose(q_proj, q)  # and it is a real change of the charges


def test_forces_frame_is_the_gradient_the_torques_are_projected_from():
    """The Cartesian forces the fit now uses, checked against what was already tested."""
    ff = SimpleFF(charge_increments=PVDF_INC, bond_terms=BOND_TERMS, angle_terms=ANGLE_TERMS,
                  charge_offsets=(("F", 0.15),))
    frame, _ = _frame_from_chain(PVDF, _CONFORMERS[1])
    F = ff.forces_frame(frame)
    assert np.abs(F.sum(axis=0)).max() < 1e-5  # no net force: the energy is translation-invariant
    assert np.allclose(np.einsum("kia,ia->k", frame.rotors(), F), ff.frame_torques(frame), atol=1e-4)


def test_relaxing_the_backbone_angles_needs_a_bend_term_and_then_relieves_strain():
    """The strain finding, and what having an angle term does to it.

    PVDC's all-trans chain is the case ``polymers.py`` records: 313 kcal/mol of
    Lennard-Jones strain when both backbone angles were frozen at 114 deg, 74 now that
    the CH2 angle carries its measured 123 deg.  Given a bend term the angles can open
    further, which is what a real chain does; the test asserts the direction and the
    mechanism, not a fitted number (those are in docs/VALENCE_FIT.md).

    The threshold moved with the angle correction, and that is the point rather than an
    inconvenience: with this potential the relaxation used to buy 249 kcal/mol from a
    114/114 start and now buys 28 from a 123/114 one, because most of what it used to
    recover was the error in the starting angles.  It still opens the same angle in the
    same direction and still lands on the same relaxed geometry -- 125.8 deg at CH2 and
    108.3 at CCl2, from either start.
    """
    from polyfind.forcefield import relax_backbone_angles

    pvdc = get_polymer("pvdc")
    ff = SimpleFF(bond_terms=BOND_TERMS, angle_terms=ANGLE_TERMS + (("Cl-C-Cl", 120.0, 109.5),
                                                                   ("C-C-Cl", 120.0, 109.0)))
    r = relax_backbone_angles(pvdc, ff, np.full(10, 180.0))
    assert r["energy"] < r["energy0"] - 20.0  # the strain was real and relaxing removes much of it
    assert r["lj"] < r["lj0"]
    assert r["angles"][0] > r["angles0"][0]  # the CH2 angle is still the one that opens
    assert r["angles"][1] < r["angles0"][1] - 3.0  # the CCl2 angle closes to compensate


# --------------------------------------------------------------------- charge flux
#
# The four things a geometry-dependent charge has to get right before anything downstream
# of it is worth reading: it reduces to the fixed model at zero coefficient, its analytic
# derivative is the derivative of its value, it stays neutral (or the dipole it feeds is not
# a dipole at all), and it refuses the cases where it is ill-defined.

def _flux_ff(**kw):
    """``pvdf-dft-valence`` with the given extra keyword arguments."""
    import polyfind.fitting  # noqa: F401 - registers the preset

    base = SimpleFF.from_preset("pvdf-dft-valence")
    return SimpleFF(**{**{k: v for k, v in base.__dict__.items() if k != "_cache"}, **kw})


def _pvdf_geometry():
    from polyfind.forcefield import infer_bonds

    st = build_chain(PVDF, [180.0] * 6)
    el, X = list(st.elements), np.asarray(st.coords, dtype=float)
    return el, X, infer_bonds(el, X)


def test_zero_charge_flux_is_exactly_the_fixed_increment_model():
    from polyfind.forcefield import bci_charges, offset_sites, project_offset_charges

    el, X, bonds = _pvdf_geometry()
    base = _flux_ff()
    q = bci_charges(el, bonds, base.increments())
    q = project_offset_charges(X, q, offset_sites(el, bonds, base.offset_table())) * base.charge_scale
    ff = _flux_ff(charge_flux=(("C", "H", 0.0, 0.0), ("C", "F", 0.0, 0.0)))
    assert np.abs(ff.flux_topology(el, bonds).charges(X)[0] - q).max() == 0.0
    assert not ff.has_flux()  # all-zero coefficients are not a flux


def test_charge_flux_derivatives_match_a_finite_difference():
    el, X, bonds = _pvdf_geometry()
    ff = _flux_ff(charge_flux=(("C", "H", -0.7, 0.31), ("C", "F", 0.45, -0.22)))
    top = ff.flux_topology(el, bonds)
    q, dq, _ = top.charges_and_grad(X, 0.0)
    assert np.abs(q - ff.charges_at(el, bonds, X)).max() == 0.0
    h = 1e-6
    num = np.empty_like(dq)
    for b in range(len(el)):
        for d in range(3):
            Xp, Xm = X.copy(), X.copy()
            Xp[b, d] += h
            Xm[b, d] -= h
            num[:, b, d] = (top.charges(Xp)[0] - top.charges(Xm)[0]) / (2 * h)
    assert np.abs(dq - num).max() < 1e-7
    assert abs(q.sum()) < 1e-12  # neutral whatever the flux does
    assert np.abs(dq.sum(axis=0)).max() < 1e-12  # and neutral for every displacement


def test_charge_flux_refuses_what_it_cannot_orient_or_evaluate():
    el, X, bonds = _pvdf_geometry()
    with pytest.raises(ValueError, match="homonuclear"):
        _flux_ff(charge_flux=(("C", "C", 0.3, 0.0),)).flux_topology(el, bonds)
    # a bond channel needs a reference length, which comes from the stretch terms
    with pytest.raises(ValueError, match="reference length"):
        SimpleFF(charge_increments=(("C", "H", -0.1),),
                 charge_flux=(("C", "H", 0.0, 0.5),)).flux_topology(el, bonds)
    # and the flux needs increments to perturb at all
    with pytest.raises(ValueError, match="charge_increments"):
        SimpleFF(charge_flux=(("C", "H", 0.3, 0.0),)).flux_topology(el, bonds)


def test_the_oligomer_energy_path_refuses_a_fluxing_potential():
    """Loud rather than silent: one fixed charge per topology cannot carry a flux.

    Using the base increments there would make the energy disagree with the dipole, which
    is exactly the failure mode this whole change exists to avoid.
    """
    ff = _flux_ff(charge_flux=(("C", "H", -0.7, 0.0),))
    with pytest.raises(ValueError, match="charge flux"):
        ff.energy(build_chain(PVDF, [180.0] * 6))


# ------------------------------------------------------------ backbone-angle relaxation
#
# docs/NITRILE_LANDSCAPE.md: with the backbone angles frozen, VDCN's one-bond profile has
# wells at +/-120 deg, 8.7 kcal/mol below trans, that a chain free to bend does not have;
# _basins hands them to the trans label, so the rigid fit gave VDCN a T state at 180 deg
# carrying the energy of a well at 120.  The relaxed scan removes that, the rigid scan is
# kept bit-for-bit, and which one a chemistry gets by default is measured, not assigned.

def _val():
    import polyfind.fitting  # noqa: F401 - registers the preset

    return SimpleFF.from_preset("pvdf-dft-valence")


def test_rigid_alias_and_explicit_rigid_are_bit_for_bit_the_default_fit():
    """PVDF is recorded rigid, so the default, ``scan="rigid"`` and ``angles="rigid"`` are
    the same fit to the last bit -- every PVDF number in every table depends on it."""
    ff = SimpleFF()
    a = fit_ris(PVDF, ff, step=20.0, n_monomers=5)
    b = fit_ris(PVDF, ff, step=20.0, n_monomers=5, scan="rigid")
    c = fit_ris(PVDF, ff, step=20.0, n_monomers=5, angles="rigid")
    assert a.angles == b.angles == c.angles == "rigid"
    assert a.discriminator is None and a.n_relaxation_energies == 0
    for x in (b, c):
        assert np.array_equal(a.model.first_order, x.model.first_order)
        assert np.array_equal(a.model.second_order, x.model.second_order)
        assert a.model.states.angles == x.model.states.angles
        assert a.n_evaluations == x.n_evaluations
        for k in a.scan1:
            assert np.array_equal(a.scan1[k], x.scan1[k]) and np.array_equal(a.scan2[k], x.scan2[k])
    with pytest.raises(ValueError, match="angles must be"):
        fit_ris(PVDF, ff, step=45.0, n_monomers=3, angles="bogus")


def test_relax_angles_batch_matches_per_row_scipy_and_only_lowers_the_energy():
    """The batched projected BFGS finds the same minima as an independent per-row L-BFGS-B
    on the same objective, and a relaxed conformer never costs more than the rigid one."""
    from scipy.optimize import minimize

    from polyfind.forcefield import relax_angles_batch

    ff = _val()
    rng = np.random.default_rng(3)
    dihs = rng.choice([180.0, 60.0, -60.0], size=(4, 8)) + rng.uniform(-10, 10, size=(4, 8))
    tpl, coords = build_chain_batch(PVDF, dihs)
    rigid = ff.energy_batch((tpl, coords))
    E, ang, status, n = relax_angles_batch(PVDF, ff, dihs)
    assert E.shape == (4,) and ang.shape == (4, 11) and n > 0
    assert np.all(E <= rigid + 1e-9)
    assert np.all(status <= 1)  # converged or stalled at the finite-difference floor, never out of iterations
    assert np.all((ang >= 95.0) & (ang <= 135.0))
    nominal = np.array([PVDF.backbone[k % 2].backbone_angle for k in range(11)])
    for m in range(4):
        def f(a):
            t, c = build_chain_batch(PVDF, dihs[m][None], bond_angles=a[None])
            return float(ff.energy_batch((t, c))[0])
        res = minimize(f, nominal, method="L-BFGS-B", bounds=[(95.0, 135.0)] * 11,
                       options={"maxiter": 500, "ftol": 1e-12})
        assert E[m] == pytest.approx(res.fun, abs=2e-3)


def test_local_minima_of_a_cyclic_profile():
    from polyfind.forcefield import _local_minima

    E = np.array([0.0, 1.0, 2.0, 1.0, 0.5, 1.5, 3.0, 2.0])  # minima at 0 (cyclic) and 4
    found = _local_minima(E)
    assert [i for i, _ in found] == [0, 4]
    # the barrier tops are the maxima at index 2 (2.0) and index 6 (3.0); index 3 is on a slope
    assert dict(found)[0] == pytest.approx(2.0)  # min(2.0, 3.0) - 0.0
    assert dict(found)[4] == pytest.approx(1.5)  # min(2.0, 3.0) - 0.5
    assert [i for i, _ in _local_minima(E, min_depth=1.8)] == [0]
    assert _local_minima(np.zeros(5)) == []


def test_the_discriminator_orphans_vdcn_s_frozen_angle_well_and_clears_pvdf():
    """The measurement behind ANGLE_RELAXATION_DEFAULTS, on the two rows that anchor it:
    VDCN's rigid +/-120 wells have no relaxed minimum within 30 deg, PVDF's rigid double
    gauche (+/-70, +/-40) sits within 30 of the relaxed +/-50 on both sides."""
    from polyfind.forcefield import ANGLE_RELAXATION_DEFAULTS, angle_relaxation_check

    ff = _val()
    vdcn = angle_relaxation_check(get_polymer("vdcn"), ff, step=10.0, n_monomers=6)
    assert vdcn["needs_relaxation"] is True and ANGLE_RELAXATION_DEFAULTS["vdcn"] is True
    for b in (0, 1):
        orphans = sorted(a for a, _, _ in vdcn["bonds"][b]["orphans"])
        assert orphans == [-120.0, 120.0]
        assert all(e < -8.0 for _, e, _ in vdcn["bonds"][b]["orphans"])  # the well is deep, not grid noise
        relaxed = sorted(a for a, _, _ in vdcn["bonds"][b]["relaxed_minima"])
        assert relaxed == [-180.0, -40.0, 40.0]
    pvdf = angle_relaxation_check(PVDF, ff, step=10.0, n_monomers=6)
    assert pvdf["needs_relaxation"] is False and ANGLE_RELAXATION_DEFAULTS["pvdf"] is False
    assert all(not pvdf["bonds"][b]["orphans"] for b in (0, 1))
    assert pvdf["restraint_k"] == 0.0  # the preset carries bend terms, so no fallback restraint


def test_auto_measures_a_polymer_it_has_no_record_for_and_falls_back_to_a_restraint():
    """A chemistry not in the table gets the check run on the fly (kept in the report), and
    a calculator without bend terms relaxes against the stated harmonic fallback."""
    from dataclasses import replace

    from polyfind.forcefield import FALLBACK_ANGLE_K, _angle_restraint_for

    unknown = replace(PE, name="pe-unrecorded")
    rep = fit_ris(unknown, SimpleFF(), step=30.0, n_monomers=3)
    assert rep.discriminator is not None and rep.discriminator["needs_relaxation"] is False
    assert rep.angles == "rigid"
    assert rep.discriminator["restraint_k"] == FALLBACK_ANGLE_K
    assert _angle_restraint_for(SimpleFF(), None) == FALLBACK_ANGLE_K
    assert _angle_restraint_for(_val(), None) == 0.0
    assert _angle_restraint_for(SimpleFF(), 7.5) == 7.5


def test_relaxed_fit_removes_vdcn_s_trans_artefact_and_the_third_order_guard_still_holds():
    """Rigid, VDCN's T state carries the energy of the +/-120 well (about -8.7 kcal/mol at
    the screen's step; here on a coarse grid that still samples 120).  Relaxed, T is a
    fraction of a kcal/mol from zero, the gauche state sits near +/-40 rather than +/-30,
    and no third-order term is a spurious inclusion-exclusion across an overlap: every
    unguarded triple is well inside the cap, and the guard fires only where the triple
    itself still overlaps after relaxing."""
    ff = _val()
    vdcn = get_polymer("vdcn")
    rigid = fit_ris(vdcn, ff, step=30.0, n_monomers=3, third_order=True, angles="rigid")
    relaxed = fit_ris(vdcn, ff, step=30.0, n_monomers=3, third_order=True, angles="relaxed")
    T = THREE_STATE.index("T")
    assert rigid.angles == "rigid" and relaxed.angles == "relaxed"
    assert np.all(rigid.model.first_order[:, T] < -5.0)  # the artefact
    assert np.all(relaxed.model.first_order[:, T] > -1.5)  # gone
    assert 25.0 <= abs(relaxed.model.states.angles[1]) <= 60.0  # on a 30 deg grid the +/-40 well reads as 30
    assert relaxed.relaxed_reference_angles.shape == (13,)
    assert np.all((relaxed.relaxed_reference_angles >= 95.0) & (relaxed.relaxed_reference_angles <= 135.0))
    assert relaxed.n_relaxation_energies > relaxed.n_evaluations
    st = relaxed.relaxation_status
    assert st["out_of_iterations"] == 0 and st["converged"] > 0
    # the third-order guard, relaxed: unguarded terms are moderate, guarded ones are +cap or 0
    g = relaxed.third_order_guard
    e3 = relaxed.model.third_order
    free = ~(g["capped"] | g["zeroed"])
    assert free.any()
    assert np.all(np.abs(e3[free]) < 25.0)
    assert np.all(e3[g["capped"]] == 50.0) and np.all(e3[g["zeroed"]] == 0.0)
    # and it never lets a difference between overlaps through as a bonus
    assert np.all(e3 > -25.0)
    # relaxing relieves contacts, so fewer triples overlap than under the rigid scan
    assert g["capped"].sum() < rigid.third_order_guard["capped"].sum()
