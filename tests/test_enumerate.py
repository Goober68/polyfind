import warnings

import numpy as np

from polyfind.enumerate import enumerate_periodic, table
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.helix import canonical_sequence
from polyfind.polymers import PE, PVDF, get_polymer
from polyfind.ris import polyethylene_like_model, RISModel
from polyfind.polymers import THREE_STATE


def _pvdf_synthetic_model():
    """The mirror-symmetric random model the PVDF enumeration tests share."""
    rng = np.random.default_rng(0)
    e1 = np.abs(rng.normal(size=(2, 3))) * 0.3
    e1[:, 0] = 0
    e2 = rng.normal(size=(2, 3, 3)) * 0.5
    e2[:, 0, 0] = 0
    mir = np.array(THREE_STATE.mirror)
    e1 = 0.5 * (e1 + e1[:, mir])
    e2 = 0.5 * (e2 + e2[:, mir][:, :, mir])
    return RISModel(THREE_STATE, 2, e1, e2)


def test_pe_textbook_model_ranks_all_trans_first():
    m = polyethylene_like_model()
    cands = enumerate_periodic(PE, m, max_period=6)
    assert cands[0].name == "T" and cands[0].known_as is not None
    names = [c.name for c in cands]
    assert len(names) == len(set(names))
    # the TG helix should be there and cost 0.5 kcal/mol per monomer (period 2 -> 0.25 per bond... per monomer = per bond for PE)
    tg = next(c for c in cands if c.name in ("TG+", "TG-"))
    assert abs(tg.energy_per_monomer - 0.25) < 1e-9
    assert tg.helix.periods_per_repeat == 3


def test_pvdf_enumeration_contains_polymorph_chains_once():
    m = _pvdf_synthetic_model()
    cands = enumerate_periodic(PVDF, m, max_period=8, k_per_period=80)
    known = [c.known_as for c in cands if c.known_as]
    assert "beta (TTTT)" in known and "alpha/delta (TGTG')" in known and "gamma/epsilon (T3GT3G')" in known
    assert len(known) == len(set(known))
    # all candidates propagate and are commensurate
    assert all(c.helix.rise_per_bond >= 0.3 for c in cands)
    # energies are per monomer and sorted
    ems = [c.energy_per_monomer for c in cands]
    assert ems == sorted(ems)
    assert "beta" in table(cands)


def test_an_achiral_polymer_still_collapses_mirror_images():
    """PVDF's candidate list must be exactly what it was before chirality was threaded in.

    Every candidate is the lexicographic representative of its class under the *full*
    group (shift, reversal, mirror and their compositions), and no two candidates are in
    the same class -- which is the statement that the achiral path is untouched.  The
    count is pinned as a regression on the list as a whole.
    """
    m = _pvdf_synthetic_model()
    cands = enumerate_periodic(PVDF, m, max_period=8, k_per_period=80)
    assert not PVDF.is_chiral
    assert len(cands) == 90
    keys = [canonical_sequence(c.seq, THREE_STATE, 2) for c in cands]  # achiral group
    assert keys == [c.seq for c in cands]
    assert len(set(keys)) == len(cands)
    # in particular no candidate's mirror image is also in the list, and the mirror
    # image of any of them has exactly the same energy: mirroring is a real symmetry here
    mir = np.array(THREE_STATE.mirror)
    worst = 0.0
    for c in cands:
        s = np.array(c.seq)
        mirrored = tuple(int(x) for x in mir[s])
        assert canonical_sequence(mirrored, THREE_STATE, 2) == c.seq  # same candidate
        e = m.energy(s, periodic=True) / (len(s) / m.B)
        em = m.energy(mirrored, periodic=True) / (len(s) / m.B)
        worst = max(worst, abs(e - em))
    assert worst == 0.0  # exactly, for an achiral polymer


def test_a_chiral_polymer_gains_the_mirror_partners_and_they_differ_in_energy():
    """CFE's G+ and G- helices are different molecules, so both must be enumerated.

    The old rule collapsed each pair into one candidate; with the chiral group the
    partners come back, and their RIS energies genuinely differ -- which is exactly the
    asymmetry that makes an isotactic chain pick one handedness.
    """
    cfe = get_polymer("cfe")
    assert cfe.is_chiral
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # angles="rigid": the partner count and energy floor below were calibrated on the
        # rigid scan; CFE's relaxed default keeps the partners distinct but with a smaller
        # G+/G- gap (3 rather than 24 kcal/mol), which puts one pair under the floor
        model = fit_ris(cfe, SimpleFF(), step=45.0, n_monomers=3, angles="rigid").model
    cands = enumerate_periodic(cfe, model, max_period=4, k_per_period=20)
    names = {c.name for c in cands}
    mir = np.array(model.states.mirror)

    partners, worst, worst_pair = 0, 0.0, None
    for c in cands:
        s = np.array(c.seq)
        partner = "".join(model.states.names[x] for x in mir[s])
        partners += partner in names
        e = model.energy(s, periodic=True) / (len(s) / model.B)
        em = model.energy(mir[s], periodic=True) / (len(s) / model.B)
        if abs(e - em) > worst:
            worst, worst_pair = abs(e - em), (c.name, partner)
    # the old (achiral) rule kept one representative per pair; the pairs are back
    old_classes = {canonical_sequence(c.seq, THREE_STATE, 2) for c in cands}
    assert len(cands) > len(old_classes)
    assert partners >= 10, partners
    assert {"TG+", "TG-"} <= names
    # ... and they are not degenerate: mirror partners differ by tens of kcal/mol
    assert worst > 5.0, (worst, worst_pair)
