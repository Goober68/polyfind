"""A potential tilted to another's gradient at one geometry (examples/mlip_train.py relax --dft-gradient)."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ase = pytest.importorskip("ase")
SPEC = importlib.util.spec_from_file_location("mlip_train", Path(__file__).resolve().parents[1] / "examples" / "mlip_train.py")
MT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MT)


def _start(shake=0.05, sigma=2.3, epsilon=0.03):
    """A strained copper cell with a vacancy, labelled by one potential (the "DFT"), to be relaxed in another."""
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.calculators.lj import LennardJones

    a = bulk("Cu", "fcc", a=3.7, cubic=True).repeat((2, 1, 1))
    del a[0]  # a vacancy: the atoms around it relax, by an amount that is each potential's own
    a.positions += np.random.default_rng(3).normal(0.0, shake, a.positions.shape) if shake else 0.0
    ref = a.copy()
    ref.calc = EMT()
    a.arrays["REF_forces"] = ref.get_forces()
    a.info["REF_stress"] = ref.get_stress(voigt=False)
    a.info["REF_energy"] = ref.get_potential_energy()
    # the potential to relax in: the labelling one with a pair term added, so that its slope is off and its curvature nearly right
    from ase.calculators.mixing import SumCalculator

    return a, SumCalculator([EMT(), LennardJones(sigma=sigma, epsilon=epsilon, rc=6.0, smooth=True)]), EMT()


def test_at_the_start_the_tilted_potential_has_the_labels_forces_and_stress():
    start, base, _ = _start()
    a = start.copy()
    a.calc = MT.tilted_calculator(base, start)
    assert a.get_forces() == pytest.approx(start.arrays["REF_forces"], abs=1e-10)
    assert a.get_stress(voigt=False) == pytest.approx(start.info["REF_stress"], abs=1e-10)
    plain = start.copy()
    plain.calc = base
    assert a.get_potential_energy() == pytest.approx(plain.get_potential_energy(), abs=1e-10)
    assert np.abs(plain.get_forces() - start.arrays["REF_forces"]).max() > 0.01  # the two potentials do differ


def test_its_forces_and_stress_are_its_energys_own_derivatives_away_from_the_start():
    start, base, _ = _start()
    a = start.copy()
    a.set_cell(start.cell.array @ (np.eye(3) + np.diag([0.004, -0.003, 0.006])), scale_atoms=True)
    a.positions += np.random.default_rng(5).normal(0.0, 0.02, a.positions.shape)
    a.calc = MT.tilted_calculator(base, start)
    f, sig, h = a.get_forces(), a.get_stress(voigt=False), 1e-5
    for i, k in ((0, 0), (2, 1), (3, 2)):
        e = []
        for sgn in (+1, -1):
            b = a.copy()
            b.positions[i, k] += sgn * h
            b.calc = MT.tilted_calculator(base, start)
            e.append(b.get_potential_energy())
        assert -(e[0] - e[1]) / (2 * h) == pytest.approx(f[i, k], abs=1e-6)
    e = []
    for sgn in (+1, -1):
        b = a.copy()
        b.set_cell(a.cell.array @ (np.eye(3) + np.diag([0.0, 0.0, sgn * h])), scale_atoms=True)
        b.calc = MT.tilted_calculator(base, start)
        e.append(b.get_potential_energy())
    assert (e[0] - e[1]) / (2 * h) / a.get_volume() == pytest.approx(sig[2, 2], abs=1e-6)


def test_one_relaxation_in_it_lands_closer_to_the_labels_minimum_than_the_plain_potential_does():
    """The Newton step the tilt amounts to, where it is meant for: the potential's slope is off (a
    pair term pushes the vacancy's neighbours in) and the start is near the labels' minimum.  The
    labelling potential's forces after one relaxation are then a fraction of what relaxing in the
    plain potential leaves; what is left is the curvature error times the step taken."""
    from ase.optimize import BFGS

    start, base, label = _start(shake=0.0, sigma=2.5, epsilon=0.05)
    left = {}
    for name, calc in (("plain", base), ("tilted", MT.tilted_calculator(base, start))):
        a = start.copy()
        a.calc = calc
        BFGS(a, logfile=None).run(fmax=1e-4, steps=400)
        b = a.copy()
        b.calc = label
        left[name] = float(np.abs(b.get_forces()).max())
    assert left["plain"] > 0.05 and left["tilted"] < 0.25 * left["plain"]
