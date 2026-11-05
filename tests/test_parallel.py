"""Independent relaxations in parallel processes give the serial results."""
import os

import numpy as np

from polyfind import electrostriction as ES
from polyfind import mechanics as M
from polyfind.polymers import PVDF

T = 0


def _beta():
    ref, _ = M.refined_reference(PVDF, [T, T], label="beta", pack_kw={"table_cache_dir": None}, lj_cutoff="force")
    return M.relax_reference(ref)


def test_parallel_map_runs_in_workers_and_falls_back_serially():
    assert M.parallel_map(abs, [-1, -2, -3], workers=3) == [1, 2, 3]
    assert M.parallel_map(abs, [-1, -2, -3], workers=1) == [1, 2, 3]
    assert M.parallel_map(abs, [], workers=4) == []


def test_strain_states_and_susceptibilities_agree_between_serial_and_parallel(monkeypatch):
    ref = _beta()
    monkeypatch.setenv("POLYFIND_WORKERS", "1")
    serial_el = M.elastic_constants(ref)
    serial_di = M.dielectric_tensor(ref)
    serial_chi = ES.nominal_susceptibility(ref, np.zeros(6))
    serial_eps = ES.true_dielectric(ref, np.array([0.01, 0, 0, 0, 0, 0]))
    monkeypatch.setenv("POLYFIND_WORKERS", "4")
    par_el = M.elastic_constants(ref)
    par_di = M.dielectric_tensor(ref)
    par_chi = ES.nominal_susceptibility(ref, np.zeros(6))
    par_eps = ES.true_dielectric(ref, np.array([0.01, 0, 0, 0, 0, 0]))
    np.testing.assert_allclose(par_el.block, serial_el.block, rtol=0, atol=1e-9)
    np.testing.assert_allclose(par_di.relaxed, serial_di.relaxed, rtol=0, atol=1e-10)
    np.testing.assert_allclose(par_chi, serial_chi, rtol=0, atol=1e-10)
    np.testing.assert_allclose(par_eps, serial_eps, rtol=0, atol=1e-10)
    monkeypatch.setenv("POLYFIND_WORKERS", "1")
    serial_pz = M.piezoelectric(ref, serial_el)
    monkeypatch.setenv("POLYFIND_WORKERS", "4")
    par_pz = M.piezoelectric(ref, serial_el)
    np.testing.assert_allclose(par_pz.e, serial_pz.e, rtol=0, atol=1e-12)
    np.testing.assert_allclose(par_pz.d_direct, serial_pz.d_direct, rtol=0, atol=1e-9)
    assert par_pz.converse_iterations == serial_pz.converse_iterations


def test_electrostriction_agrees_between_serial_and_parallel(monkeypatch):
    ref = _beta()
    el = M.elastic_constants(ref)
    monkeypatch.setenv("POLYFIND_WORKERS", "1")
    serial = ES.electrostriction(ref, el)
    monkeypatch.setenv("POLYFIND_WORKERS", "4")
    par = ES.electrostriction(ref, el)
    np.testing.assert_allclose(par.M_strain, serial.M_strain, rtol=0, atol=1e-30, equal_nan=True)
    np.testing.assert_allclose(par.M_polarization, serial.M_polarization, rtol=0, atol=1e-30, equal_nan=True)
    np.testing.assert_allclose(par.chi_free, serial.chi_free, rtol=0, atol=1e-12)
