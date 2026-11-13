"""The Ewald sum and the induced-dipole solve on a device array module, tested without a device.

``device_double`` is NumPy underneath and as strict as CuPy on top: an array of it is not a
NumPy array, will not become one by itself and will not mix with one.  A routine that runs on
it keeps host and device apart everywhere, and gives the host's numbers wherever it does the
host's arithmetic.  What only a real device can show -- that CuPy's own functions do what
their documentation says, and how long they take -- is ``tests/test_gpu.py`` and
``examples/ewald_device_benchmark.py``.
"""
import pickle

import device_double as dd
import numpy as np
import pytest

from polyfind import backend as bk
from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.ewald import Ewald, EwaldSpec, Thole
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.linegroup import repeat_chains
from polyfind.polarizability import Polarizable
from polyfind.polymers import THREE_STATE, pvdf_belt

T, GP, GM = 0, 1, 2
BELT4 = "pvdf-dft-valence-belt4"
EWALD = dict(lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
ARRAYS = ("grad_coords", "grad_lattice", "grad_charges", "grad_dipoles", "field", "tensor")
ENERGIES = ("real", "recip", "self_", "background", "surface")


def _cell(n=40, seed=3):
    """A skewed cell with atoms outside it, so several images and both shifts matter."""
    rng = np.random.default_rng(seed)
    h = np.array([[5.1, 0.3, 0.0], [-1.2, 8.3, 0.4], [0.2, -0.6, 6.7]])
    X = rng.uniform(-0.7, 1.6, (n, 3)) @ h
    q = rng.standard_normal(n)
    q -= q.mean()
    p = 0.1 * rng.standard_normal((n, 3))
    scale = rng.uniform(1.0, 1.6, (n, n))
    return X, q, h, p, 0.5 * (scale + scale.T)


def _wide_cell(n=90, seed=2):
    """A cell wide enough that its phases are worth keeping (sites times k vectors over a million)."""
    rng = np.random.default_rng(seed)
    h = np.diag([18.0, 19.0, 20.0])
    X = rng.uniform(0, 1, (n, 3)) @ h
    q = rng.normal(size=n)
    return X, q - q.mean(), h, 0.1 * rng.standard_normal((n, 3))


def _on(xp, spec=None, min_sites=0):
    ew = Ewald(8.0, spec, xp=xp)
    ew.device_min_sites = min_sites
    return ew


def _same(host, device, rel=1e-12):
    for name in ENERGIES:
        assert getattr(device, name) == pytest.approx(getattr(host, name), rel=rel, abs=1e-12), name
    for name in ARRAYS:
        a, b = getattr(host, name), getattr(device, name)
        assert (a is None) == (b is None), name
        if a is not None:
            assert isinstance(b, np.ndarray), f"{name} came back as {type(b)}"
            np.testing.assert_allclose(b, a, rtol=0, atol=rel * max(1.0, float(np.abs(a).max())), err_msg=name)


# ------------------------------------------------------------------ the double has teeth
def test_the_double_refuses_what_a_device_array_refuses():
    a = dd.asarray(np.arange(6.0).reshape(2, 3))
    with pytest.raises(TypeError, match="Implicit conversion"):
        np.asarray(a)
    for mix in (lambda: a + np.ones(3), lambda: np.ones(3) + a, lambda: a @ np.ones(3), lambda: np.ones(2) @ a,
                lambda: dd.concatenate([a, np.ones((2, 3))]), lambda: dd.sqrt(np.ones(3)), lambda: a[np.array([0, 1])],
                lambda: dd.bincount(np.array([0, 1])), lambda: dd.linalg.solve(np.eye(2), dd.asarray(np.ones(2)))):
        with pytest.raises(TypeError, match="Unsupported type"):
            mix()
    assert isinstance(a.sum(), dd.DeviceArray) and a.sum().ndim == 0  # a reduction is not a scalar
    assert float(a.sum()) == 15.0 and float(2.0 * a.sum()) == 30.0
    assert isinstance(np.float64(2.0) * a, dd.DeviceArray) and isinstance(np.exp(a), dd.DeviceArray)
    assert not hasattr(a, "setflags")
    view = a[:, 1]
    view += 10.0  # a slice is a view, as it is on a device
    assert a.get()[:, 1].tolist() == [11.0, 14.0]


# ------------------------------------------------------------------ the sum
def test_the_host_is_the_default_and_a_device_is_asked_for(monkeypatch):
    monkeypatch.delenv("POLYFIND_EWALD_DEVICE", raising=False)
    assert bk.ewald_backend() is np and Ewald(8.0).xp is np
    monkeypatch.setenv("POLYFIND_EWALD_DEVICE", "cpu")
    assert bk.ewald_backend() is np
    monkeypatch.setenv("POLYFIND_EWALD_DEVICE", "gpu")
    with pytest.raises(ValueError, match="POLYFIND_EWALD_DEVICE"):
        bk.ewald_backend()
    # asked for, it is CuPy or an error -- never a silent fall back to the host
    monkeypatch.setenv("POLYFIND_EWALD_DEVICE", "cuda")
    try:
        import cupy
    except Exception:
        with pytest.raises(ImportError):
            bk.ewald_backend()
    else:
        assert bk.ewald_backend() is cupy


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
@pytest.mark.parametrize("damping", [None, 4.0])
def test_every_term_on_a_device_is_the_hosts(boundary, damping):
    X, q, h, p, scale = _cell()
    spec = EwaldSpec(boundary=boundary, charge_damping=damping)
    for th in (Thole(), Thole(a=0.39, form="cubic"), None):
        kw = {} if th is None else dict(thole=th, thole_scale=scale)
        for flags in (dict(dipoles=p, grad=True, charge_grad=True, dipole_grad=True, field=True, tensor=True),
                      dict(dipoles=p, grad=True), dict(field=True, tensor=True), dict(dipole_grad=True)):
            host = _on(np, spec).terms(X, q, h, **flags, **kw)
            _same(host, _on(dd, spec).terms(X, q, h, **flags, **kw))
    for flags in (dict(), dict(grad=True), dict(grad=True, charge_grad=True)):  # the charge-only path
        _same(_on(np, spec).terms(X, q, h, **flags), _on(dd, spec).terms(X, q, h, **flags))
    # a charged cell: the neutralising background and its derivatives
    _same(_on(np, spec).terms(X, q + 0.01, h, grad=True, charge_grad=True),
          _on(dd, spec).terms(X, q + 0.01, h, grad=True, charge_grad=True))


def test_only_the_tensor_can_stay_on_the_device_and_only_when_asked():
    X, q, h, p, scale = _cell()
    kw = dict(thole=Thole(), thole_scale=scale, field=True, tensor=True)
    back = _on(dd).terms(X, q, h, **kw)
    kept = _on(dd).terms(X, q, h, device_tensor=True, **kw)
    assert isinstance(back.tensor, np.ndarray) and isinstance(kept.tensor, dd.DeviceArray)
    assert isinstance(kept.field, np.ndarray)
    np.testing.assert_array_equal(kept.tensor.get(), back.tensor)
    # on the host the flag changes nothing
    assert isinstance(_on(np).terms(X, q, h, device_tensor=True, **kw).tensor, np.ndarray)


def test_nothing_of_the_cells_size_squared_crosses_in_a_call():
    """What goes to the device is the geometry and the k vectors; what comes back is per site.

    The (N, N) table of Thole lengths is the one large thing a caller owns: handed over as a
    device array it is not sent again."""
    n = 120
    X, q, h, p, scale = _cell(n)
    ew = _on(dd)
    scale_dev = dd.asarray(scale)
    kw = dict(thole=Thole(), thole_scale=scale_dev)
    dd.reset_traffic()
    t = ew.terms(X, q, h, field=True, tensor=True, device_tensor=True, **kw)
    n_k = t.n_kvec
    per_site, vectors, squared = 8 * 3 * n, 8 * 3 * n_k, 8 * n * n
    assert max(per_site, vectors) < squared  # so the two bounds below say something
    assert dd.traffic["largest_to_host"] <= per_site  # a field, never the tensor or a list of pairs
    assert dd.traffic["largest_to_device"] <= max(per_site, vectors)  # the coordinates, the k vectors
    dd.reset_traffic()
    ew.terms(X, q, h, dipoles=p, grad=True, charge_grad=True, **kw)
    assert dd.traffic["largest_to_host"] <= per_site
    # in all: X, q, p, the lattice and its inverse, the image integers, and k, k^2, A_k
    assert dd.traffic["to_device_bytes"] <= 8 * (7 * n + 18 + 3 * t.n_image + 5 * n_k)


def test_a_small_cell_stays_on_the_host_whatever_the_module():
    X, q, h, p, scale = _cell(24)
    ew = Ewald(8.0, xp=dd)
    assert ew.device_min_sites == 256 and ew.array_module(24) is np and ew.array_module(256) is dd
    dd.reset_traffic()
    t = ew.terms(X, q, h, dipoles=p, thole=Thole(), thole_scale=scale, grad=True, tensor=True, device_tensor=True)
    assert dd.traffic["to_device"] == 0 and isinstance(t.tensor, np.ndarray)
    ref = Ewald(8.0).terms(X, q, h, dipoles=p, thole=Thole(), thole_scale=scale, grad=True, tensor=True)
    assert t.total == ref.total and np.array_equal(t.grad_coords, ref.grad_coords) and np.array_equal(t.tensor, ref.tensor)


def test_the_phases_are_kept_per_module_and_dropped_on_request():
    X, q, h, p = _wide_cell()
    ew = _on(dd)
    first = ew.terms(X, q, h)
    assert first.n_kvec * len(q) >= 1_000_000 and ew._phase_cache is not None
    kept = ew._phase_cache[2]
    assert isinstance(kept, dd.DeviceArray) and ew._phase_cache[4] is dd
    ew.terms(X, q, h, dipoles=p, grad=True)
    assert ew._phase_cache[2] is kept  # the second call of an evaluation reuses them
    ew.xp = np  # the same geometry on the host: the device's arrays are not handed to NumPy
    host = ew.terms(X, q, h)
    assert isinstance(ew._phase_cache[2], np.ndarray)
    assert host.total == pytest.approx(first.total, rel=1e-12)
    ew.release()
    assert ew._phase_cache is None


def test_an_ewald_object_pickles_without_its_module_or_its_phases():
    X, q, h, p = _wide_cell()
    ew = Ewald(8.0)
    e0 = ew.terms(X, q, h).total
    assert ew._phase_cache is not None
    twin = pickle.loads(pickle.dumps(ew))
    assert twin.xp is np and twin._phase_cache is None and twin.device_min_sites == ew.device_min_sites
    assert twin.terms(X, q, h).total == e0
    assert ew._phase_cache is not None  # pickling took nothing from the original


# ------------------------------------------------------------------ the solve and the whole evaluation
def _model(nx, ny):
    polymer = pvdf_belt(0.7, "pbe")
    seq = [T, GP, T, GM]
    tors = np.array([THREE_STATE.angles[s] for s in seq], dtype=float)
    ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float)
    name = "".join(THREE_STATE.names[s] for s in seq)
    ch = repeat_chains(polymer, name, tors[None], ang0)[0]
    nc = 2 * nx * ny
    pk = pack_mod.CrystalPacker(ch, n_chains=nc, **EWALD)
    st, flips = MC.supercell_state([5.0, 9.5, 91.0, 30.0, 215.0, 2.2, 1.0], ch.c, tors, ang0, nx, ny)
    return MC.MixedModel(polymer, name, pk, flips, [ch] * nc, angle_stiffness=105.0), st.vector()


def test_a_mixed_cell_refuses_a_packer_of_the_batched_device_backend():
    """With a GPU in the machine the packer's batched kernels go to it in float32 unasked; a mixed
    cell built on such a packer would run its float64 host loops through them.  It says so."""
    with PRESET_PARAMETERS[BELT4].applied():
        m, x = _model(1, 1)
        for attr, value in (("_dt", np.float32), ("xp", dd)):
            old = getattr(m.pk, attr)
            setattr(m.pk, attr, value)
            try:
                with pytest.raises(ValueError, match="POLYFIND_DEVICE=cpu"):
                    MC.MixedModel(m.polymer, m.name, m.pk, m.flips, [m.refs[0]] * m.nc, angle_stiffness=105.0)
            finally:
                setattr(m.pk, attr, old)


def test_the_induced_dipole_evaluation_on_a_device_is_the_hosts():
    """The whole energy and gradient of a cell of four chains: sums, solve and gradient on the
    device, and neither the tensor nor the matrix of the solve brought back."""
    with PRESET_PARAMETERS[BELT4].applied():
        m, x = _model(2, 1)
        rng = np.random.default_rng(0)
        x = x + 0.02 * rng.standard_normal(x.size) * m.scale(x)
        e_host, g_host = m.value_and_grad(x)
        ew = m.pk._ewald
        na = m.pk.N if m.pk._es_idx is None else len(m.pk._es_idx)
        assert ew.xp is np and na < ew.device_min_sites
        ew.xp, ew.device_min_sites = dd, 0
        try:
            dd.reset_traffic()
            e_dev, g_dev = m.value_and_grad(x)
            first = dict(dd.traffic)
            assert isinstance(m.pk._thole_scale_dev[2], dd.DeviceArray)
            dd.reset_traffic()
            m.value_and_grad(x)
            second = dict(dd.traffic)
        finally:
            ew.xp, ew.device_min_sites = np, 256
            ew.release()
    assert e_dev == pytest.approx(e_host, rel=1e-12)
    np.testing.assert_allclose(g_dev, g_host, rtol=0, atol=1e-10 * max(1.0, float(np.abs(g_host).max())))
    assert isinstance(g_dev, np.ndarray)
    # the largest thing brought back is a per-site array, not the (3 Na)^2 matrix or the tensor
    assert first["largest_to_host"] <= 8 * 3 * na < 8 * 9 * na * na
    # the Thole table went over once, with the first evaluation, and not again
    assert first["largest_to_device"] == 8 * na * na and second["largest_to_device"] < 8 * na * na
