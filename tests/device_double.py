"""A stand-in for a device array module, for testing the device path on a machine without one.

``DeviceArray`` is a NumPy array underneath, wrapped so that it behaves the way CuPy's arrays
do in the places where a port from NumPy goes wrong (each rule below was read off CuPy 14.2's
own source, ``cupy/_core/core.pyx``, and its guide "Differences between CuPy and NumPy"):

* it does not turn into a NumPy array by itself: ``np.asarray(a)`` raises, as CuPy's
  ``__array__`` does ("Implicit conversion to a NumPy array is not allowed"); the only way to
  the host is ``a.get()``, or ``float(a)`` / ``int(a)`` / ``bool(a)`` for a single number;
* arithmetic with a NumPy array on either side raises ("Unsupported type <class
  'numpy.ndarray'>"); Python and NumPy *scalars* are fine;
* functions of the module take device arrays only -- no NumPy arrays, no lists;
* a reduction returns a zero-dimensional device array, not a scalar;
* there is no ``setflags``;
* a NumPy ufunc or function applied to one is handed to the module's function of that name.

It is stricter than CuPy in one place: a NumPy array is refused as an *index* as well.  Code
that passes here mixes host and device nowhere, which is the point.

Because the arithmetic underneath is NumPy's, a routine written against an array module gives
the same numbers here as on the host wherever it does the same operations, and the module
counts what crosses between "host" and "device" (:data:`traffic`), so a test can assert that
nothing of the size of the cell's large arrays is carried across.

Use the module itself as the array module: ``import device_double as xp``.
"""
from __future__ import annotations

import numpy as np
from scipy import special as _special

__all__ = ["DeviceArray", "traffic", "reset_traffic"]

traffic = {"to_device": 0, "to_device_bytes": 0, "to_host": 0, "to_host_bytes": 0, "largest_to_host": 0,
           "largest_to_device": 0}


def reset_traffic() -> None:
    for key in traffic:
        traffic[key] = 0


_SCALARS = (bool, int, float, np.generic)


def _refuse(x):
    raise TypeError(f"Unsupported type {type(x)}")


def _raw(x):
    """What an operand of a device operation is underneath; a host array or a list is refused."""
    if isinstance(x, DeviceArray):
        return x._a
    if isinstance(x, _SCALARS):
        return x
    _refuse(x)


def _wrap(r):
    if isinstance(r, tuple):
        return tuple(_wrap(v) for v in r)
    return DeviceArray(np.asarray(r))


def _index(idx):
    def one(i):
        if isinstance(i, DeviceArray):
            return i._a
        if i is None or i is Ellipsis or isinstance(i, (slice, int, np.integer)):
            return i
        _refuse(i)

    return tuple(one(i) for i in idx) if isinstance(idx, tuple) else one(idx)


class DeviceArray:
    __array_priority__ = 100
    __slots__ = ("_a",)

    def __init__(self, a):
        self._a = a

    # --- the ways to the host -----------------------------------------------------------
    def get(self):
        traffic["to_host"] += 1
        traffic["to_host_bytes"] += self._a.nbytes
        traffic["largest_to_host"] = max(traffic["largest_to_host"], self._a.nbytes)
        return np.array(self._a)

    def item(self):
        return self.get().item()

    def __array__(self, dtype=None, copy=None):
        raise TypeError("Implicit conversion to a NumPy array is not allowed. "
                        "Please use `.get()` to construct a NumPy array explicitly.")

    def __float__(self):
        return float(self.get())

    def __int__(self):
        return int(self.get())

    def __bool__(self):
        return bool(self.get())

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented
        fn = globals().get(ufunc.__name__)
        if fn is None:
            raise TypeError(f"the device double has no {ufunc.__name__}")
        if "out" in kwargs:
            (kwargs["out"],) = kwargs["out"]
        return fn(*inputs, **kwargs)

    def __array_function__(self, func, types, args, kwargs):
        fn = globals().get(func.__name__)
        if fn is None:
            raise TypeError(f"the device double has no {func.__name__}")
        return fn(*args, **kwargs)

    # --- shape ---------------------------------------------------------------------------
    shape = property(lambda self: self._a.shape)
    ndim = property(lambda self: self._a.ndim)
    size = property(lambda self: self._a.size)
    dtype = property(lambda self: self._a.dtype)
    nbytes = property(lambda self: self._a.nbytes)
    T = property(lambda self: DeviceArray(self._a.T))

    def __len__(self):
        return len(self._a)

    def __repr__(self):
        return f"DeviceArray(shape={self._a.shape}, dtype={self._a.dtype})"

    def reshape(self, *shape):
        return DeviceArray(self._a.reshape(*shape))

    def transpose(self, *axes):
        return DeviceArray(self._a.transpose(*axes))

    def ravel(self):
        return DeviceArray(self._a.ravel())

    def copy(self):
        return DeviceArray(self._a.copy())

    def astype(self, dtype, copy=True):
        return DeviceArray(self._a.astype(dtype, copy=copy))

    def fill(self, value):
        self._a.fill(_raw(value))

    # --- reductions: zero-dimensional device arrays, never scalars ------------------------
    def sum(self, axis=None):
        return _wrap(self._a.sum(axis=axis))

    def max(self, axis=None):
        return _wrap(self._a.max(axis=axis))

    def min(self, axis=None):
        return _wrap(self._a.min(axis=axis))

    def any(self, axis=None):
        return _wrap(self._a.any(axis=axis))

    def all(self, axis=None):
        return _wrap(self._a.all(axis=axis))

    # --- indexing -------------------------------------------------------------------------
    def __getitem__(self, idx):
        return DeviceArray(np.asarray(self._a[_index(idx)]))

    def __setitem__(self, idx, value):
        self._a[_index(idx)] = _raw(value)

    # --- arithmetic -------------------------------------------------------------------------
    def __neg__(self):
        return DeviceArray(-self._a)

    def __pos__(self):
        return self

    def __abs__(self):
        return DeviceArray(np.abs(self._a))

    def __invert__(self):
        return DeviceArray(~self._a)


def _binary(name, fn):
    def op(self, other):
        return _wrap(fn(self._a, _raw(other)))

    def rop(self, other):
        return _wrap(fn(_raw(other), self._a))

    def iop(self, other):
        self._a = fn(self._a, _raw(other), out=self._a) if isinstance(fn, np.ufunc) else fn(self._a, _raw(other))
        return self

    setattr(DeviceArray, f"__{name}__", op)
    setattr(DeviceArray, f"__r{name}__", rop)
    setattr(DeviceArray, f"__i{name}__", iop)


for _name, _fn in (("add", np.add), ("sub", np.subtract), ("mul", np.multiply), ("truediv", np.true_divide),
                   ("pow", np.power), ("and", np.bitwise_and), ("or", np.bitwise_or), ("matmul", np.matmul)):
    _binary(_name, _fn)
# a product of matrices has no in-place form worth the name
DeviceArray.__imatmul__ = lambda self, other: _wrap(np.matmul(self._a, _raw(other)))
for _name, _fn in (("lt", np.less), ("le", np.less_equal), ("gt", np.greater), ("ge", np.greater_equal),
                   ("eq", np.equal), ("ne", np.not_equal)):
    setattr(DeviceArray, f"__{_name}__", (lambda fn: lambda self, other: _wrap(fn(self._a, _raw(other))))(_fn))
DeviceArray.__hash__ = None


# ------------------------------------------------------------------ the module's functions
def asarray(a, dtype=None):
    """Host to device (an array already there is returned as it is, as CuPy's does)."""
    if isinstance(a, DeviceArray):
        return a if dtype is None or a._a.dtype == np.dtype(dtype) else DeviceArray(a._a.astype(dtype))
    host = np.array(a, dtype=dtype)
    traffic["to_device"] += 1
    traffic["to_device_bytes"] += host.nbytes
    traffic["largest_to_device"] = max(traffic["largest_to_device"], host.nbytes)
    return DeviceArray(host)


def asnumpy(a):
    return a.get() if isinstance(a, DeviceArray) else np.asarray(a)


def zeros(shape, dtype=float):
    return DeviceArray(np.zeros(shape, dtype=dtype))


def eye(n, dtype=float):
    return DeviceArray(np.eye(n, dtype=dtype))


def arange(*args, dtype=None):
    return DeviceArray(np.arange(*args, dtype=dtype))


def empty_like(a):
    return DeviceArray(np.empty_like(_raw(a)))


def _unary(fn):
    return lambda a: _wrap(fn(_raw(a)))


sqrt, exp, cos, sin, absolute, negative = (_unary(f) for f in (np.sqrt, np.exp, np.cos, np.sin, np.abs, np.negative))
round = _unary(np.round)  # noqa: A001  (the name the array modules use)
erfc = _unary(_special.erfc)


def _binary_function(fn):
    def call(a, b, out=None):
        if out is None:
            return _wrap(fn(_raw(a), _raw(b)))
        fn(_raw(a), _raw(b), out=out._a)
        return out

    return call


add, subtract, multiply, true_divide, power, matmul = (
    _binary_function(f) for f in (np.add, np.subtract, np.multiply, np.true_divide, np.power, np.matmul))
divide = true_divide
less, less_equal, greater, greater_equal, equal, not_equal, bitwise_and, bitwise_or = (
    _binary_function(f) for f in (np.less, np.less_equal, np.greater, np.greater_equal, np.equal, np.not_equal,
                                  np.bitwise_and, np.bitwise_or))


def einsum(subscripts, *operands):
    return _wrap(np.einsum(subscripts, *(_raw(o) for o in operands)))


def nonzero(a):
    return _wrap(np.nonzero(_raw(a)))


def bincount(x, weights=None, minlength=None):
    if not isinstance(x, DeviceArray) or (weights is not None and not isinstance(weights, DeviceArray)):
        _refuse(x if not isinstance(x, DeviceArray) else weights)
    if x._a.dtype.kind == "f":
        raise TypeError("x must be int array")
    return DeviceArray(np.bincount(x._a, weights=None if weights is None else weights._a,
                                   minlength=0 if minlength is None else int(minlength)))


def stack(arrays, axis=0):
    return DeviceArray(np.stack([_raw(a) for a in arrays], axis=axis))


def concatenate(arrays, axis=0):
    for a in arrays:
        if not isinstance(a, DeviceArray):
            _refuse(a)
    return DeviceArray(np.concatenate([a._a for a in arrays], axis=axis))


class linalg:  # noqa: N801  (a namespace, as numpy.linalg is)
    @staticmethod
    def solve(a, b):
        if not isinstance(a, DeviceArray) or not isinstance(b, DeviceArray):
            _refuse(a if not isinstance(a, DeviceArray) else b)
        return DeviceArray(np.linalg.solve(a._a, b._a))
