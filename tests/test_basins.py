"""examples/copolymer_basins.py: the packing screen's minima refined under the full potential (docs/CRYSTAL_STAGE.md 5.17)."""
import importlib.util
from pathlib import Path

import numpy as np

from polyfind import fitting as F
from polyfind import mechanics as M
from polyfind.pack import pack, periodic_chain
from polyfind.polymers import PVDF, THREE_STATE

SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "copolymer_basins.py"
SPEC = importlib.util.spec_from_file_location("copolymer_basins", SCRIPT)
CB = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CB)


def test_beta_has_two_basins_and_the_first_is_the_e1_reference():
    """Under the pre-registered potential the screen finds two distinct minima for beta: the DFT
    setting (dipoles along a, the cell every set stands on) and the chains turned by about 60 deg;
    refined under the full potential they are within 0.01 kcal/mol per monomer, and refining the
    first through the script's tail reproduces mechanics.refined_reference's cell."""
    with F.FITTED_VALENCE.applied():
        chain = periodic_chain(PVDF, [0, 0], THREE_STATE)
        cands = pack(chain, n_chains=2, n_refine=5, table_cache_dir=None)
        assert len(cands) >= 2
        ref0, rr0 = CB.refine_candidate(PVDF, cands[0], "beta basin 0")
        ref1, rr1 = CB.refine_candidate(PVDF, cands[1], "beta basin 1")
        d0, d1 = CB.describe(ref0, rr0, cands[0]), CB.describe(ref1, rr1, cands[1])
        # the E1 path, same potential: the same cell to 1e-3 A
        e1, _ = M.refined_reference(PVDF, [0, 0], label="beta", refine_kw={k: CB.KW[k] for k in CB.KW},
                                    pack_kw={"table_cache_dir": None}, **CB.KW)
        e1 = M.relax_reference(e1)
    f0, f1 = d0["refined"], d1["refined"]
    assert abs(f0["a"] - e1.params[0]) < 1e-3 and abs(f0["b"] - e1.params[1]) < 1e-3
    assert abs(abs(f0["tilt_from_a_deg"]) - 180.0) < 2.0  # dipoles along a (the sign is the cell's orientation)
    assert 100.0 < abs(f1["tilt_from_a_deg"]) < 140.0  # the chains turned by about 60 deg
    assert abs(f0["total_per_monomer"] - f1["total_per_monomer"]) < 0.01
    assert abs(f0["P_abs"] - f1["P_abs"]) < 1e-3  # the same chains, the same dipole per cell
