"""Enumerate, deduplicate and rank periodic torsion sequences (crystalline chains).

Uses the cyclic k-best dynamic program of :class:`polyfind.ris.RISModel` per
period, then collapses symmetry-equivalent sequences (shift by a repeat unit,
and -- for an achiral repeat -- reversal and mirror; for a chiral one only their
composition, see :func:`polyfind.helix.sequence_images`), discards non-primitive
repetitions of shorter periods,
computes helix descriptors, and filters out chains that do not propagate
(rings, rise per bond below ``min_rise_per_bond``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .helix import HelixParams, canonical_sequence, helix_parameters, is_primitive
from .polymers import Polymer
from .ris import RISModel

# Known polymorph chain conformations (canonical forms are compared, so any
# equivalent writing of the sequence matches).
KNOWN_CHAINS: dict[str, dict[str, str]] = {
    "pvdf": {"beta (TTTT)": "TT", "alpha/delta (TGTG')": "TG+TG-", "gamma/epsilon (T3GT3G')": "TTTG+TTTG-"},
    "pe": {"orthorhombic (all-trans)": "T"},
}


@dataclass
class Candidate:
    seq: tuple[int, ...]
    name: str
    period: int
    energy_per_period: float
    energy_per_monomer: float
    helix: HelixParams
    known_as: str | None = None
    rank: int = -1

    def row(self) -> str:
        known = f"  <- {self.known_as}" if self.known_as else ""
        c = f"{self.helix.c:6.2f}" if self.helix.c is not None else "  n/a "
        return (
            f"{self.rank:3d}  {self.name:<18s} P={self.period:2d}  E/mon={self.energy_per_monomer:7.3f}  "
            f"c={c} A  theta={self.helix.rotation_per_period:7.1f}  r={self.helix.radius_all:4.2f}  {self.helix.label}{known}"
        )


def enumerate_periodic(
    polymer: Polymer,
    model: RISModel,
    max_period: int = 8,
    k_per_period: int = 60,
    min_rise_per_bond: float = 0.3,
    require_commensurate: bool = False,
    max_m: int = 8,
) -> list[Candidate]:
    """Ranked list of distinct periodic chain conformations up to ``max_period`` bonds.

    Helices needing more than ``max_m`` sequence periods per crystallographic repeat are
    kept (with ``helix.c = None``) unless ``require_commensurate``; they are real chain
    conformations but too large to pack cheaply.

    Deduplication uses the symmetry group of ``polymer``'s repeat: for a chiral one
    (:attr:`polyfind.polymers.Polymer.is_chiral`) a sequence and its G+/G- mirror are
    *distinct* conformations -- a right- and a left-handed helix with different
    energies -- so both are kept.  See :func:`polyfind.helix.sequence_images`.
    """
    B = polymer.bonds_per_repeat
    chiral = bool(polymer.is_chiral)
    states = model.states
    seen: dict[tuple, Candidate] = {}
    for P in range(B, max_period + 1, B):
        for e, seq in model.cyclic_k_best(P, k_per_period):
            seq = tuple(int(s) for s in seq)
            if not is_primitive(seq, B):
                continue
            key = canonical_sequence(seq, states, B, chiral=chiral)
            if key in seen:
                continue
            h = helix_parameters(polymer, np.array(key), states, max_m=max_m)
            if h.rise_per_bond < min_rise_per_bond:
                continue
            if require_commensurate and h.c is None:
                continue
            cand = Candidate(
                seq=key,
                name="".join(states.names[s] for s in key),
                period=P,
                energy_per_period=float(e),
                energy_per_monomer=float(e) / (P / B),
                helix=h,
            )
            seen[key] = cand
    known = {
        canonical_sequence(model.parse(s), states, B, chiral=chiral): lbl
        for lbl, s in KNOWN_CHAINS.get(polymer.name, {}).items()
    }
    out = sorted(seen.values(), key=lambda c: (c.energy_per_monomer, c.period))
    for i, c in enumerate(out):
        c.rank = i + 1
        c.known_as = known.get(c.seq)
    return out


def table(cands: list[Candidate], top: int | None = None) -> str:
    rows = [c.row() for c in (cands[:top] if top else cands)]
    hdr = "rank  sequence           period  E/monomer(kcal/mol)  c  rotation/period  radius  description"
    return hdr + "\n" + "\n".join(rows)
