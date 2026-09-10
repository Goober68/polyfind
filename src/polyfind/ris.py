"""Rotational isomeric state (RIS) model with exact dynamic-programming solvers.

A chain of ``N`` backbone bonds carries one discrete state per bond
(e.g. T, G+, G-).  The energy is a first/second(/third)-neighbour Markov field::

    E(s) = sum_i e1[t(i), s_i] + sum_i e2[t(i), s_i, s_{i+1}] (+ sum_i e3[t(i), s_i, s_{i+1}, s_{i+2}])

where ``t(i) = i mod B`` is the bond type inside the chemical repeat unit
(``B`` = bonds per repeat; PVDF has B = 2 so that CH2-centred and CF2-centred
dihedral pairs get different pair matrices).

Because the energy is a chain-structured Markov field, every quantity below is
exact and costs O(N * S^2) (O(N * S^4) with third-order terms, which are handled
by *state augmentation*: pairs of consecutive states become the states of an
equivalent second-order model):

* :meth:`RISModel.minimum`        - global minimum-energy sequence (Viterbi)
* :meth:`RISModel.k_best`         - the k lowest-energy sequences
* :meth:`RISModel.cyclic_k_best`  - the k lowest-energy *periodic* sequences (crystalline chains)
* :meth:`RISModel.log_partition`  - partition function / conformational free energy
* :meth:`RISModel.marginals`      - per-bond state probabilities (forward-backward)
* :meth:`RISModel.sample`         - exact Boltzmann samples, batched, optionally clamped

Energies are in kcal/mol, temperatures in K.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.special import logsumexp

from . import backend as bk
from .polymers import RISStates, THREE_STATE

KB = 0.0019872041  # kcal / (mol K)
INF = float("inf")


def _mask(n_bonds: int, S: int, clamp) -> np.ndarray:
    """(n, S) additive log-mask: 0 allowed, -inf forbidden.  clamp: {pos: int | iterable}."""
    mask = np.zeros((n_bonds, S))
    if clamp:
        for i, allowed in clamp.items():
            mask[i] = -INF
            if np.isscalar(allowed):
                mask[i, int(allowed)] = 0.0
            else:
                for s in allowed:
                    mask[i, int(s)] = 0.0
    return mask


@dataclass
class RISModel:
    states: RISStates
    bonds_per_repeat: int
    first_order: np.ndarray  # (B, S)
    second_order: np.ndarray  # (B, S, S): pair (bond type b, state s) -> (bond type b+1, state s')
    third_order: np.ndarray | None = None  # (B, S, S, S): triple starting at bond type b
    name: str = ""

    def __post_init__(self):
        self.first_order = np.asarray(self.first_order, dtype=float)
        self.second_order = np.asarray(self.second_order, dtype=float)
        B, S = self.bonds_per_repeat, self.states.n
        assert self.first_order.shape == (B, S), self.first_order.shape
        assert self.second_order.shape == (B, S, S), self.second_order.shape
        if self.third_order is not None:
            self.third_order = np.asarray(self.third_order, dtype=float)
            assert self.third_order.shape == (B, S, S, S), self.third_order.shape
            if not np.any(self.third_order):
                self.third_order = None
        self._aug_cache = None

    # ------------------------------------------------------------------ basics
    @property
    def S(self) -> int:
        return self.states.n

    @property
    def B(self) -> int:
        return self.bonds_per_repeat

    @property
    def order(self) -> int:
        return 3 if self.third_order is not None else 2

    def bond_type(self, i: int) -> int:
        return i % self.B

    def energy(self, seq: Sequence[int], periodic: bool = False) -> float:
        """Energy of a state sequence (indices). ``periodic`` closes the ring."""
        s = np.asarray(seq, dtype=int)
        n = len(s)
        if periodic and n % self.B != 0:
            raise ValueError("periodic sequence length must be a multiple of bonds_per_repeat")
        t = np.arange(n) % self.B
        e = float(self.first_order[t, s].sum())
        if periodic:
            s2 = np.concatenate([s, s[:2]])
            e += float(self.second_order[t, s2[:n], s2[1 : n + 1]].sum())
            if self.third_order is not None:
                e += float(self.third_order[t, s2[:n], s2[1 : n + 1], s2[2 : n + 2]].sum())
        else:
            if n > 1:
                e += float(self.second_order[t[:-1], s[:-1], s[1:]].sum())
            if n > 2 and self.third_order is not None:
                e += float(self.third_order[t[:-2], s[:-2], s[1:-1], s[2:]].sum())
        return e

    def names(self, seq: Sequence[int]) -> str:
        return "".join(self.states.names[i] for i in seq)

    def parse(self, text: str) -> list[int]:
        """Parse e.g. 'TG+TG-' into state indices (longest-match tokenizer)."""
        out, i = [], 0
        names = sorted(self.states.names, key=len, reverse=True)
        while i < len(text):
            for nm in names:
                if text.startswith(nm, i):
                    out.append(self.states.index(nm))
                    i += len(nm)
                    break
            else:
                raise ValueError(f"cannot parse state at {text[i:]!r}")
        return out

    # ------------------------------------------------------ third order support
    def _aug(self) -> "RISModel":
        """Equivalent second-order model over pair-states q = (s_j, s_{j+1})."""
        if self._aug_cache is not None:
            return self._aug_cache
        S, B = self.S, self.B
        Q = S * S
        e1, e2, e3 = self.first_order, self.second_order, self.third_order
        f1 = np.zeros((B, Q))
        f2 = np.full((B, Q, Q), INF)
        for t in range(B):
            for a in range(S):
                for b in range(S):
                    q = a * S + b
                    f1[t, q] = e1[t, a] + e2[t, a, b]
                    for c in range(S):
                        f2[t, q, b * S + c] = e3[t, a, b, c]
        names = tuple(f"{self.states.names[a]}|{self.states.names[b]}" for a in range(S) for b in range(S))
        mirror = tuple(self.states.mirror[a] * S + self.states.mirror[b] for a in range(S) for b in range(S))
        aug = RISModel(RISStates(names, (0.0,) * Q, mirror), B, f1, f2, name=self.name + "-aug")
        self._aug_cache = aug
        return aug

    def _aug_end_term(self, n_bonds: int) -> np.ndarray:
        """First-order energy of the last bond, indexed by pair-state (a, b) -> e1[t_{N-1}, b]."""
        S = self.S
        return np.tile(self.first_order[(n_bonds - 1) % self.B], S)

    def _aug_clamp(self, n_bonds: int, clamp):
        if not clamp:
            return None
        S = self.S
        out = {}
        for j in range(n_bonds - 1):
            if j in clamp or (j + 1) in clamp:
                allowed = [a * S + b for a in range(S) for b in range(S)
                           if (j not in clamp or a == clamp[j]) and ((j + 1) not in clamp or b == clamp[j + 1])]
                out[j] = allowed
        return out

    def _decode(self, qseq: np.ndarray) -> np.ndarray:
        """Pair-state sequence (n-1,) or (M, n-1) -> state sequence (n,) or (M, n)."""
        S = self.S
        qseq = np.asarray(qseq)
        if qseq.ndim == 1:
            return np.concatenate([qseq // S, [qseq[-1] % S]])
        return np.concatenate([qseq // S, (qseq[:, -1] % S)[:, None]], axis=1)

    # ---------------------------------------------------------------- viterbi
    def minimum(self, n_bonds: int, end_term: np.ndarray | None = None) -> tuple[float, np.ndarray]:
        """Global minimum-energy sequence of ``n_bonds`` bonds. O(N S^2)."""
        if self.third_order is not None and n_bonds >= 2:
            e, q = self._aug().minimum(n_bonds - 1, end_term=self._aug_end_term(n_bonds))
            return e, self._decode(q)
        e1, e2, S = self.first_order, self.second_order, self.S
        V = e1[0].copy()
        back = np.zeros((n_bonds, S), dtype=int)
        for i in range(1, n_bonds):
            cand = V[:, None] + e2[(i - 1) % self.B]  # (p, s)
            back[i] = cand.argmin(axis=0)
            V = cand.min(axis=0) + e1[i % self.B]
        if end_term is not None:
            V = V + end_term
        seq = np.empty(n_bonds, dtype=int)
        seq[-1] = int(V.argmin())
        for i in range(n_bonds - 1, 0, -1):
            seq[i - 1] = back[i, seq[i]]
        return float(V.min()), seq

    def k_best(self, n_bonds: int, k: int, start: int | None = None, end_pair_to: int | None = None, end_term: np.ndarray | None = None):
        """The ``k`` lowest-energy sequences (k-best Viterbi), ascending (energy, seq).

        ``start`` forces s_0; ``end_pair_to`` adds the closing pair term (cyclic use);
        ``end_term`` adds a per-state energy at the last position.
        """
        if self.third_order is not None and n_bonds >= 2 and start is None and end_pair_to is None:
            out = self._aug().k_best(n_bonds - 1, k, end_term=self._aug_end_term(n_bonds))
            return [(e, self._decode(q)) for e, q in out]
        e1, e2, S, B = self.first_order, self.second_order, self.S, self.B
        first = []
        for s in range(S):
            if start is not None and s != start:
                first.append((np.empty(0), np.empty((0, 2), dtype=int)))
            else:
                first.append((np.array([e1[0, s]]), np.array([[-1, -1]])))
        lists = [first]
        for i in range(1, n_bonds):
            t = (i - 1) % B
            cur = []
            for s in range(S):
                es, bps = [], []
                for p in range(S):
                    pe, _ = lists[i - 1][p]
                    if pe.size == 0 or not np.isfinite(e2[t, p, s]):
                        continue
                    es.append(pe + e2[t, p, s] + e1[i % B, s])
                    bps.append(np.stack([np.full(pe.size, p), np.arange(pe.size)], axis=1))
                if not es:
                    cur.append((np.empty(0), np.empty((0, 2), dtype=int)))
                    continue
                es = np.concatenate(es)
                bps = np.concatenate(bps)
                order = np.argsort(es, kind="stable")[:k]
                cur.append((es[order], bps[order]))
            lists.append(cur)
        finals = []
        for s in range(S):
            pe, _ = lists[-1][s]
            for r in range(pe.size):
                e = pe[r]
                if end_pair_to is not None:
                    e = e + e2[(n_bonds - 1) % B, s, end_pair_to]
                if end_term is not None:
                    e = e + end_term[s]
                if np.isfinite(e):
                    finals.append((float(e), s, r))
        finals.sort(key=lambda x: x[0])
        out = []
        for e, s, r in finals[:k]:
            seq = np.empty(n_bonds, dtype=int)
            i, st, rk = n_bonds - 1, s, r
            while i >= 0:
                seq[i] = st
                _, bps = lists[i][st]
                st, rk = int(bps[rk, 0]), int(bps[rk, 1])
                i -= 1
            out.append((e, seq))
        return out

    def cyclic_k_best(self, period: int, k: int):
        """Lowest-energy periodic sequences of a given period (cycles in state space).

        Each cycle is reported once per distinct starting state it passes through;
        callers should canonicalise under shift/reversal/mirror (see :mod:`polyfind.helix`).
        Energies are per period (all first-order, pair and triple terms, ring closed).
        """
        if period % self.B != 0:
            raise ValueError("period must be a multiple of bonds_per_repeat")
        if self.third_order is not None:
            S = self.S
            out = self._aug().cyclic_k_best(period, k)
            seen = {}
            for e, q in out:
                seq = tuple(int(x) for x in (np.asarray(q) // S))
                seen[seq] = e
            return sorted(((e, np.array(s)) for s, e in seen.items()), key=lambda x: x[0])
        seen: dict[tuple, float] = {}
        for s0 in range(self.S):
            if period == 1:
                e = float(self.first_order[0, s0] + self.second_order[0, s0, s0])
                if np.isfinite(e):
                    seen[(s0,)] = e
                continue
            for e, seq in self.k_best(period, k, start=s0, end_pair_to=s0):
                seen[tuple(int(x) for x in seq)] = e
        return sorted(((e, np.array(s)) for s, e in seen.items()), key=lambda x: x[0])

    def enumerate_all(self, n_bonds: int, periodic: bool = False):
        """Brute-force enumeration (testing only; S^N sequences)."""
        import itertools

        out = []
        for seq in itertools.product(range(self.S), repeat=n_bonds):
            out.append((self.energy(seq, periodic=periodic), np.array(seq)))
        out.sort(key=lambda x: x[0])
        return out

    # ---------------------------------------------------------- partition fn
    def _log_weights(self, T: float):
        beta = 1.0 / (KB * T)
        return -beta * self.first_order, -beta * self.second_order

    def forward(self, n_bonds: int, T: float, clamp=None, end_term: np.ndarray | None = None) -> np.ndarray:
        """Forward log-messages log alpha[i, s] (unnormalised log-probabilities of prefixes)."""
        lw1, lw2 = self._log_weights(T)
        S, B = self.S, self.B
        mask = _mask(n_bonds, S, clamp)
        la = np.empty((n_bonds, S))
        la[0] = lw1[0] + mask[0]
        for i in range(1, n_bonds):
            la[i] = logsumexp(la[i - 1][:, None] + lw2[(i - 1) % B], axis=0) + lw1[i % B] + mask[i]
        if end_term is not None:
            la[-1] = la[-1] - end_term / (KB * T)
        return la

    def log_partition(self, n_bonds: int, T: float, clamp=None) -> float:
        if self.third_order is not None and n_bonds >= 2:
            return self._aug().log_partition_general(n_bonds - 1, T, self._aug_clamp(n_bonds, clamp), self._aug_end_term(n_bonds))
        return self.log_partition_general(n_bonds, T, clamp, None)

    def log_partition_general(self, n_bonds, T, clamp, end_term):
        return float(logsumexp(self.forward(n_bonds, T, clamp, end_term)[-1]))

    def free_energy(self, n_bonds: int, T: float) -> float:
        """Conformational free energy -kT ln Z (kcal/mol) of an ``n_bonds`` chain."""
        return -KB * T * self.log_partition(n_bonds, T)

    def infinite_chain_free_energy_per_repeat(self, T: float) -> float:
        """-kT ln(lambda_max) of the repeat-unit transfer matrix product (Flory's result).

        Conformational free energy per chemical repeat of the infinite chain.
        """
        if self.third_order is not None:
            return self._aug().infinite_chain_free_energy_per_repeat(T)
        lw1, lw2 = self._log_weights(T)
        U = np.eye(self.S)
        for b in range(self.B):
            Ub = np.exp(lw1[b][:, None] + lw2[b])
            U = U @ Ub
        lam = np.max(np.abs(np.linalg.eigvals(U)))
        return -KB * T * float(np.log(lam))

    def marginals(self, n_bonds: int, T: float, clamp=None) -> np.ndarray:
        """Per-bond state probabilities p[i, s] via forward-backward."""
        if self.third_order is not None and n_bonds >= 2:
            S = self.S
            pq = self._aug().marginals_general(n_bonds - 1, T, self._aug_clamp(n_bonds, clamp), self._aug_end_term(n_bonds))
            p = np.empty((n_bonds, S))
            p[:-1] = pq.reshape(n_bonds - 1, S, S).sum(axis=2)
            p[-1] = pq[-1].reshape(S, S).sum(axis=0)
            return p
        return self.marginals_general(n_bonds, T, clamp, None)

    def marginals_general(self, n_bonds, T, clamp, end_term):
        lw1, lw2 = self._log_weights(T)
        S, B = self.S, self.B
        la = self.forward(n_bonds, T, clamp, end_term)
        lb = np.zeros((n_bonds, S))
        mask = _mask(n_bonds, S, clamp)
        extra = np.zeros(S) if end_term is None else -end_term / (KB * T)
        for i in range(n_bonds - 2, -1, -1):
            nxt = lw1[(i + 1) % B] + mask[i + 1] + lb[i + 1]
            if i + 1 == n_bonds - 1:
                nxt = nxt + extra
            lb[i] = logsumexp(lw2[i % B] + nxt[None, :], axis=1)
        lp = la + lb
        lp -= logsumexp(lp, axis=1, keepdims=True)
        return np.exp(lp)

    # ---------------------------------------------------------------- sampling
    def sample(self, n_bonds: int, T: float, n_samples: int, rng=None, clamp=None):
        """Exact Boltzmann samples: (n_samples, n_bonds) int array of state indices.

        Forward pass is sequential (NumPy); the backward sampling step is batched over
        samples on the active backend (CuPy on GPU), so 10^5 chains cost about the same
        wall time as 10^2 on a GPU.
        ``clamp`` = {bond index: state index} fixes states (e.g. a crystalline stem).
        """
        if self.third_order is not None and n_bonds >= 2:
            q = self._aug().sample_general(n_bonds - 1, T, n_samples, rng, self._aug_clamp(n_bonds, clamp), self._aug_end_term(n_bonds))
            return self._decode(bk.to_numpy(q))
        return self.sample_general(n_bonds, T, n_samples, rng, clamp, None)

    def sample_general(self, n_bonds, T, n_samples, rng, clamp, end_term):
        xp = bk.get_backend()
        if rng is None:
            rng = bk.default_rng()
        la = xp.asarray(self.forward(n_bonds, T, clamp, end_term))
        lw2 = xp.asarray(self._log_weights(T)[1])
        S, B = self.S, self.B
        out = xp.empty((n_samples, n_bonds), dtype=xp.int64)
        u = rng.random((n_samples, n_bonds))

        def draw(logits, ui):
            m = logits.max(axis=1, keepdims=True)
            p = xp.exp(logits - m)
            c = xp.cumsum(p, axis=1)
            c = c / c[:, -1:]
            return xp.minimum((ui[:, None] > c).sum(axis=1), S - 1)

        out[:, -1] = draw(xp.broadcast_to(la[-1], (n_samples, S)), u[:, -1])
        for i in range(n_bonds - 2, -1, -1):
            nxt = out[:, i + 1]
            logits = la[i][None, :] + lw2[i % B][:, nxt].T  # (n, S)
            out[:, i] = draw(logits, u[:, i])
        return out

    # -------------------------------------------------------------------- io
    def to_dict(self) -> dict:
        d = {
            "name": self.name,
            "states": {"names": list(self.states.names), "angles": list(self.states.angles), "mirror": list(self.states.mirror)},
            "bonds_per_repeat": self.B,
            "first_order": self.first_order.tolist(),
            "second_order": self.second_order.tolist(),
        }
        if self.third_order is not None:
            d["third_order"] = self.third_order.tolist()
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RISModel":
        st = d["states"]
        states = RISStates(tuple(st["names"]), tuple(st["angles"]), tuple(st["mirror"]))
        e3 = np.array(d["third_order"]) if "third_order" in d else None
        return cls(states, int(d["bonds_per_repeat"]), np.array(d["first_order"]), np.array(d["second_order"]), e3, d.get("name", ""))

    def save(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str) -> "RISModel":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def describe(self) -> str:
        st = self.states
        lines = [f"RIS model {self.name!r}: {self.S} states " + ", ".join(f"{n}({a:.0f} deg)" for n, a in zip(st.names, st.angles)) + f"; {self.B} bond type(s) per repeat; order {self.order}"]
        for b in range(self.B):
            lines.append(f"  bond type {b}: first-order " + ", ".join(f"{n}={e:+.2f}" for n, e in zip(st.names, self.first_order[b])))
            lines.append(f"  pair (type {b} -> type {(b + 1) % self.B}) energies (kcal/mol):")
            lines.append("        " + "".join(f"{n:>8s}" for n in st.names))
            for s, n in enumerate(st.names):
                lines.append(f"    {n:>4s}" + "".join(f"{self.second_order[b, s, sp]:8.2f}" for sp in range(self.S)))
        if self.third_order is not None:
            lines.append("  third-order corrections (|e3| > 0.05 kcal/mol):")
            for b in range(self.B):
                for a in range(self.S):
                    for c in range(self.S):
                        for d in range(self.S):
                            v = self.third_order[b, a, c, d]
                            if abs(v) > 0.05:
                                lines.append(f"    type {b}: {st.names[a]}{st.names[c]}{st.names[d]} {v:+.2f}")
        return "\n".join(lines)


# ------------------------------------------------------- copolymer parameter transfer
#
# A copolymer's explicit repeat (``polymers.copolymer``) has many more bond types than the
# homopolymers it is made of, and only a few of them are chemically new.  Rather than
# annotate by hand which copolymer bond gets which homopolymer's fitted energies, the
# assignment below is *derived* from the backbone chemistry and reported, so the answer to
# "which bonds got which parameters, and which had none" is data rather than a claim.


def term_window(polymer, bond: int, order: int) -> tuple:
    """The backbone atoms an order-``order`` RIS term at bond type ``bond`` is a function of.

    An order-1 term (``first_order[b]``) is the energy of bond ``b``'s own dihedral, whose
    four atoms are ``b, b+1, b+2, b+3``.  An order-2 term (``second_order[b]``) spans bonds
    ``b`` and ``b+1``, so five atoms; order 3 spans three bonds and six.  Indices wrap
    modulo the repeat, which is correct because the chain is infinite and periodic.
    """
    B = polymer.bonds_per_repeat
    return tuple(polymer.backbone[(bond + i) % B] for i in range(order + 3))


def term_core(polymer, bond: int, order: int) -> tuple:
    """The atoms the term's rotatable bonds actually join: ``b+1 .. b+order+1``.

    For order 1 this is the two atoms of the rotating bond itself.  A torsional preference
    is set far more by those than by the 1-4 atoms at the window's ends, so the matching
    below ranks a candidate by its core agreement first and its full-window agreement
    second -- which is what makes the bonds flanking the VDCN carbon take VDCN's
    parameters rather than VDF's, despite agreeing with each in three atoms of four.
    """
    B = polymer.bonds_per_repeat
    return tuple(polymer.backbone[(bond + 1 + i) % B] for i in range(order + 1))


@dataclass(frozen=True)
class TermTransfer:
    """Where one RIS term of a copolymer got its parameters, and how good the match is."""

    order: int  # 1, 2 or 3
    bond: int  # bond type of the target repeat
    source: str  # source polymer name
    source_bond: int  # bond type within the source repeat
    core_matches: int
    core_size: int
    window_matches: int
    window_size: int

    @property
    def exact(self) -> bool:
        """True when the source term's whole backbone window is the same chemistry.

        An exact match means the transferred number really is the fitted number for this
        environment.  Anything else is a transfer: the right bond in the wrong
        neighbourhood, and not a fitted value for the bond it is being used for.
        """
        return self.window_matches == self.window_size

    def describe(self) -> str:
        kind = "fitted" if self.exact else "TRANSFERRED"
        return (f"order {self.order} bond {self.bond:2d} <- {self.source}[{self.source_bond}]  "
                f"core {self.core_matches}/{self.core_size} window {self.window_matches}/{self.window_size}  {kind}")


def _best_source(target, bond: int, order: int, sources) -> TermTransfer:
    core_t, win_t = term_core(target, bond, order), term_window(target, bond, order)
    best = None
    for poly, model in sources:
        for t in range(model.B):
            core_s, win_s = term_core(poly, t, order), term_window(poly, t, order)
            score = (sum(a == b for a, b in zip(core_t, core_s)), sum(a == b for a, b in zip(win_t, win_s)))
            if best is None or score > best[0]:
                best = (score, poly.name, t)
    (nc, nw), src, t = best
    return TermTransfer(order=order, bond=bond, source=src, source_bond=t, core_matches=nc,
                        core_size=len(core_t), window_matches=nw, window_size=len(win_t))


def transfer_ris(target, sources, name: str = "", third_order: bool | None = None) -> tuple[RISModel, list[TermTransfer]]:
    """Assemble a copolymer's RIS model from its comonomers' fitted homopolymer models.

    ``sources`` is a sequence of ``(Polymer, RISModel)`` pairs -- the homopolymer and the
    model fitted for it.  For every bond type of ``target`` and every order the model
    carries, the source term whose backbone environment best matches is used: core
    agreement (the atoms the rotating bonds join) first, then whole-window agreement, with
    ties going to the earlier source.  Returns the assembled model and one
    :class:`TermTransfer` per term, whose :attr:`TermTransfer.exact` flag separates the
    bonds that carry a genuinely fitted value for their own environment from those that
    carry a transferred one.

    This is the cheap half of ``docs/CHEMISTRY_EXTENSION.md`` section 3's copolymer bullet:
    it reuses the background fits and makes the junction terms visible so they can be
    fitted later, rather than pretending they were.
    """
    if not sources:
        raise ValueError("transfer_ris needs at least one (polymer, model) source")
    # The assembled model uses the TARGET's state set: a chain has one set of state
    # angles, and ``chain.build_chain`` reads them from the polymer, not per bond type.
    states = target.states
    for poly, model in sources:
        if model.states.names != states.names or model.states.mirror != states.mirror:
            raise ValueError(f"source model for {poly.name!r} has a different state set "
                             f"({model.states.names} vs {states.names})")
        if model.states.angles != states.angles:
            warnings.warn(
                f"source model for {poly.name!r} was fitted with adapted state angles "
                f"{tuple(round(a, 1) for a in model.states.angles)}, but the copolymer chain is built at "
                f"{tuple(round(a, 1) for a in states.angles)}: its energies are basin minima located at "
                "angles the copolymer will not use. Fit the sources with adapt_angles=False, or give the "
                "target a states= carrying the adapted angles, rather than mixing the two.",
                UserWarning, stacklevel=2,
            )
        if model.B != poly.bonds_per_repeat:
            raise ValueError(f"source model for {poly.name!r} has B={model.B} but the polymer has "
                             f"{poly.bonds_per_repeat} bonds per repeat")
    if third_order is None:
        third_order = all(m.third_order is not None for _, m in sources)
    elif third_order and any(m.third_order is None for _, m in sources):
        raise ValueError("third_order=True needs every source model to carry third-order terms")
    by_name = {poly.name: (poly, model) for poly, model in sources}
    B, S = target.bonds_per_repeat, states.n
    e1 = np.zeros((B, S))
    e2 = np.zeros((B, S, S))
    e3 = np.zeros((B, S, S, S)) if third_order else None
    report: list[TermTransfer] = []
    for order, dest in ((1, e1), (2, e2), (3, e3)):
        if dest is None:
            continue
        for b in range(B):
            tr = _best_source(target, b, order, sources)
            src = by_name[tr.source][1]
            arr = {1: src.first_order, 2: src.second_order, 3: src.third_order}[order]
            dest[b] = arr[tr.source_bond]
            report.append(tr)
    return RISModel(states, B, e1, e2, e3, name=name or f"{target.name}-transferred"), report


def transfer_summary(report: list[TermTransfer]) -> str:
    """A per-order table of :func:`transfer_ris`'s provenance, exact matches collapsed."""
    lines = []
    for order in sorted({t.order for t in report}):
        rows = [t for t in report if t.order == order]
        n_exact = sum(t.exact for t in rows)
        lines.append(f"order {order}: {n_exact}/{len(rows)} bond types carry a fitted value for their own environment")
        for t in rows:
            if not t.exact:
                lines.append("  " + t.describe())
    return "\n".join(lines)


def polyethylene_like_model(e_gauche: float = 0.5, e_pentane: float = 2.0, name="pe-textbook") -> RISModel:
    """Textbook three-state PE model: sigma = exp(-e_gauche/kT), omega = exp(-e_pentane/kT)."""
    S = 3
    e1 = np.array([[0.0, e_gauche, e_gauche]])
    e2 = np.zeros((1, S, S))
    e2[0, 1, 2] = e2[0, 2, 1] = e_pentane  # G+G- and G-G+ (pentane effect)
    return RISModel(THREE_STATE, 1, e1, e2, name=name)
