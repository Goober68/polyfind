"""Independent Cartesian atoms and six-strain relaxation of complete cells.

Reference lattice/atoms transform by a symmetric positive stretch exp(L).
L uses xx,yy,zz,2yz,2xz,2xy logarithmic strain. This clamps the macroscopic
polar rotation relative to the reference, not individual chain rotations.
Fixed log components are explicit constraints, not fixed lattice lengths
when free shears couple them. No chain-helix symmetry is imposed on atoms.
Physical Cauchy and log-conjugate stresses are reported separately; a solver
stop is not stationarity, stability, calibration or a physical trajectory.
"""
from dataclasses import dataclass
from enum import Enum

import numpy as np
from scipy.linalg import expm, expm_frechet
from scipy.optimize import minimize

from .mechanics import KCAL_MOL_A3_TO_GPA, strain_tensor
from .pack import MASS
from .periodic_geometry import ChainLayout, PlacedCell
from .topology import Cell

_BASIS = np.stack([strain_tensor(e) for e in np.eye(6)])
_BASIS.setflags(write=False)
KCAL_MOL_A3_TO_MPA = 1000.*KCAL_MOL_A3_TO_GPA


def _readonly(value):
    value = np.array(value,copy=True)
    value.setflags(write=False)
    return value


class TranslationGauge:
    """O(N) orthonormal Householder chart with mass-centre translation removed."""

    def __init__(self,masses):
        masses = np.asarray(masses)
        if (masses.ndim != 1 or not len(masses) or masses.dtype.kind not in "iuf"
                or not np.isfinite(masses).all() or np.any(masses <= 0)):
            raise ValueError("one positive finite real mass is required per atom")
        self.masses = _readonly(masses.astype(float))
        scaled = self.masses/self.masses.max()
        direction = -scaled/np.linalg.norm(scaled)
        direction[-1] += 1.
        norm = np.linalg.norm(direction)
        self.reflector = _readonly(direction/norm if norm else direction)

    def _reflect(self,value):
        return value-2.*self.reflector[:,None]*(self.reflector@value)[None,:]

    def displacements(self,coordinates):
        coordinates = np.asarray(coordinates,dtype=float)
        if coordinates.shape != (len(self.masses)-1,3) or not np.isfinite(coordinates).all():
            raise ValueError("translation-free coordinates must have shape (N-1,3)")
        return self._reflect(np.vstack([coordinates,np.zeros((1,3))]))

    def pullback(self,gradient):
        gradient = np.asarray(gradient,dtype=float)
        if gradient.shape != (len(self.masses),3) or not np.isfinite(gradient).all():
            raise ValueError("atomic gradient must have shape (N,3) and be finite")
        return self._reflect(gradient)[:-1]


class CartesianChart:
    """Own explicit strain constraints, translation gauge and exact pullback.

    Optimizer strain coordinates are eta*strain_scale; conditioning never
    changes energy or force/stress acceptance. A restart reuses this chart
    and accepted parameters, not an inferred canonical helix/setting angle.
    """

    def __init__(self,cell,free_strain=None,log_strain=None,strain_scale=None):
        self.layout = ChainLayout.from_cell(cell)
        self.reference = PlacedCell(cell.coords,cell.lattice)
        if len(self.reference.coords) != self.layout.atoms.size:
            raise ValueError("reference geometry must match the declared cell layout")
        free = np.ones(6,dtype=bool) if free_strain is None else np.asarray(free_strain)
        if free.shape != (6,) or free.dtype.kind != "b":
            raise ValueError("free log-strain mask must contain six booleans")
        eta = np.zeros(6) if log_strain is None else np.asarray(log_strain)
        if eta.shape != (6,) or eta.dtype.kind not in "iuf" or not np.isfinite(eta).all():
            raise ValueError("log strain must be a finite real six-vector")
        scale = np.asarray(len(cell.coords) if strain_scale is None else strain_scale)
        if scale.shape != () or scale.dtype.kind not in "iuf" or not np.isfinite(scale) or scale <= 0:
            raise ValueError("strain conditioning scale must be positive and finite")
        self.free_strain,self.log_strain = _readonly(free),_readonly(eta.astype(float))
        self.strain_scale = float(scale)
        try:
            self.gauge = TranslationGauge([MASS[e] for e in self.layout.elements])
        except KeyError as error:
            raise ValueError("declared atom element has no model mass") from error
        self.atomic_dof = 3*(len(cell.coords)-1)
        self.n_dof = self.atomic_dof+int(free.sum())
        self.x0 = _readonly(np.r_[np.zeros(self.atomic_dof),eta[free]*self.strain_scale])

    def decode(self,x):
        x = np.asarray(x)
        if x.shape != (self.n_dof,) or x.dtype.kind not in "iuf" or not np.isfinite(x).all():
            raise ValueError("Cartesian chart parameters must be a finite real vector")
        eta = self.log_strain.copy()
        eta[self.free_strain] = x[self.atomic_dof:]/self.strain_scale
        stretch = expm(strain_tensor(eta))
        undeformed = self.reference.coords+self.gauge.displacements(x[:self.atomic_dof].reshape(-1,3))
        placed = PlacedCell(undeformed@stretch,self.reference.lattice@stretch)
        return placed,undeformed,stretch,eta

    def cell(self,x):
        placed = self.decode(x)[0]
        nc,n = self.layout.atoms.shape
        chains,local = np.empty(nc*n,dtype=int),np.empty(nc*n,dtype=int)
        chains[self.layout.atoms] = np.arange(nc)[:,None]
        local[self.layout.atoms] = np.arange(n)[None,:]
        return Cell(list(self.layout.elements),placed.coords.copy(),placed.lattice.copy(),
                    chains,local,n,self.layout.reversed_of.copy())

    def pullback(self,x,evaluation):
        placed,undeformed,stretch,eta = self.decode(x)
        gP,gH = evaluation.grad_coords,evaluation.grad_lattice
        atomic = self.gauge.pullback(gP@stretch.T).ravel()
        gF = undeformed.T@gP+self.reference.lattice.T@gH
        log_gradient = np.array([np.sum(gF*expm_frechet(strain_tensor(eta),direction,compute_expm=False))
                                 for direction in _BASIS])
        gradient = np.r_[atomic,log_gradient[self.free_strain]/self.strain_scale]
        virial = placed.coords.T@gP+placed.lattice.T@gH
        stress = KCAL_MOL_A3_TO_MPA*np.einsum("ij,kij->k",virial,_BASIS)/np.linalg.det(placed.lattice)
        conjugate = KCAL_MOL_A3_TO_MPA*log_gradient/np.linalg.det(self.reference.lattice)
        return gradient,stress,conjugate


@dataclass(frozen=True)
class CartesianTolerance:
    atomic_force_kcal_mol_A: float = 1e-4
    stress_MPa: float = 1.

    def __post_init__(self):
        for value in (self.atomic_force_kcal_mol_A,self.stress_MPa):
            value = np.asarray(value)
            if value.shape != () or value.dtype.kind not in "iuf" or not np.isfinite(value) or value <= 0:
                raise ValueError("positive finite force and stress tolerances required")

    def assess(self,chart,x,evaluation):
        _,stress,conjugate = chart.pullback(x,evaluation)
        force = float(np.linalg.norm(evaluation.grad_coords,axis=1).max())
        free = chart.free_strain
        physical = float(np.max(abs(stress[free]))) if free.any() else 0.
        projected = float(np.max(abs(conjugate[free]))) if free.any() else 0.
        if not np.isfinite([force,physical,projected]).all():
            raise ValueError("finite force and full stress observations required")
        return CellStationarity(force,_readonly(stress),_readonly(conjugate),physical,projected,
                                bool(force <= self.atomic_force_kcal_mol_A and
                                     physical <= self.stress_MPa and projected <= self.stress_MPa))


@dataclass(frozen=True,eq=False)
class CellStationarity:
    max_atomic_force_kcal_mol_A: float
    stress_MPa: np.ndarray
    log_conjugate_stress_MPa: np.ndarray
    max_selected_physical_stress_MPa: float
    max_free_log_stress_MPa: float
    converged: bool


class RelaxationPhase(str,Enum):
    CREATED = "created"
    RUNNING = "running"
    CONVERGED = "converged"
    UNCONVERGED = "unconverged"
    FAILED = "failed"


@dataclass(frozen=True,eq=False)
class RelaxationCheckpoint:
    """Accepted optimizer state; a persistence owner can store/recover it.

    Parameters have meaning only with the SAME chart and model declaration.
    Source/protocol-authenticated durable campaign ownership is separate.
    """
    phase: RelaxationPhase
    parameters: np.ndarray
    accepted_steps: int
    energy_kcal_mol: float
    stationarity: CellStationarity


@dataclass(frozen=True,eq=False)
class CartesianRelaxation:
    checkpoint: RelaxationCheckpoint
    geometry: PlacedCell
    optimizer_status: int
    message: str
    evaluations: int


class CartesianRelaxer:
    """Single-run lifecycle with accepted checkpoints and terminal guards.

    CREATED -> RUNNING -> CONVERGED / UNCONVERGED / FAILED. An iteration
    limit or line-search stop never supplies a convergence claim. Checkpoint
    publication happens only at accepted states, not trial geometries.
    Recovery starts a new owner with the original chart and accepted x.
    No path trace is labelled molecular dynamics or a switching mechanism.
    """
    def __init__(self,packer,chart,tolerance=None,checkpoint_sink=None):
        self.packer,self.chart = packer,chart
        self.tolerance = tolerance or CartesianTolerance()
        self.checkpoint_sink = checkpoint_sink
        self.phase = RelaxationPhase.CREATED
        self.checkpoint = None
        self.evaluations = 0
        self.failure_diagnosis = None

    def _transition(self,phase):
        allowed = {RelaxationPhase.CREATED:{RelaxationPhase.RUNNING},
                   RelaxationPhase.RUNNING:{RelaxationPhase.CONVERGED,RelaxationPhase.UNCONVERGED,RelaxationPhase.FAILED}}
        if phase not in allowed.get(self.phase,set()):
            raise RuntimeError(f"illegal Cartesian relaxation transition {self.phase} -> {phase}")
        self.phase = phase

    def _evaluate(self,x):
        placed = self.chart.decode(x)[0]
        evaluation = self.packer.evaluate_chain_cell(placed.coords,placed.lattice,self.chart.layout)
        self.evaluations += 1
        gradient = self.chart.pullback(x,evaluation)[0]
        return evaluation,gradient

    def _accept(self,x,step):
        evaluation,_ = self._evaluate(x)
        stationarity = self.tolerance.assess(self.chart,x,evaluation)
        self.checkpoint = RelaxationCheckpoint(self.phase,_readonly(x),step,evaluation.terms.total,stationarity)
        if self.checkpoint_sink is not None:
            self.checkpoint_sink(self.checkpoint)
        return stationarity

    def run(self,maxiter=2000,initial=None):
        if isinstance(maxiter,bool) or not isinstance(maxiter,(int,np.integer)) or maxiter < 0:
            raise ValueError("iteration limit must be a nonnegative integer")
        self._transition(RelaxationPhase.RUNNING)
        x = self.chart.x0.copy() if initial is None else np.array(initial,copy=True)
        step = 0
        def fg(value):
            evaluation,gradient = self._evaluate(value)
            return evaluation.terms.total,gradient
        def accepted(value):
            nonlocal step
            step += 1
            if self._accept(value,step).converged:
                raise StopIteration
        try:
            initial_stationarity = self._accept(x,0)
            if initial_stationarity.converged or maxiter == 0 or self.chart.n_dof == 0:
                status,message = 0,"initial assessment" if initial_stationarity.converged else "iteration limit"
            else:
                result = minimize(fg,x,jac=True,method="L-BFGS-B",callback=accepted,
                                  options={"gtol":0.,"ftol":0.,"maxiter":int(maxiter),"maxcor":30})
                x = result.x
                status,message = int(result.status),str(result.message)
                self._accept(x,int(result.nit))
            terminal = RelaxationPhase.CONVERGED if self.checkpoint.stationarity.converged else RelaxationPhase.UNCONVERGED
            checkpoint = self.checkpoint
            completed = RelaxationCheckpoint(terminal,checkpoint.parameters,checkpoint.accepted_steps,
                                             checkpoint.energy_kcal_mol,checkpoint.stationarity)
            if self.checkpoint_sink is not None:
                self.checkpoint_sink(completed)
            self._transition(terminal)
            self.checkpoint = completed
            return CartesianRelaxation(self.checkpoint,self.chart.decode(x)[0],status,message,self.evaluations)
        except Exception as error:
            self.failure_diagnosis = f"{type(error).__name__}: {error}"
            if self.phase == RelaxationPhase.RUNNING:
                self._transition(RelaxationPhase.FAILED)
            if self.checkpoint is not None:
                checkpoint = self.checkpoint
                self.checkpoint = RelaxationCheckpoint(self.phase,checkpoint.parameters,checkpoint.accepted_steps,
                                                        checkpoint.energy_kcal_mol,checkpoint.stationarity)
            raise


def relax_cartesian_cell(packer,cell,*,free_strain=None,log_strain=None,tolerance=None,maxiter=2000,
                         checkpoint_sink=None):
    """Full independent-atom/six-strain model relaxation, or explicit clamps."""
    chart = CartesianChart(cell,free_strain,log_strain)
    return CartesianRelaxer(packer,chart,tolerance,checkpoint_sink).run(maxiter)
