"""Force-label diagnosis preserves the source and the squared-error identity."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("reference_force_audit",Path(__file__).resolve().parents[1]/"examples/reference_force_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_translation_is_orthogonal_and_source_unchanged():
    forces = np.array([[1.,2.,3.],[-3.,2.,1.],[2.,-1.,-4.]])
    original = forces.copy()
    result = audit.translation_components(forces)
    assert result["squared"] == pytest.approx(result["translation_squared"]+result["internal_squared"])
    assert result["net_force_norm"] == pytest.approx(np.linalg.norm(forces.sum(axis=0)))
    np.testing.assert_array_equal(forces,original)


def test_zero_sum_model_cannot_remove_label_translation_error():
    labels = np.array([[1.,2.,3.],[-2.,3.,1.]])
    optimum = labels-labels.mean(axis=0)
    bound = audit.translation_components(labels)["translation_squared"]
    assert np.sum((optimum-labels)**2) == pytest.approx(bound)
    perturbation = np.array([[.2,-.3,.4],[-.2,.3,-.4]])
    assert np.sum((optimum+perturbation-labels)**2) == pytest.approx(bound+np.sum(perturbation**2))


def test_summary_weights_force_components_not_frames():
    records = [audit.translation_components(np.ones((2,3))),audit.translation_components(np.zeros((4,3)))]
    result = audit.translation_summary(records)
    assert result["component_rms_kcal_mol_A"] == pytest.approx(np.sqrt(1/3))
    assert result["translation_fraction_of_squared_norm"] == 1.
    assert result["orthogonal_squared_norm_closure_error"] == 0.
    assert audit.translation_summary([]) is None


@pytest.mark.parametrize("forces",[[],np.zeros((2,2)),np.array([[np.nan,0.,0.]]),np.array([[np.inf,0.,0.]])])
def test_invalid_force_arrays_rejected(forces):
    with pytest.raises(ValueError,match="N-by-3"):
        audit.translation_components(forces)
