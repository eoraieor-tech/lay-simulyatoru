"""Phase 2 — "Full Tensor Permeability Implementation" test suite.

Audit tapşırığı §15-də tələb olunan MINIMUM 12 test, dəqiq adlandırılmış
ssenarilərlə. Bəzi aspektlər (TPFA xəbərdarlığı) artıq `test_discretization.
py`-də əlavə örtülür — bura TƏKRAR salınmır, yalnız istinad edilir.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from imex2d.domain.properties import PermeabilityTensor, PropertyMap
from imex2d.domain.validation import ValidationResult


def _uniform_tensor(kxx, kyy, kzz, kxy=None, kxz=None, kyz=None, n=1):
    def _maybe(name, value):
        return None if value is None else PropertyMap.uniform(name, value, n)
    return PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", kxx, n),
        kyy=PropertyMap.uniform("KYY", kyy, n),
        kzz=PropertyMap.uniform("KZZ", kzz, n),
        kxy=_maybe("KXY", kxy), kxz=_maybe("KXZ", kxz), kyz=_maybe("KYZ", kyz))


# ── Test 1: izotrop ───────────────────────────────────────────────────────
def test_1_isotropic_tensor():
    tensor = _uniform_tensor(100.0, 100.0, 100.0)
    result = tensor.validate()
    assert result.ok
    assert np.allclose(np.sort(tensor.eigenvalues()[0]), [100.0, 100.0, 100.0])
    assert np.isclose(tensor.anisotropy_ratio()[0], 1.0)


# ── Test 2: diaqonal anizotropluq ──────────────────────────────────────────
def test_2_diagonal_anisotropy():
    tensor = _uniform_tensor(1000.0, 100.0, 10.0)
    assert tensor.validate().ok
    assert np.allclose(np.sort(tensor.eigenvalues()[0]), [10.0, 100.0, 1000.0])
    assert np.isclose(tensor.anisotropy_ratio()[0], 100.0)


# ── Test 3: fırlanmış anizotropluq ─────────────────────────────────────────
def _rotation_z(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def test_3_rotated_anisotropy_preserves_eigenvalues_and_positive_definiteness():
    original = _uniform_tensor(1000.0, 100.0, 10.0)
    original_eig = np.sort(original.eigenvalues()[0])

    rotated = original.rotate(_rotation_z(np.deg2rad(37.0)))
    rotated_eig = np.sort(rotated.eigenvalues()[0])

    assert np.allclose(rotated_eig, original_eig)
    assert rotated.validate().ok                       # müsbət-müəyyənlik qorunub
    assert rotated.has_off_diagonal()                   # off-diaqonal ORTAYA ÇIXDI

    # 1) izotrop tenzor izotrop QALIR
    iso = _uniform_tensor(100.0, 100.0, 100.0)
    iso_rotated = iso.rotate(_rotation_z(np.deg2rad(53.0)))
    assert not iso_rotated.has_off_diagonal()
    assert np.allclose(iso_rotated.as_matrices()[0], 100.0 * np.eye(3))

    # qeyri-ortoqonal "fırlanma" rədd edilməlidir
    with pytest.raises(ValueError):
        original.rotate(np.array([[2.0, 0, 0], [0, 1, 0], [0, 0, 1]]))


# ── Test 4: off-diaqonal tenzor (müsbət-müəyyən) ──────────────────────────
def test_4_off_diagonal_tensor_is_positive_definite():
    tensor = _uniform_tensor(100.0, 100.0, 50.0, kxy=20.0)
    result = tensor.validate()
    assert result.ok
    eig = np.sort(tensor.eigenvalues()[0])
    assert np.allclose(eig, [50.0, 80.0, 120.0])


# ── Test 5: etibarsız (indefinit) tenzor ──────────────────────────────────
def test_5_indefinite_tensor_fails_validation():
    tensor = _uniform_tensor(100.0, 100.0, 50.0, kxy=150.0)   # eig: -50, 50, 250
    result = tensor.validate()
    assert not result.ok
    assert any("müsbət-müəyyən" in e for e in result.errors)


# ── Test 6/7: NaN / Inf ────────────────────────────────────────────────────
def test_6_nan_fails_validation():
    tensor = _uniform_tensor(np.nan, 100.0, 50.0)
    result = tensor.validate()
    assert not result.ok


def test_7_inf_fails_validation():
    tensor = _uniform_tensor(np.inf, 100.0, 50.0)
    result = tensor.validate()
    assert not result.ok


# ── Test 8: güclü anizotropluq (süni kəsilmə OLMADAN) ─────────────────────
@pytest.mark.parametrize("ratio", [10.0, 100.0, 1000.0, 10000.0])
def test_8_strong_anisotropy_without_artificial_clipping(ratio):
    tensor = _uniform_tensor(ratio, 1.0, 1.0)
    result = tensor.validate()
    assert result.ok
    assert np.isclose(tensor.anisotropy_ratio()[0], ratio, rtol=1e-9)
    # dəyərlərin ÖZÜ dəyişməyib (heç bir "təmir"/kəsilmə yoxdur)
    assert tensor.kxx.values[0] == ratio


# ── Test 9: heterogen PropertyMap (hüceyrə-hüceyrə fərqli) ────────────────
def test_9_heterogeneous_property_map_validates_cell_by_cell():
    n = 4
    kxx = PropertyMap("KXX", np.array([100.0, 1000.0, 100.0, np.nan]))
    kyy = PropertyMap("KYY", np.array([100.0, 100.0, 100.0, 100.0]))
    kzz = PropertyMap("KZZ", np.array([100.0, 10.0, 50.0, 100.0]))
    kxy = PropertyMap("KXY", np.array([0.0, 0.0, 150.0, 0.0]))   # hüceyrə 2 indefinit
    tensor = PermeabilityTensor(kxx=kxx, kyy=kyy, kzz=kzz, kxy=kxy)

    assert tensor.ncell == n
    eig = tensor.eigenvalues()
    assert eig.shape == (n, 3)
    # hüceyrə 0: izotrop, müsbət-müəyyən
    assert np.all(eig[0] > 0)
    # hüceyrə 3: NaN → validate() bunu AYRICA tutmalıdır (aşağıda)
    result = tensor.validate()
    assert not result.ok
    assert any("NaN" in e or "sonsuz" in e for e in result.errors)


# ── Test 10: serializasiya round-trip ─────────────────────────────────────
def test_10_serialization_round_trip_preserves_every_component():
    from imex2d.application.serialization import (_permeability_tensor_from_dict,
                                                   _permeability_tensor_to_dict)
    n = 5
    rng = np.random.default_rng(0)
    tensor = PermeabilityTensor(
        kxx=PropertyMap("KXX", rng.uniform(50, 500, n), unit="mD"),
        kyy=PropertyMap("KYY", rng.uniform(50, 500, n), unit="mD"),
        kzz=PropertyMap("KZZ", rng.uniform(5, 50, n), unit="mD"),
        kxy=PropertyMap("KXY", rng.uniform(-50, 50, n), unit="mD"),
        kxz=PropertyMap("KXZ", rng.uniform(-10, 10, n), unit="mD"),
        kyz=PropertyMap("KYZ", rng.uniform(-10, 10, n), unit="mD"))

    payload = _permeability_tensor_to_dict(tensor)
    restored = _permeability_tensor_from_dict(payload)

    for name in ("kxx", "kyy", "kzz", "kxy", "kxz", "kyz"):
        original_component = getattr(tensor, name)
        restored_component = getattr(restored, name)
        assert np.allclose(restored_component.values, original_component.values)
        assert restored_component.unit == original_component.unit
    assert np.allclose(restored.as_matrices(), tensor.as_matrices())


def test_10b_full_project_round_trip_preserves_tensor():
    """Test 10-un tam layihə save/load səviyyəsində təkrarı."""
    import os
    import sys
    import tempfile

    sys.path.insert(0, os.path.dirname(__file__))
    from helpers import default_scal, five_spot_model
    from imex2d.application.project import Project
    from imex2d.application.serialization import ProjectSerializer

    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, scal=default_scal())
    n = model.ncell
    model.rock.permeability_tensor = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 200.0, n),
        kyy=PropertyMap.uniform("KYY", 150.0, n),
        kzz=PropertyMap.uniform("KZZ", 20.0, n),
        kxy=PropertyMap.uniform("KXY", 30.0, n))

    project = Project("Tensor testi")
    project.add_reservoir_model(model)

    serializer = ProjectSerializer()
    handle, path = tempfile.mkstemp(suffix=".imx")
    os.close(handle)
    try:
        serializer.save(project, path)
        restored_project = serializer.load(path)
    finally:
        os.unlink(path)

    restored_model = restored_project.reservoir_models[model.name]
    restored_tensor = restored_model.rock.permeability_tensor
    assert restored_tensor is not None
    assert np.allclose(restored_tensor.as_matrices(), model.rock.permeability_tensor.as_matrices())
    assert restored_tensor.validate().ok


# ── Test 11: geriyə uyğunluq (skalyar-K köhnə model) ──────────────────────
def test_11_backward_compatible_scalar_model_still_runs():
    import os
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from helpers import default_scal, five_spot_model, make_service, short_config

    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    assert model.rock.permeability_tensor is None   # köhnə model: tenzor YOXDUR

    result = make_service(scal).run(model, short_config(end_time=100.0))
    assert result.converged
    assert np.all(np.isfinite(result.series.average_pressure))


# ── Test 12: TPFA xəbərdarlığı (bax `test_discretization.py`-də əlavə) ────
def test_12_full_tensor_with_off_diagonal_triggers_explicit_tpfa_warning():
    import os
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from helpers import five_spot_model
    from imex2d.simulation.discretization import TwoPointFluxDiscretization

    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    n = model.ncell
    model.rock.permeability_tensor = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 200.0, n),
        kyy=PropertyMap.uniform("KYY", 200.0, n),
        kzz=PropertyMap.uniform("KZZ", 20.0, n),
        kxy=PropertyMap.uniform("KXY", 50.0, n))

    grid = TwoPointFluxDiscretization().build(model)
    assert any("MPFA" in w for w in grid.warnings)
    assert any("Kxy" in w for w in grid.warnings)   # YENİ: HANSI komponent olduğu göstərilir


# ── §7: vahid çevirməsi BÜTÜN 6 komponentə EYNİ tətbiq olunmalıdır ────────
def test_unit_conversion_applies_consistently_to_all_six_components():
    tensor = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 1.0, 2, unit="D"),
        kyy=PropertyMap.uniform("KYY", 2.0, 2, unit="D"),
        kzz=PropertyMap.uniform("KZZ", 0.1, 2, unit="D"),
        kxy=PropertyMap.uniform("KXY", 0.5, 2, unit="D"),
        kxz=PropertyMap.uniform("KXZ", 0.05, 2, unit="D"),
        kyz=PropertyMap.uniform("KYZ", 0.02, 2, unit="D"))

    converted = tensor.convert_units("D", "mD")   # 1 D = 1000 mD
    for name, original_value in (("kxx", 1.0), ("kyy", 2.0), ("kzz", 0.1),
                                 ("kxy", 0.5), ("kxz", 0.05), ("kyz", 0.02)):
        component = getattr(converted, name)
        assert np.allclose(component.values, original_value * 1000.0)
        assert component.unit == "mD"

    # diaqonal VƏ off-diaqonal AYNI faktorla çevrilib — nisbətlər DƏYİŞMƏYİB
    # (yəni eigenvalue/anisotropy strukturu invariant qalır)
    assert np.allclose(np.sort(tensor.eigenvalues()[0]) * 1000.0,
                       np.sort(converted.eigenvalues()[0]))


def test_property_map_rejects_invalid_unit_for_tensor_components():
    """`PROPERTY_QUANTITY`-ə KXY/KXZ/KYZ əlavə olunub — PERMX-lə EYNİ
    qorumaya malikdirlər (bax audit §7/§3)."""
    with pytest.raises(ValueError):
        PropertyMap("KXY", np.array([10.0]), unit="psi")   # təzyiq vahidi, keçiricilik YOX


# ── §17: performans — validasiya vektorlaşdırılıb, O(N) davranmalıdır ────
def test_validation_scales_linearly_not_quadratically():
    def _time_validate(n):
        rng = np.random.default_rng(1)
        tensor = PermeabilityTensor(
            kxx=PropertyMap("KXX", rng.uniform(50, 500, n)),
            kyy=PropertyMap("KYY", rng.uniform(50, 500, n)),
            kzz=PropertyMap("KZZ", rng.uniform(5, 50, n)),
            kxy=PropertyMap("KXY", rng.uniform(-20, 20, n)))
        start = time.perf_counter()
        tensor.validate()
        return time.perf_counter() - start

    _time_validate(1_000)          # isinmə (JIT/keş effektlərini aradan qaldırır)
    small = _time_validate(5_000)
    large = _time_validate(100_000)
    ratio = large / max(small, 1e-9)
    # 20x hüceyrə artımı; O(N²) olsaydı ~400x, O(N) isə ~20x gözlənilir —
    # 60x həddi hər ikisini aydın ayırd edir, ölçmə səs-küyünə dözümlüdür
    assert ratio < 60, f"Validasiya vaxtı N-dən SUPERXƏTTİ artır (nisbət={ratio:.1f})"
