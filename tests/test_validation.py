"""Phase 1 — fiziki yoxlama qatı: sərt xəta vs. xəbərdarlıq ayrımı."""

from __future__ import annotations

import numpy as np

from imex2d.domain import validation as v


def test_porosity_rejects_negative_and_over_one():
    result = v.validate_porosity([0.1, -0.05, 0.2])
    assert not result.ok
    result = v.validate_porosity([0.1, 1.2, 0.2])
    assert not result.ok


def test_porosity_accepts_typical_range_without_warnings():
    result = v.validate_porosity([0.10, 0.18, 0.25, 0.30])
    assert result.ok
    assert result.warnings == []


def test_porosity_warns_but_does_not_reject_unusual_high_value():
    """0.42 fiziki cəhətdən mümkündür (boşluqlu qum) — RƏDD EDİLMİR."""
    result = v.validate_porosity([0.20, 0.42])
    assert result.ok
    assert result.warnings


def test_porosity_rejects_nan_and_inf():
    result = v.validate_porosity([0.2, float("nan"), 0.3])
    assert not result.ok
    assert "NaN" in result.errors[0]
    result = v.validate_porosity([0.2, float("inf")])
    assert not result.ok


def test_saturation_bounds():
    assert v.validate_saturation([0.0, 0.5, 1.0]).ok
    assert not v.validate_saturation([-0.01, 0.5]).ok
    assert not v.validate_saturation([0.5, 1.5]).ok


def test_permeability_rejects_zero_and_negative():
    assert not v.validate_permeability([100.0, 0.0]).ok
    assert not v.validate_permeability([100.0, -5.0]).ok


def test_permeability_warns_on_extreme_but_valid_values():
    high = v.validate_permeability([50000.0])
    assert high.ok and high.warnings
    low = v.validate_permeability([0.001])
    assert low.ok and low.warnings
    typical = v.validate_permeability([10.0, 150.0, 500.0])
    assert typical.ok and not typical.warnings


def test_viscosity_rejects_non_positive():
    assert not v.validate_viscosity(0.0).ok
    assert not v.validate_viscosity(-1.0).ok
    assert v.validate_viscosity(3.2).ok


def test_viscosity_warns_on_heavy_oil_range_but_accepts():
    result = v.validate_viscosity(80000.0)
    assert result.ok and result.warnings


def test_density_rejects_non_positive_and_warns_on_extreme():
    assert not v.validate_density(-10.0).ok
    normal = v.validate_density(850.0)
    assert normal.ok and not normal.warnings
    extreme = v.validate_density(50.0)
    assert extreme.ok and extreme.warnings


def test_compressibility_rejects_non_positive():
    assert not v.validate_compressibility(0.0).ok
    assert not v.validate_compressibility(-1e-5).ok
    assert v.validate_compressibility(4.5e-5).ok


def test_pressure_rejects_non_positive():
    assert not v.validate_pressure([0.0, 100.0]).ok
    assert not v.validate_pressure([-5.0]).ok
    assert v.validate_pressure([200.0]).ok


def test_pressure_warns_above_fracture_gradient_estimate_but_does_not_reject():
    # 2500 m dərinlikdə defolt qradiyentlə (0.160 bar/m) hədd = 400 bar
    result = v.validate_pressure([450.0], depth_m=2500.0)
    assert result.ok
    assert result.warnings
    result_ok = v.validate_pressure([350.0], depth_m=2500.0)
    assert result_ok.ok and not result_ok.warnings


def test_grid_dimensions_reject_non_positive():
    assert not v.validate_grid_dimensions(0, 10, 1, 20.0, 20.0).ok
    assert not v.validate_grid_dimensions(10, 10, 1, -5.0, 20.0).ok
    assert v.validate_grid_dimensions(41, 41, 3, 20.0, 20.0).ok


def test_thickness_rejects_non_positive():
    assert not v.validate_thickness([10.0, 0.0, 5.0]).ok
    assert not v.validate_thickness([10.0, -5.0]).ok
    assert v.validate_thickness([10.0, 8.0, 12.0]).ok


def test_cell_volumes_rejects_zero_and_negative_and_warns_on_extreme():
    assert not v.validate_cell_volumes([100.0, 0.0, 50.0]).ok
    assert not v.validate_cell_volumes([100.0, -1.0]).ok
    tiny = v.validate_cell_volumes([1e-6, 100.0])
    assert tiny.ok and tiny.warnings
    huge = v.validate_cell_volumes([1e12, 100.0])
    assert huge.ok and huge.warnings
    normal = v.validate_cell_volumes([100.0, 200.0, 4000.0])
    assert normal.ok and not normal.warnings


def test_well_rate_rejects_negative_warns_on_zero_and_extreme():
    assert not v.validate_well_rate(-10.0).ok
    zero = v.validate_well_rate(0.0)
    assert zero.ok and zero.warnings
    huge = v.validate_well_rate(500000.0)
    assert huge.ok and huge.warnings
    normal = v.validate_well_rate(60.0)
    assert normal.ok and not normal.warnings


def test_validate_query_range_accepts_in_range_values():
    result = v.validate_query_range([150.0, 200.0], 1.0, 400.0, "PVT")
    assert result.ok and result.warnings == []


def test_validate_query_range_warns_on_mild_extrapolation():
    # diapazon eni 399, severe_factor=0.5 -> 199.5-ə qədər kənar YÜNGÜLDÜR
    result = v.validate_query_range([450.0], 1.0, 400.0, "PVT")
    assert result.ok
    assert result.warnings


def test_validate_query_range_hard_errors_on_severe_extrapolation():
    """3000 (məs. bar sanılan psi dəyəri) 1-400 diapazonundan HƏDDİNDƏN
    ARTIQ kənardadır — bu, VAHİD QARIŞIQLIĞI əlamətidir, sərt xətadır."""
    result = v.validate_query_range([3000.0], 1.0, 400.0, "PVT")
    assert not result.ok
    assert "QARIŞIQLIĞI" in result.errors[0]


def test_check_extrapolation_range_reports_out_of_bounds_queries():
    warnings = v.check_extrapolation_range([50.0, 150.0, 500.0], table_min=100.0,
                                           table_max=400.0, label="PVT təzyiqi")
    assert len(warnings) == 2   # 50 (aşağı) və 500 (yuxarı)
    assert v.check_extrapolation_range([150.0, 200.0], 100.0, 400.0) == []
