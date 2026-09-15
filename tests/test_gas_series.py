"""Sahə qaz seriyaları — GOR və vurulan qaz (B7 addım 1-in qalıqları).

ÖLÇÜLMÜŞ SƏHV: üç fazalı mühərrik sahə qaz debitini `-min(rates.gas.sum(), 0)`
ilə hesablayırdı. Hüceyrə cəmi vurulan (+) və hasil olunan (−) qazı
qarışdırırdı; vurucu güclü olanda cəm müsbət çıxır və `min` onu sıfırlayırdı.
Qaz vuran modeldə 1500 gün boyu sahə qaz debiti və GOR = 0 idi, halbuki
istismarçının öz GOR-u (həll olmuş qaz daxil) 109 sm³/sm³-dən başlayırdı.

Vurulan qaz isə heç bir seriyada saxlanılmırdı — qrafikdə yalnız «Vurulan su»
var idi.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.figure import Figure

from test_gas_injection import _model, _run
from imex2d.domain.wells import Phase
from imex2d.rendering import renderers as R
from imex2d.reporting.results_export import _columns


@pytest.fixture(scope="module")
def gas_result():
    result = _run(_model(Phase.GAS))
    assert result.converged, result.message
    return result


def test_field_gas_rate_is_the_producers_gas_not_net_of_injection(gas_result):
    series = gas_result.series
    producer = np.asarray(gas_result.well_gas_rate["PROD"])
    assert np.allclose(series.gas_rate, producer, rtol=1e-12, atol=0.0)
    assert min(series.gas_rate) > 0.0, "istismarçı qaz verir — sahə debiti sıfır ola bilməz"


def test_gor_includes_dissolved_gas_while_gas_is_injected(gas_result):
    """Sərbəst qaz gəlməmişdən ƏVVƏL də GOR > 0: neftdə həll olmuş qaz (Rs)."""
    series = gas_result.series
    gor = np.asarray(series.gas_oil_ratio)
    expected = np.asarray(gas_result.well_gas_rate["PROD"]) / np.asarray(
        gas_result.well_oil_rate["PROD"])
    assert gor[0] > 10.0, gor[:3]
    assert np.allclose(gor, expected, rtol=1e-9)


def test_injected_gas_is_recorded_as_a_series(gas_result):
    series = gas_result.series
    assert len(series.gas_injection_rate) == len(series.time)
    assert min(series.gas_injection_rate) > 0.0
    assert not any(series.water_injection_rate), "qaz vurucusu su vurmur"


def test_water_injection_leaves_the_gas_injection_series_at_zero():
    result = _run(_model(Phase.WATER))
    assert result.converged, result.message
    assert result.series.gas_injection_rate
    assert not any(result.series.gas_injection_rate)
    assert min(result.series.gas_rate) > 0.0


def test_two_phase_result_has_no_gas_injection_series():
    result = _run(_model(Phase.WATER, include_gas=False))
    assert result.converged, result.message
    assert result.series.gas_injection_rate == []


def test_export_writes_the_injected_gas_column(gas_result):
    headers = [header for header, _ in _columns(gas_result)]
    assert "q_vurulan_qaz [m³/gün]" in headers


def test_rate_panel_shows_the_injected_gas(gas_result):
    figure = Figure()
    R.ProductionCurveRenderer().draw(figure.subplots(2, 2), gas_result)
    labels = {line.get_label() for ax in figure.axes for line in ax.get_lines()}
    assert {"Vurulan qaz", "Qaz", "Neft"} <= labels, labels
