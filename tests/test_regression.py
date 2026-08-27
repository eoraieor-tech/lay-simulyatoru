"""REQRESSİYA TESTİ — refaktorinqin nəticəni dəyişmədiyini qoruyur.

Etalon rəqəmlər refaktorinqdən ƏVVƏLKİ core.py-dən götürülüb.
Bu test pozulursa, hesablama davranışı dəyişib — dərhal araşdırılmalıdır.

Yavaşdır (~20 san). Keçmək üçün:  IMEX_SKIP_SLOW=1 pytest
"""

from helpers import (REFERENCE_FIVE_SPOT, SKIP_SLOW, default_scal,
                     five_spot_model, make_service)
from imex2d.application.config import SimulationConfig


def test_five_spot_reference_case_reproduces_legacy_results():
    if SKIP_SLOW:
        return
    scal = default_scal()
    model = five_spot_model(scal=scal)
    result = make_service(scal).run(model, SimulationConfig(end_time=1500.0))

    assert abs(result.ooip - REFERENCE_FIVE_SPOT["ooip"]) < 1.0, \
        f"OOIP dəyişib: {result.ooip:.1f}"
    assert abs(result.final_recovery_factor
               - REFERENCE_FIVE_SPOT["recovery_factor"]) < 0.01, \
        f"RF dəyişib: {result.final_recovery_factor:.3f} %"
    assert result.steps == REFERENCE_FIVE_SPOT["steps"], \
        f"Addım sayı dəyişib: {result.steps}"
    assert result.converged


def test_result_series_lengths_are_consistent():
    scal = default_scal()
    result = make_service(scal).run(five_spot_model(nx=15, ny=15, scal=scal),
                                    SimulationConfig(end_time=200.0))
    series = result.series
    n = len(series.time)
    for values in (series.oil_rate, series.water_rate, series.water_cut,
                   series.cumulative_oil, series.average_pressure,
                   series.recovery_factor):
        assert len(values) == n
    for rates in result.well_oil_rate.values():
        assert len(rates) == n
