"""SPE1CASE2 modeli (Seans 38) — qurucu, vahidlər, ilkin vəziyyət, qısa qaçış.

OPM deck faylı repoda deyil; burada onun `PROPS` cədvəlləri (SPE1 benchmark
məlumatı — Odeh 1981) yazılıb ki, qurucu faylsız sınansın. Etalonla tam
müqayisə 80+ saniyədir və real OPM faylları tələb edir —
`tools/spe1_compare.py` ilə aparılır.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.benchmarks.spe1 import (REFERENCE_BLOCKS, build_spe1case2_model,
                                    first_time_above, first_time_below,
                                    gas_front_arrivals, simulated_series,
                                    spe1case2_config)
from imex2d.domain.unit_conversions import convert
from imex2d.domain.wells import ControlMode, Phase, RateBasis
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

PROPS = """PROPS
PVTW
    4017.55 1.038 3.22E-6 0.318 0.0 /

ROCK
    14.7 3E-6 /

SWOF
0.12  0                      1       0
0.18  4.64876033057851E-008  1       0
0.24  0.000000186            0.997   0
0.3   4.18388429752066E-007  0.98    0
0.36  7.43801652892562E-007  0.7     0
0.42  1.16219008264463E-006  0.35    0
0.48  1.67355371900826E-006  0.2     0
0.54  2.27789256198347E-006  0.09    0
0.6   2.97520661157025E-006  0.021   0
0.66  3.7654958677686E-006   0.01    0
0.72  4.64876033057851E-006  0.001   0
0.78  0.000005625            0.0001  0
0.84  6.69421487603306E-006  0       0
0.91  8.05914256198347E-006  0       0
1     0.00001                0       0 /

SGOF
0     0      1      0
0.001 0      1      0
0.02  0      0.997  0
0.05  0.005  0.980  0
0.12  0.025  0.700  0
0.2   0.075  0.350  0
0.25  0.125  0.200  0
0.3   0.190  0.090  0
0.4   0.410  0.021  0
0.45  0.60   0.010  0
0.5   0.72   0.001  0
0.6   0.87   0.0001 0
0.7   0.94   0.000  0
0.85  0.98   0.000  0
0.88  0.984  0.000  0 /

DENSITY
    53.66 64.49 0.0533 /

PVDG
14.700  166.666  0.008000
264.70  12.0930  0.009600
514.70  6.27400  0.011200
1014.7  3.19700  0.014000
2014.7  1.61400  0.018900
2514.7  1.29400  0.020800
3014.7  1.08000  0.022800
4014.7  0.81100  0.026800
5014.7  0.64900  0.030900
9014.7  0.38600  0.047000 /

PVTO
0.0010  14.7    1.0620  1.0400 /
0.0905  264.7   1.1500  0.9750 /
0.1800  514.7   1.2070  0.9100 /
0.3710  1014.7  1.2950  0.8300 /
0.6360  2014.7  1.4350  0.6950 /
0.7750  2514.7  1.5000  0.6410 /
0.9300  3014.7  1.5650  0.5940 /
1.2700  4014.7  1.6950  0.5100
        9014.7  1.5790  0.7400 /
1.6180  5014.7  1.8270  0.4490
        9014.7  1.7370  0.6310 /
/

SOLUTION
"""


@pytest.fixture(scope="module")
def deck_path():
    handle, path = tempfile.mkstemp(suffix=".DATA")
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(PROPS)
    yield path
    os.unlink(path)


@pytest.fixture(scope="module")
def model(deck_path):
    return build_spe1case2_model(deck_path)


def _engine(model, end_time):
    config = spe1case2_config()
    config.end_time = end_time
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(model.scal_parameters),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    return service.create_engine(model, config)


# ═══════════════════════ quruluş və vahidlər ═════════════════════════

def test_grid_and_layer_properties(model):
    assert (model.grid.nx, model.grid.ny, model.grid.nz) == (10, 10, 3)
    assert model.geometry.dz == pytest.approx([6.096, 9.144, 15.24])
    layers = model.rock.permx.values.reshape(3, 100)
    assert layers[:, 0] == pytest.approx([500.0, 50.0, 200.0])
    assert np.array_equal(model.rock.permz.values, model.rock.permx.values)
    assert model.rock.porosity.values == pytest.approx(0.3)


def test_rock_and_fluid_constants_are_converted(model):
    assert model.rock.compressibility == pytest.approx(
        convert(3e-6, "bar", "psi", "pressure"), rel=1e-9)  # 1/psi → 1/bar
    assert model.rock.compressibility_reference_pressure == pytest.approx(1.01353,
                                                                          rel=1e-5)
    assert model.fluids.oil_density == pytest.approx(859.55, rel=1e-5)
    assert model.fluids.gas_density == pytest.approx(0.85378, rel=1e-4)


def test_wells_follow_the_deck(model):
    producer, injector = model.wells
    assert (producer.perforations[0].i, producer.perforations[0].j,
            producer.perforations[0].k) == (9, 9, 2)
    assert producer.control.mode is ControlMode.RATE
    assert producer.control.rate_basis is RateBasis.SURFACE
    assert producer.control.target == pytest.approx(3179.75, rel=1e-5)
    assert producer.control.bhp_limit == pytest.approx(68.948, rel=1e-4)
    assert injector.control.injected_phase is Phase.GAS
    assert injector.control.target == pytest.approx(2.8317e6, rel=1e-4)
    assert injector.control.bhp_limit == pytest.approx(621.47, rel=1e-4)
    assert producer.radius == pytest.approx(0.0762)


def test_initial_conditions_and_tables(model):
    ic = model.initial_conditions
    assert ic.use_equilibration
    assert ic.datum_pressure == pytest.approx(330.948, rel=1e-5)
    assert ic.solution_gor == pytest.approx(226.197, rel=1e-5)
    assert ic.water_saturation == pytest.approx(0.12)
    assert model.scal_parameters.sor == pytest.approx(0.16)
    assert len(model.pvt_oil_branches) == 9
    assert model.gas_scal_tables is not None and model.scal_tables is not None


def test_model_passes_diagnostics_without_errors(model):
    assert model.validate() == []


# ═══════════════════════ ilkin vəziyyət və qısa qaçış ════════════════

def test_initial_state_is_undersaturated_at_the_datum_pressure(model):
    engine = _engine(model, 10.0)
    state = engine.state
    assert not np.any(state.is_saturated), "SPE1-də ilkin sərbəst qaz yoxdur"
    assert state.third_variable == pytest.approx(226.197, rel=1e-5)
    layer_mean = [float(np.mean(state.pressure[k * 100:(k + 1) * 100]))
                  for k in range(3)]
    # datum (8400 ft) 3-cü layın mərkəzidir; yuxarı laylar hidrostatik az
    assert convert(layer_mean[2], "bar", "psi", "pressure") == pytest.approx(4800.0,
                                                                             rel=1e-6)
    assert layer_mean[0] < layer_mean[1] < layer_mean[2]


def test_short_run_holds_the_oil_rate_target(model):
    engine = _engine(model, 30.0)
    result = engine.run()
    assert result.converged, result.message
    series = simulated_series(model, result)
    time, oil, unit = series["FOPR"]
    assert unit == "STB/DAY"
    assert oil[-1] == pytest.approx(20000.0, rel=1e-3)
    assert "WBHP:INJ" in series and "BPR:300" in series


def test_first_time_below():
    assert first_time_below(np.array([1.0, 2.0, 3.0]), np.array([5.0, 4.0, 1.0]),
                            3.0) == 3.0
    assert first_time_below(np.array([1.0]), np.array([5.0]), 3.0) is None


# ═══════════════ qaz cəbhəsi — etalonun BGSAT blokları ════════════════

def test_reference_blocks_match_eclipse_ordering(model):
    """Etalonun blok nömrələri bizim hüceyrə indeksi ilə EYNİ sıralamadadır.

    Eclipse təbii sıralaması `i`-ni ən sürətli dəyişir. Bizim
    `CartesianGrid.index` də belədir, ona görə blok N ↔ hüceyrə N−1.
    Bu yoxlama olmasa BGSAT müqayisəsi səssizcə BAŞQA hüceyrəni oxuyar
    və nəticə yanlış olduğu bilinməzdi.
    """
    grid = model.grid
    for block in REFERENCE_BLOCKS:
        assert 1 <= block <= grid.ncell
        i, j, k = grid.ijk(block - 1)
        assert grid.index(i, j, k) == block - 1
    assert grid.ijk(0) == (0, 0, 0), "blok 1 — vurucunun hüceyrəsi"
    assert grid.ijk(299) == (9, 9, 2), "blok 300 — istismarçının hüceyrəsi"


def test_first_time_above():
    assert first_time_above(np.array([1.0, 2.0, 3.0]),
                            np.array([0.0, 0.0, 0.5]), 0.01) == 3.0
    assert first_time_above(np.array([1.0]), np.array([0.0]), 0.01) is None


def test_gas_saturation_blocks_are_reported(model):
    """BGSAT seriyaları çıxarılır, fiziki hədlərdədir və məkanca doğrudur.

    30 gündə vurucu blokuna (1) vurulan qazın həcmi həmin blokun məsamə
    həcmindən böyükdür, ona görə orada qaz OLMALIDIR.

    İstismarçının blokunda (300) isə 30 gündə KİÇİK miqdarda qaz olur və
    bu, vurulan qazın çatması DEYİL — 9 hüceyrə uzaqdan çata bilməz.
    Səbəb başqadır: təzyiq müvəqqəti olaraq doyma nöqtəsindən (4014.7
    psia) aşağı düşür və həll olmuş qaz ayrılır. ETALON DA eynisini
    göstərir — `SPE1CASE2.UNSMRY`-dən ölçüldü:

        gün      31      59     120     181     212
        BGSAT  0.0131  0.0222  0.0196  0.0054  0.0000
        BPR     3934    3869    3905    4012    4079  psia

    yəni təzyiq bərpa olunanda qaz GERİ HƏLL OLUR. Ona görə hədd
    "sıfır" deyil, "kiçik"dir.
    """
    engine = _engine(model, 30.0)
    result = engine.run()
    assert result.converged, result.message
    series = simulated_series(model, result)
    for block in REFERENCE_BLOCKS:
        label = f"BGSAT:{block}"
        assert label in series, label
        values = series[label][1]
        assert np.all(values >= -1e-12) and np.all(values <= 1.0 + 1e-12)
    assert series["BGSAT:1"][1][-1] > 0.01, "vurucu blokunda qaz olmalıdır"
    assert 0.0 < series["BGSAT:300"][1][-1] < 0.05, "kiçik, keçici ayrılma"


def test_gas_front_arrivals_skips_missing_labels(model):
    """Etalonda olmayan etiket sükutla buraxılır (süni etalonla yoxlanılır)."""

    class _Reference:
        def __contains__(self, label):
            return label in ("TIME", "BGSAT:1")

        def series(self, label):
            return (np.array([0.0, 10.0, 20.0]) if label == "TIME"
                    else np.array([0.0, 0.0, 0.4]))

    engine = _engine(model, 30.0)
    result = engine.run()
    rows = gas_front_arrivals(model, result, _Reference())
    assert [row.block for row in rows] == [1]
    assert rows[0].reference == 20.0
    assert rows[0].ijk == (0, 0, 0)
