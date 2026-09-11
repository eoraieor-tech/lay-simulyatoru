"""B4 (A variantı) — quyu başı təzyiqi (THP) və şaquli axın traversi.

ƏHATƏ (sahibkarın seçimi, `ISH_HESABATI.md` → Seans 11):
  * BHP məlumdur → yuxarı traverse → THP (mühərrik TOXUNULMUR);
  * sürüşmə YOXDUR (no-slip), interfeys Beggs-Brill üçün açıqdır;
  * yalnız BHP rejimli istismarçılar;
  * lülə tam şaquli (TVD = MD).

Testlərin məntiqi: traversin hər komponentini ANALİTİK həddlə
tutuşdurmaq. Sürtünməsiz hədddə nəticə hidrostatik sütuna, sıfır
uzunluqda isə BHP-nin özünə bərabər olmalıdır — bunlar uydurma
etalon deyil, düsturun öz sərhəd halıdır.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.properties import FluidProperties
from imex2d.domain.scal import GasCoreyParameters
from imex2d.domain.tubing import TubingGeometry
from imex2d.domain.wells import ControlMode, Well, WellControl, WellType
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.wellbore import (NoSlipHoldup, WellStream,
                                        chen_friction_factor,
                                        pressure_traverse, reynolds)
from imex2d.simulation.wellbore.hydraulics import WellboreHydraulics

GRAVITY = 9.80665


def _fluids():
    return FluidProperties(water_fvf=1.0, oil_fvf=1.0)


def _gas_pvt(bubble_point=240.0):
    return BlackOilPVTProvider(build_pvt_table(
        pressure_min=1.0, pressure_max=400.0, n_points=40,
        bubble_point_bar=bubble_point, include_gas=True))


# ═══════════════════════ sürtünmə əmsalı ═══════════════════════════════

def test_laminar_branch_is_the_analytic_value():
    """Re < 2000 → `f = 64/Re`, yaxınlaşma deyil, DƏQİQ düstur."""
    for re in (100.0, 500.0, 1500.0):
        assert chen_friction_factor(re, 0.001) == pytest.approx(64.0 / re)


def test_chen_matches_colebrook_within_half_a_percent():
    """Chen (1979) öz məqaləsində Colebrook-dan < 0.5 % fərq iddia edir.

    Colebrook burada İTERASİYA ilə həll olunur — yəni etalon müstəqildir.
    """
    def colebrook(re, eps):
        f = 0.02
        for _ in range(200):
            rhs = -2.0 * math.log10(eps / 3.7 + 2.51 / (re * math.sqrt(f)))
            f = 1.0 / rhs ** 2
        return f

    for re in (5e3, 1e4, 1e5, 1e6, 1e7):
        for eps in (0.0, 1e-5, 1e-4, 1e-3, 1e-2):
            chen = chen_friction_factor(re, eps)
            exact = colebrook(re, eps)
            assert abs(chen - exact) / exact < 0.005, (re, eps, chen, exact)


def test_friction_factor_is_zero_without_flow():
    """Dayanmış mayedə sürtünmə itkisi yoxdur — `64/Re` partlamamalıdır."""
    assert chen_friction_factor(0.0, 0.001) == 0.0


def test_transition_band_has_no_jump():
    """2000 → 4000 aralığında sıçrayış olmamalıdır.

    Sıçrayış Nyutonun yığılmasını pozardı (B variantında `f` qalığa
    girəcək), ona görə keçid xəttidir.
    """
    values = [chen_friction_factor(re, 1e-4)
              for re in np.linspace(1900, 4100, 60)]
    steps = np.abs(np.diff(values))
    assert steps.max() < 0.02, f"ən böyük sıçrayış {steps.max():.4f}"


def test_reynolds_converts_centipoise():
    """μ cP-dədir; 1 cP = 1e-3 Pa·s."""
    assert reynolds(1000.0, 1.0, 0.1, 1.0) == pytest.approx(1000.0 * 0.1 / 1e-3)


# ═══════════════════════ traverse — analitik həddlər ═══════════════════

def test_frictionless_limit_is_the_hydrostatic_column():
    """ƏSAS SƏRHƏD HALI: sürtünmə sıfıra gedəndə nəticə `ρ·g·h`.

    Diametr böyük seçilir ki, sürət (deməli sürtünmə) sıfıra yaxınlaşsın.
    """
    tubing = TubingGeometry(diameter=2.0, roughness=0.0, segments=40)
    stream = WellStream(water=100.0, water_density=1000.0)
    result = pressure_traverse(bhp=250.0, perforation_depth=2000.0,
                               tubing=tubing, stream=stream,
                               pvt=None, fluids=_fluids())

    analytic = 1000.0 * GRAVITY * 2000.0 / 1e5
    assert result.gravity_drop == pytest.approx(analytic, rel=1e-9)
    assert result.friction_drop < 1e-6
    assert result.thp == pytest.approx(250.0 - analytic, rel=1e-9)


def test_zero_length_returns_the_bottom_hole_pressure():
    """Perforasiya quyu başındadırsa düşgü yoxdur."""
    tubing = TubingGeometry(wellhead_depth=1500.0)
    stream = WellStream(water=100.0)
    result = pressure_traverse(200.0, 1500.0, tubing, stream, None, _fluids())
    assert result.thp == pytest.approx(200.0)


def test_friction_grows_with_rate():
    """Monotonluq — debit artdıqca sürtünmə itkisi artmalıdır."""
    tubing = TubingGeometry(diameter=0.062, segments=20)
    drops = []
    for rate in (50.0, 200.0, 500.0, 1000.0):
        stream = WellStream(water=rate, water_density=1000.0)
        drops.append(pressure_traverse(250.0, 2000.0, tubing, stream,
                                       None, _fluids()).friction_drop)
    assert all(b > a for a, b in zip(drops, drops[1:])), drops


def test_narrower_tubing_loses_more_pressure():
    """Eyni debitdə dar boru daha çox sürtünmə itkisi verir."""
    stream = WellStream(water=300.0, water_density=1000.0)
    wide = pressure_traverse(250.0, 2000.0, TubingGeometry(diameter=0.100),
                             stream, None, _fluids())
    narrow = pressure_traverse(250.0, 2000.0, TubingGeometry(diameter=0.050),
                               stream, None, _fluids())
    assert narrow.friction_drop > wide.friction_drop


def test_segments_converge():
    """Seqment sayı artdıqca nəticə SABİTLƏŞMƏLİDİR.

    Defolt 20 seqmentin kifayət etdiyini bu test kilidləyir: 20 ilə
    200 arasındakı fərq 0.05 bar-dan kiçikdir.
    """
    stream = WellStream(oil=80.0, water=20.0, gas=80.0 * 150.0)
    thp = {}
    for n in (20, 200):
        tubing = TubingGeometry(diameter=0.062, segments=n)
        thp[n] = pressure_traverse(250.0, 2000.0, tubing, stream,
                                   _gas_pvt(), _fluids()).thp
    assert abs(thp[20] - thp[200]) < 0.05, thp


def test_single_segment_is_measurably_worse():
    """Tək seqmentin YETƏRSİZ olduğunu sənədləşdirir.

    Bu test planın ilkin eskizinə (tək seqment) qarşı verilmiş qərarın
    ölçülmüş əsasıdır — bax `QARARLAR.md` → Q-10.
    """
    stream = WellStream(oil=80.0, water=20.0, gas=80.0 * 150.0)
    one = pressure_traverse(250.0, 2000.0, TubingGeometry(segments=1),
                            stream, _gas_pvt(), _fluids()).thp
    many = pressure_traverse(250.0, 2000.0, TubingGeometry(segments=200),
                             stream, _gas_pvt(), _fluids()).thp
    assert abs(one - many) > 1.0, (one, many)


# ═══════════════════════ qazın təsiri ═════════════════════════════════

def test_more_gas_lightens_the_column():
    """FİZİKİ MONOTONLUQ: GOR artdıqca hidrostatik sütun yüngülləşir.

    Bu test bir SƏHVİ tutmuşdu: Bo doymuş cədvəldən oxunanda (neftin
    həqiqi Rs-i nəzərə alınmadan) nəticə qeyri-monoton çıxırdı —
    B3-B-dəki səhvin lülədəki analoqu. Bax `traverse._segment_gradient`.
    """
    pvt = _gas_pvt()
    gravity = []
    for gor in (0.0, 50.0, 150.0, 400.0):
        stream = WellStream(oil=80.0, water=20.0, gas=80.0 * gor)
        gravity.append(pressure_traverse(
            250.0, 2000.0, TubingGeometry(segments=20), stream,
            pvt, _fluids()).gravity_drop)
    assert all(b < a for a, b in zip(gravity, gravity[1:])), gravity


def test_dissolved_gas_uses_the_undersaturated_branch():
    """Doymamış neftdə Bo doymuş qoldan OXUNMAMALIDIR.

    GOR doyma həddindən aşağı olanda sərbəst qaz yoxdur, lakin neft
    yüngülləşir — sütun çəkisi GOR = 0 halından AZ olmalıdır.
    """
    pvt = _gas_pvt()
    dead = pressure_traverse(250.0, 2000.0, TubingGeometry(segments=20),
                             WellStream(oil=100.0, gas=0.0), pvt, _fluids())
    live = pressure_traverse(250.0, 2000.0, TubingGeometry(segments=20),
                             WellStream(oil=100.0, gas=100.0 * 30.0),
                             pvt, _fluids())
    assert live.gravity_drop < dead.gravity_drop


# ═══════════════════════ sərhəd halları ═══════════════════════════════

def test_shut_well_returns_nan_not_a_made_up_number():
    """Dayanmış quyuda axan traverse TƏYİN OLUNMAYIB."""
    result = pressure_traverse(250.0, 2000.0, TubingGeometry(),
                               WellStream(), None, _fluids())
    assert math.isnan(result.thp)
    assert not result.converged
    assert "axmır" in result.message


def test_well_that_cannot_reach_surface_returns_nan():
    """Sütun BHP-ni üstələyirsə `nan` — sıfıra "qısaldılmış" dəyər YOX.

    Sıfır qrafikdə real ölçmə kimi görünərdi; `nan` isə boşluq buraxır
    və mesajda səbəb yazılır.
    """
    stream = WellStream(water=10.0, water_density=1000.0)
    result = pressure_traverse(bhp=50.0, perforation_depth=3000.0,
                               tubing=TubingGeometry(segments=20),
                               stream=stream, pvt=None, fluids=_fluids())
    assert math.isnan(result.thp)
    assert not result.converged
    assert "süni qaldırma" in result.message


def test_no_slip_holdup_equals_the_volume_fraction():
    """V1-in açıq fərziyyəsi: `H_L = λ_L`."""
    correlation = NoSlipHoldup()
    for lam in (0.0, 0.25, 0.5, 1.0):
        assert correlation.liquid_holdup(lam, 1.0, 1.0, 0.062) == pytest.approx(lam)


def test_tubing_validation_catches_bad_geometry():
    assert TubingGeometry(diameter=-1.0).validate()
    assert TubingGeometry(diameter=0.05, roughness=0.06).validate()
    assert TubingGeometry(segments=0).validate()
    assert not TubingGeometry().validate()


def test_few_segments_warns():
    assert any("seqment" in w for w in TubingGeometry(segments=2).validate_warnings())


# ═══════════════════════ uc-uca: simulyasiya ilə ══════════════════════

def _model(top_depth=1200.0, with_tubing=True, mode=ControlMode.BHP):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0, porosity=0.22,
        permx_base=150.0, top_depth=top_depth)
    wells = five_spot(geology.grid)
    for well in wells:
        if well.well_type is WellType.PRODUCER:
            well.control = WellControl(mode, well.control.target)
            if with_tubing:
                well.tubing = TubingGeometry(diameter=0.062, segments=20)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=default_scal(),
        gas_scal=GasCoreyParameters(),
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=True),
        name="THP sınağı")


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)


def test_simulation_records_thp_and_bhp():
    """ƏSAS TƏLƏB: qaçışdan sonra hər iki təzyiq sırası mövcuddur."""
    result = _service().run(_model(), SimulationConfig(end_time=300.0))
    assert result.converged, result.message

    assert result.well_thp, "THP sırası yazılmadı"
    assert result.well_bhp, "BHP sırası yazılmadı"
    for name, thp in result.well_thp.items():
        assert len(thp) == len(result.well_bhp[name])
        finite = [v for v in thp if math.isfinite(v)]
        assert finite, f"{name}: bütün THP dəyərləri nan"
        assert all(v < result.well_bhp[name][0] for v in finite), \
            "THP quyu dibi təzyiqindən kiçik olmalıdır"


def test_no_tubing_means_no_thp_and_no_behaviour_change():
    """GERİYƏ UYĞUNLUQ: lülə verilməyəndə heç nə dəyişmir."""
    with_tubing = _service().run(_model(with_tubing=True),
                                 SimulationConfig(end_time=300.0))
    without = _service().run(_model(with_tubing=False),
                             SimulationConfig(end_time=300.0))

    assert not without.well_thp and not without.well_bhp
    # Mühərrik toxunulmadığı üçün nəticə BİTƏ-BİT eyni olmalıdır.
    assert without.series.recovery_factor[-1] == pytest.approx(
        with_tubing.series.recovery_factor[-1], rel=0, abs=0)


def test_rate_controlled_wells_are_skipped_in_v1():
    """V1 məhdudiyyəti: RATE rejimində mühərrik BHP-ni hesablamır."""
    result = _service().run(_model(mode=ControlMode.RATE),
                            SimulationConfig(end_time=200.0))
    assert not result.well_thp


def test_injectors_are_skipped_in_v1():
    model = _model()
    hydraulics = WellboreHydraulics()
    injectors = [w for w in model.wells if w.well_type is WellType.INJECTOR]
    assert injectors
    for well in injectors:
        well.tubing = TubingGeometry()
        assert hydraulics._skip_reason(well, {well.name: 1000.0}) is not None


def test_very_deep_well_cannot_flow_at_all():
    """5000 m-də BHP 150 bar sütunu heç vaxt qaldıra bilmir → hamısı `nan`.

    Ölçüldü (dərinliyə görə axan addımların sayı, 24 addımdan):
    1200 m → 24, 2000 m → 20, 3000 m → 11, 3500 m → 1, 5000 m → 0.
    Yəni davranış kəskin deyil, TƏDRİCİDİR: qaz ayrıldıqca sütun
    yüngülləşir və bəzi addımlarda quyu axır, bəzilərində yox.
    """
    result = _service().run(_model(top_depth=5000.0),
                            SimulationConfig(end_time=200.0))
    assert result.converged, result.message
    assert result.well_thp
    for thp in result.well_thp.values():
        assert all(math.isnan(v) for v in thp)


def test_intermediate_depth_flows_only_part_of_the_time():
    """Sərhəd rejimi: bəzi addımlar axır, bəziləri yox — ikisi də qanunidir."""
    result = _service().run(_model(top_depth=3000.0),
                            SimulationConfig(end_time=200.0))
    assert result.converged, result.message
    for thp in result.well_thp.values():
        finite = [v for v in thp if math.isfinite(v)]
        assert 0 < len(finite) < len(thp), (len(finite), len(thp))


def test_tubing_survives_the_imx_round_trip(tmp_path):
    original = TubingGeometry(diameter=0.08, roughness=1e-4,
                              wellhead_depth=25.0, segments=33)
    data = ProjectSerializer._well_to_dict(Well("P1", tubing=original))
    restored = ProjectSerializer._well_from_dict(data)
    assert restored.tubing == original


def test_old_imx_without_the_key_has_no_tubing():
    """B4-dən əvvəlki fayllar açılmalı və THP-siz qalmalıdır."""
    data = ProjectSerializer._well_to_dict(Well("P1"))
    data.pop("tubing", None)
    assert ProjectSerializer._well_from_dict(data).tubing is None
