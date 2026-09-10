"""B1 — MPFA-O-nun istifadəçi tərəfindən seçilə bilməsi.

Bu fayl MPFA-O-nun RİYAZİYYATINI təkrar yoxlamır — o, artıq
`tests/test_mpfa_o.py`, `tests/test_mpfa_o_global_assembly.py` və
`tests/test_phase_d_mpfa_integration.py`-də dərin doğrulanıb.

Burada yoxlanılan YEGANƏ şey çatışmayan halqadır (bax `ICRA_PLANI.md`
→ B1): konfiqurasiyadakı seçim → servis → mühərrikə ötürülən
diskretizasiya → `.imx` faylı → UI paneli.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal, five_spot_model
from imex2d.application.config import (FLUX_SCHEMES, MPFA_O, TPFA,
                                       SimulationConfig)
from imex2d.application.project import Project
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import (ModelValidationError,
                                                   SimulationService)
from imex2d.discretization import MPFAODiscretization
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


@pytest.fixture(scope="module")
def qt_app():
    """Real (ekransız) Qt tətbiqi — `test_ui_units.py`-dəki ilə eyni üsul."""
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
    yield QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _model(nx=4, ny=4):
    return five_spot_model(nx=nx, ny=ny)


def _service(engine_factory=FullyImplicitEngine):
    return SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(),
        engine_factory=engine_factory)


def _config(scheme=TPFA, end_time=10.0, **kwargs):
    return SimulationConfig(end_time=end_time, flux_scheme=scheme, **kwargs)


# ══════════════════════════════ konfiqurasiya ═══════════════════════════
def test_default_scheme_is_tpfa_so_existing_models_do_not_change():
    """Defolt DƏYİŞMİR — köhnə modellərin nəticəsi eyni qalmalıdır."""
    assert SimulationConfig().flux_scheme == TPFA
    assert not SimulationConfig().uses_multipoint_flux


def test_both_schemes_are_offered():
    assert FLUX_SCHEMES == (TPFA, MPFA_O)


def test_unknown_scheme_is_rejected_with_a_readable_message():
    issues = _config(scheme="MPFA-L").validate()
    assert issues
    assert "MPFA-L" in issues[0]
    assert MPFA_O in issues[0]        # nə seçmək olar — mesajda yazılır


def test_valid_schemes_pass_validation():
    for scheme in FLUX_SCHEMES:
        assert _config(scheme=scheme).validate() == []


# ══════════════════════════════ servis qatı ═════════════════════════════
def test_service_builds_tpfa_when_scheme_is_tpfa():
    engine = _service().create_engine(_model(), _config(TPFA))
    assert isinstance(engine.flux_discretization, TwoPointFluxDiscretization)
    assert not engine.flux_discretization.supports_multipoint_stencil()


def test_service_builds_mpfa_when_scheme_is_mpfa():
    """B1-in ƏSAS testi — bu, əvvəl MÜMKÜN DEYİLDİ.

    MPFA-O nüvəsi və onun Jakobiana qoşulması (Phase 5B-2) mövcud idi,
    amma servis `flux_discretization` ötürmədiyi üçün istifadəçi ona
    çata bilmirdi.
    """
    engine = _service().create_engine(_model(), _config(MPFA_O))
    assert isinstance(engine.flux_discretization, MPFAODiscretization)
    assert engine.flux_discretization.supports_multipoint_stencil()


def test_mpfa_engine_actually_uses_the_multipoint_residual_path():
    """Seçim mühərrikin İÇİNƏ qədər gedir — sadəcə saxlanmır."""
    engine = _service().create_engine(_model(), _config(MPFA_O))
    assert engine.residual_assembler._multipoint is True

    tpfa_engine = _service().create_engine(_model(), _config(TPFA))
    assert tpfa_engine.residual_assembler._multipoint is False


def test_impes_with_mpfa_gives_a_user_facing_message_not_a_crash():
    """IMPES çoxnöqtəli stensili qəsdən rədd edir (`_reject_multipoint_impes`).

    İstifadəçi `NotImplementedError` görməməlidir — nə etməli olduğunu
    deyən mesaj görməlidir.
    """
    with pytest.raises(ModelValidationError) as excinfo:
        _service(ImpesEngine).create_engine(_model(), _config(MPFA_O))

    message = str(excinfo.value)
    assert "MPFA-O" in message
    assert "implicit" in message.lower()      # həlli göstərilir
    assert "TPFA" in message                  # alternativ də göstərilir


def test_impes_with_tpfa_still_works_unchanged():
    engine = _service(ImpesEngine).create_engine(_model(), _config(TPFA))
    assert isinstance(engine.flux_discretization, TwoPointFluxDiscretization)


# ══════════════════════════════ `.imx` faylı ════════════════════════════
def _round_trip(config: SimulationConfig) -> SimulationConfig:
    serializer = ProjectSerializer()
    return serializer._config_from_dict(serializer._config_to_dict(config))


def test_scheme_survives_the_imx_round_trip():
    for scheme in FLUX_SCHEMES:
        assert _round_trip(_config(scheme)).flux_scheme == scheme


def test_old_imx_files_without_the_key_open_as_tpfa():
    """GERİYƏ UYĞUNLUQ — B1-dən əvvəlki fayllarda açar YOXDUR.

    Həmin layihələr TPFA ilə hesablanıb; açılanda da TPFA olmalıdır,
    yoxsa istifadəçinin köhnə nəticəsi səssizcə başqa ədədlərlə
    əvəzlənərdi.
    """
    serializer = ProjectSerializer()
    data = serializer._config_to_dict(_config(MPFA_O))
    del data["flux_scheme"]
    assert serializer._config_from_dict(data).flux_scheme == TPFA


def test_run_config_keeps_the_scheme_inside_a_project(tmp_path):
    project = Project("B1 sınağı")
    model = _model()
    project.add_reservoir_model(model)
    project.new_run(model.name, _config(MPFA_O))

    path = str(tmp_path / "b1.imx")
    ProjectSerializer().save(project, path)
    reopened = ProjectSerializer().load(path)
    assert reopened.latest_run().config.flux_scheme == MPFA_O


# ══════════════════════════════ UI paneli ═══════════════════════════════
def test_numerical_panel_offers_both_schemes_and_defaults_to_tpfa(qt_app):
    from imex2d.ui.panels import NumericalPanel

    panel = NumericalPanel()
    values = [panel.flux_scheme.itemData(i)
              for i in range(panel.flux_scheme.count())]
    assert values == list(FLUX_SCHEMES)
    assert panel.flux_scheme_choice() == TPFA
    assert panel.simulation_config().flux_scheme == TPFA


def test_numerical_panel_choice_reaches_the_config(qt_app):
    from imex2d.ui.panels import NumericalPanel

    panel = NumericalPanel()
    panel.set_flux_scheme(MPFA_O)
    assert panel.simulation_config().flux_scheme == MPFA_O


def test_unknown_scheme_falls_back_to_tpfa_in_the_panel(qt_app):
    """Köhnə/zədəli fayl paneli sındırmamalıdır."""
    from imex2d.ui.panels import NumericalPanel

    panel = NumericalPanel()
    panel.set_flux_scheme("MPFA-L")
    assert panel.flux_scheme_choice() == TPFA


# ══════════════════════════ uc-uca: sxem nəticəyə təsir edir ════════════
def test_choosing_mpfa_changes_results_under_anisotropy_end_to_end():
    """Seçim ƏDƏDLƏRƏ təsir edir — yəni həqiqətən işləyir.

    Fırladılmış (off-diagonal) tenzorda TPFA və MPFA-O EYNİ nəticəni
    verə BİLMƏZ; eyni çıxsaydı, seçim heç yerə çatmırdı deməkdir.
    """
    from imex2d.domain.properties import PermeabilityTensor, PropertyMap

    def anisotropic_model():
        model = _model(nx=5, ny=5)
        n = model.ncell
        angle = np.deg2rad(30.0)
        c, s = np.cos(angle), np.sin(angle)
        k1, k2 = 500.0, 50.0
        model.rock.permeability_tensor = PermeabilityTensor(
            kxx=PropertyMap("KXX", np.full(n, k1 * c * c + k2 * s * s)),
            kyy=PropertyMap("KYY", np.full(n, k1 * s * s + k2 * c * c)),
            kzz=PropertyMap("KZZ", np.full(n, 10.0)),
            kxy=PropertyMap("KXY", np.full(n, (k1 - k2) * c * s)))
        return model

    config = _config(TPFA, end_time=30.0)
    tpfa = _service().create_engine(anisotropic_model(), config).run()

    config_mpfa = _config(MPFA_O, end_time=30.0)
    mpfa = _service().create_engine(anisotropic_model(), config_mpfa).run()

    assert tpfa.series.cumulative_oil and mpfa.series.cumulative_oil
    assert not np.isclose(tpfa.series.cumulative_oil[-1],
                          mpfa.series.cumulative_oil[-1], rtol=1e-6)
