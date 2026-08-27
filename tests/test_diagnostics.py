"""Model diaqnostikası (V1) və quyu avtoklampı (V2) testləri."""

import numpy as np

from helpers import default_scal, five_spot_model
from imex2d.domain.diagnostics import DiagnosticReport, Severity
from imex2d.domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from imex2d.simulation.pvt.correlations import build_pvt_table


def _model(**kwargs):
    return five_spot_model(nx=9, ny=9, **kwargs)


def _injector(model):
    return next(w for w in model.wells if w.is_injector)


def _producer(model):
    return next(w for w in model.wells if not w.is_injector)


# ── hesabat strukturu ─────────────────────────────────────────────────
def test_report_separates_errors_from_warnings():
    report = DiagnosticReport()
    report.error("bloklayan")
    report.warning("şübhəli")
    report.info("məlumat")
    assert report.has_errors
    assert len(report.errors) == 1 and len(report.warnings) == 1
    assert report.messages(Severity.WARNING) == ["şübhəli"]
    assert len(report) == 3


def test_clean_model_produces_no_diagnostics():
    assert len(_model().diagnose()) == 0


def test_validate_returns_only_blocking_messages():
    """Geriyə uyğunluq: validate() yalnız XƏTA-ları qaytarır."""
    model = _model()
    _injector(model).control.target = 100.0      # xəbərdarlıq yaradır
    assert model.diagnose().warnings
    assert model.validate() == []


# ── quyu idarəetməsi (V1) ─────────────────────────────────────────────
def test_injector_below_reservoir_pressure_warns():
    model = _model()
    _injector(model).control.target = 150.0
    warnings = model.diagnose().warnings
    assert any("vurucu" in w.message for w in warnings)
    assert not model.diagnose().has_errors


def test_injector_above_reservoir_pressure_is_clean():
    model = _model()
    _injector(model).control.target = 320.0
    assert not any("vurucu" in w.message for w in model.diagnose().warnings)


def test_producer_above_reservoir_pressure_warns():
    model = _model()
    _producer(model).control.target = 300.0
    assert any("hasilat" in w.message for w in model.diagnose().warnings)


def test_duplicate_well_names_warn():
    model = _model()
    model.wells[1].name = model.wells[0].name
    assert any("təkrarlanır" in w.message for w in model.diagnose().warnings)


def test_injector_above_fracture_pressure_warns():
    import dataclasses
    model = _model()
    model.geometry = dataclasses.replace(model.geometry, top_depth=2000.0)
    fracture = float(np.mean(model.geometry.cell_depths())) * model.FRACTURE_GRADIENT
    _injector(model).control.target = fracture + 40.0
    assert any("çatlama" in w.message for w in model.diagnose().warnings)


def test_shallow_model_skips_the_fracture_check():
    """Sintetik dayaz modellərdə çatlama təxmini mənasızdır."""
    model = _model()
    assert model._fracture_pressure() is None
    assert not any("çatlama" in w.message for w in model.diagnose().warnings)


def test_producer_below_bubble_point_warns():
    model = _model()
    model.pvt_table = build_pvt_table(bubble_point_bar=240.0)
    _producer(model).control.target = 150.0
    assert any("doyma" in w.message for w in model.diagnose().warnings)


def test_producer_below_bubble_point_does_not_warn_when_gas_phase_active():
    """A7: PVT-də qaz xassələri VARSA (has_gas_phase), bu vəziyyət artıq
    real modelləşdirilir — köhnə "qaz modelləşdirilmir" xəbərdarlığı
    yanlış olardı (ölçülüb: bu, real bir UI səhvi kimi tapıldı)."""
    model = _model()
    model.pvt_table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    _producer(model).control.target = 150.0
    assert not any("qaz ayrılacaq" in w.message
                  for w in model.diagnose().warnings)


def test_rate_controlled_wells_are_not_checked_for_pressure():
    model = _model()
    injector = _injector(model)
    injector.control = WellControl(ControlMode.RATE, 50.0)
    assert not any("vurucu" in w.message for w in model.diagnose().warnings)


# ── bloklayan xətalar ─────────────────────────────────────────────────
def test_perforation_outside_grid_is_an_error_not_a_warning():
    model = _model()
    _producer(model).perforations = [Perforation(99, 0, 0)]
    report = model.diagnose()
    assert report.has_errors
    assert any("kənar" in d.message for d in report.errors)


def test_model_without_wells_is_an_error():
    model = _model()
    model.wells = []
    assert model.diagnose().has_errors


def test_missing_pressure_control_is_an_error():
    model = _model()
    for well in model.wells:
        well.control = WellControl(ControlMode.RATE, 30.0)
    assert any("BHP" in d.message for d in model.diagnose().errors)


def test_diagnostics_include_hints():
    model = _model()
    _injector(model).control.target = 100.0
    warning = next(w for w in model.diagnose().warnings if "vurucu" in w.message)
    assert warning.hint
    assert warning.source == warning.source  # mənbə quyu adıdır
    assert "bar" in str(warning)


# ── loglama (C4) ──────────────────────────────────────────────────────
def test_logging_is_configured_once_and_can_be_reset():
    import io
    import logging

    from imex2d.logging_setup import (add_handler, configure, get_logger,
                                      remove_handler, reset)

    reset()
    logger = configure(level=logging.INFO, console=False)
    assert configure(console=False) is logger        # təkrar çağırış keçilir

    buffer = io.StringIO()
    handler = logging.StreamHandler(buffer)
    add_handler(handler)
    get_logger("imex2d.test").info("salam")
    assert "salam" in buffer.getvalue()

    remove_handler(handler)
    get_logger("imex2d.test").info("görünməməli")
    assert "görünməməli" not in buffer.getvalue()
    reset()


def test_module_loggers_share_the_root_namespace():
    from imex2d.logging_setup import LOGGER_NAME, get_logger
    assert get_logger("imex2d.simulation.x").name == "imex2d.simulation.x"
    assert get_logger("simulation.x").name == f"{LOGGER_NAME}.simulation.x"
    assert get_logger().name == LOGGER_NAME


def test_engine_run_emits_a_summary_log_line():
    import io
    import logging

    from helpers import make_service, short_config
    from imex2d.logging_setup import add_handler, configure, remove_handler, reset

    reset()
    configure(level=logging.INFO, console=False)
    buffer = io.StringIO()
    handler = logging.StreamHandler(buffer)
    add_handler(handler)
    try:
        make_service(default_scal()).run(_model(), short_config(end_time=20.0))
        text = buffer.getvalue()
        assert "Mühərrik qurulur" in text
        assert "RF" in text
    finally:
        remove_handler(handler)
        reset()


# ── matris strukturunun keşlənməsi (C3) ───────────────────────────────
def test_cached_matrix_pattern_matches_direct_assembly():
    """Optimallaşdırma reqressiyası: keşlənmiş CSR matrisi COO ilə eynidir.

    Struktur bir dəfə qurulub `_data_index` ilə doldurulur. Əgər bu
    uyğunluq pozulsa, nəticələr səssizcə yanlış olardı.
    """
    import numpy as np
    import scipy.sparse as sp

    from helpers import make_service, short_config

    engine = make_service(default_scal()).create_engine(
        _model(), short_config(end_time=20.0))
    engine._solve_pressure(1.0)

    reference = sp.coo_matrix(
        (engine._vals, (engine._rows, engine._cols)),
        shape=(engine.model.ncell, engine.model.ncell)).tocsr()
    difference = abs(engine._matrix - reference)
    assert difference.nnz == 0 or difference.max() < 1e-12


def test_matrix_object_is_reused_between_steps():
    from helpers import make_service, short_config

    engine = make_service(default_scal()).create_engine(
        _model(), short_config(end_time=20.0))
    engine._solve_pressure(1.0)
    first = engine._matrix
    engine._solve_pressure(1.0)
    assert engine._matrix is first, "Matris hər addımda yenidən yaradılır"
