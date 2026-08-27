"""PDF hesabat generatoru (B6)."""

import os
import tempfile

import matplotlib
matplotlib.use("Agg")

import numpy as np
from pypdf import PdfReader

from helpers import default_scal, make_service, short_config
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.observations import (ObservationSet, ObservedQuantity,
                                        ObservedSeries)
from imex2d.domain.structure import FaultReference
from imex2d.history.mismatch import MismatchCalculator
from imex2d.reporting.report import (ReportContext, ReportGenerator,
                                     ReportSections)


def _geology(nx=10, ny=8, nz=1):
    return SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.21,
        permx_base=180.0, nz=nz, top_depth=2000.0)


def _model(faults=None):
    geology = _geology()
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        fault_references=faults)


def _path() -> str:
    handle, path = tempfile.mkstemp(suffix=".pdf")
    os.close(handle)
    os.unlink(path)          # ReportGenerator özü yaratmalıdır
    return path


def _page_count(path: str) -> int:
    return len(PdfReader(path).pages)


# ── minimal hesabat (model tək başına) ─────────────────────────────────
def test_report_can_be_built_from_a_model_alone():
    path = _path()
    try:
        ReportGenerator().write(ReportContext(model=_model()), path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 1000
        assert _page_count(path) >= 2          # başlıq + xülasə minimum
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_report_pdf_has_correct_metadata():
    path = _path()
    try:
        model = _model()
        ReportGenerator().write(ReportContext(model=model, author="Sınaq"),
                                path)
        reader = PdfReader(path)
        assert model.name in (reader.metadata.title or "")
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_summary_page_lists_wells_and_faults():
    model = _model(faults=[FaultReference(
        name="F1", source_id="F1", axis="I", plane_index=3,
        transmissibility_multiplier=0.2)])
    lines = ReportGenerator._summary_lines(model)
    text = "\n".join(lines)
    assert "INJ" in text or "PROD" in text
    assert "F1" in text


def test_summary_page_reports_no_faults_when_empty():
    lines = ReportGenerator._summary_lines(_model(faults=[]))
    assert "Fault yoxdur." in "\n".join(lines)


# ── nəticə əlavə olunanda ────────────────────────────────────────────
def test_report_grows_with_simulation_result():
    scal = default_scal()
    without_result = _path()
    with_result = _path()
    try:
        ReportGenerator().write(ReportContext(model=_model()), without_result)
        model = _model()
        result = make_service(scal).run(model, short_config(end_time=200.0))
        ReportGenerator().write(ReportContext(model=model, result=result),
                                with_result)
        assert _page_count(with_result) > _page_count(without_result)
    finally:
        for path in (without_result, with_result):
            if os.path.exists(path):
                os.unlink(path)


def test_title_page_includes_the_recovery_factor_when_available():
    scal = default_scal()
    model = _model()
    result = make_service(scal).run(model, short_config(end_time=200.0))
    path = _path()
    try:
        ReportGenerator().write(ReportContext(model=model, result=result), path)
        text = PdfReader(path).pages[0].extract_text()
        assert "Recovery Factor" in text
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_history_match_page_is_skipped_without_observations():
    path = _path()
    try:
        ReportGenerator().write(ReportContext(model=_model(), mismatch=None),
                                path)
        # heç bir xəta atılmamalıdır — bölmə sadəcə keçilir
        assert os.path.exists(path)
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_history_match_page_appears_with_observations():
    scal = default_scal()
    model = _model()
    result = make_service(scal).run(model, short_config(end_time=300.0))
    time = np.asarray(result.series.time, float)
    sample_times = np.linspace(time[1], time[-1] * 0.9, 8)
    observations = ObservationSet(series=[ObservedSeries(
        "", ObservedQuantity.OIL_RATE, sample_times,
        np.interp(sample_times, time, result.series.oil_rate))])
    mismatch = MismatchCalculator().evaluate(result, observations)

    with_match = _path()
    without_match = _path()
    try:
        ReportGenerator().write(
            ReportContext(model=model, result=result, mismatch=mismatch),
            with_match)
        ReportGenerator().write(
            ReportContext(model=model, result=result, mismatch=None),
            without_match)
        assert _page_count(with_match) > _page_count(without_match)
    finally:
        for path in (with_match, without_match):
            if os.path.exists(path):
                os.unlink(path)


# ── bölmə seçimi ────────────────────────────────────────────────────
def test_sections_can_be_disabled():
    minimal = ReportSections(summary=False, diagnostics=False, maps=False,
                             scal=False, pvt=False, results=False,
                             history_match=False)
    path = _path()
    try:
        ReportGenerator().write(
            ReportContext(model=_model(), sections=minimal), path)
        assert _page_count(path) == 1          # yalnız başlıq
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_maps_section_can_be_disabled_independently():
    sections = ReportSections(diagnostics=False, scal=False, pvt=False,
                              results=False, history_match=False)
    with_maps = _path()
    sections_no_maps = ReportSections(diagnostics=False, maps=False,
                                      scal=False, pvt=False, results=False,
                                      history_match=False)
    without_maps = _path()
    try:
        ReportGenerator().write(
            ReportContext(model=_model(), sections=sections), with_maps)
        ReportGenerator().write(
            ReportContext(model=_model(), sections=sections_no_maps),
            without_maps)
        assert _page_count(with_maps) > _page_count(without_maps)
    finally:
        for path in (with_maps, without_maps):
            if os.path.exists(path):
                os.unlink(path)


# ── açıq tema tətbiqi ──────────────────────────────────────────────────
def test_print_friendly_forces_white_background_and_dark_text():
    from matplotlib.figure import Figure

    from imex2d.rendering.renderers import ScalRenderer
    from imex2d.domain.scal import CoreyParameters

    figure = Figure()
    axes = figure.subplots(1, 2)
    ScalRenderer().draw(axes, CoreyParameters(), 0.5, 3.0)

    ReportGenerator._print_friendly(figure)
    for ax in figure.get_axes():
        assert ax.get_facecolor() == (1.0, 1.0, 1.0, 1.0)
        if ax.get_title():
            assert ax.title.get_color() == "#111827"


def test_output_directory_is_created_if_missing():
    directory = tempfile.mkdtemp()
    path = os.path.join(directory, "nested", "report.pdf")
    try:
        ReportGenerator().write(ReportContext(model=_model()), path)
        assert os.path.exists(path)
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_diagnostics_page_appears_only_when_there_are_issues():
    from imex2d.domain.structure import FaultReference

    clean_path = _path()
    broken_path = _path()
    try:
        ReportGenerator().write(ReportContext(model=_model()), clean_path)
        broken = _model(faults=[FaultReference(
            name="BAD", source_id="BAD", axis="I", plane_index=999)])
        ReportGenerator().write(ReportContext(model=broken), broken_path)
        assert _page_count(broken_path) > _page_count(clean_path)
    finally:
        for path in (clean_path, broken_path):
            if os.path.exists(path):
                os.unlink(path)
