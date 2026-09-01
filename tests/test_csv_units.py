"""Phase 1 (giriş boru xətti) — CSV sütun-adı vahid metadatası (`NAME[vahid]`).

Format: `well_data_io.py` modul docstring-i. Qısaca: `PERMX[D]` sütunu
`PERMX` kimi saxlanılır, amma dəyərlər mühərrik vahidinə (mD) çevrilir.
Vahid göstərilməyəndə (köhnə format) DƏYƏR DƏYİŞMİR — yalnız
keçiricilik üçün AÇIQ xəbərdarlıq əlavə olunur (`dataset.warnings`).

Qeyd: `WellDataset.validate()` hər xassə üçün ən azı 2 nöqtə tələb edir
(mövcud, bu fazadan əvvəlki qayda) — ona görə bütün nümunə fayllarda
ən azı 2 sətir var.
"""

from __future__ import annotations

import os
import tempfile

import pytest

from imex2d.geology.well_data_io import read_well_csv


def _write(text: str) -> str:
    handle, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


def test_permx_in_darcy_converts_to_engine_md():
    path = _write("well,x,y,PORO,PERMX[D]\nW-1,10,20,0.20,0.15\nW-2,90,80,0.25,0.40\n")
    try:
        dataset = read_well_csv(path)
        values = [s.values["PERMX"] for s in dataset.samples]
        assert values == pytest.approx([150.0, 400.0])   # 1 D = 1000 mD
    finally:
        os.unlink(path)


def test_permx_in_md_bracket_is_explicit_and_unchanged():
    path = _write("well,x,y,PORO,PERMX[mD]\nW-1,10,20,0.20,150\nW-2,90,80,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.samples[0].values["PERMX"] == 150.0
        assert dataset.warnings == []   # vahid AÇIQ göstərilib — xəbərdarlıq yoxdur
    finally:
        os.unlink(path)


def test_legacy_permx_without_unit_is_unchanged_but_warns():
    """GERİYƏ UYĞUNLUQ: köhnə format (vahidsiz) DƏYƏRİ DƏYİŞMİR."""
    path = _write("well,x,y,PORO,PERMX\nW-1,10,20,0.20,150\nW-2,90,80,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.samples[0].values["PERMX"] == 150.0   # DƏYİŞMƏYİB
        assert dataset.samples[1].values["PERMX"] == 400.0
        assert any("PERMX" in w and "mD" in w for w in dataset.warnings)
    finally:
        os.unlink(path)


def test_legacy_poro_column_never_warns():
    """PORO ölçüsüz kəmiyyətdir — vahidsiz olanda xəbərdarlıq YOXDUR
    (bax modul docstring-i: yalnız keçiricilik üçün xəbərdarlıq)."""
    path = _write("well,x,y,PORO,PERMX[mD]\nW-1,10,20,0.20,150\nW-2,90,80,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.warnings == []
    finally:
        os.unlink(path)


def test_xy_coordinates_in_feet_convert_to_engine_metres():
    from imex2d.domain import unit_conversions as uc
    path = _write("well,x[ft],y[ft],PORO,PERMX\n"
                  "W-1,328.084,656.168,0.20,150\nW-2,984.252,328.084,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        sample = dataset.samples[0]
        assert sample.x == pytest.approx(uc.ft_to_m(328.084))
        assert sample.y == pytest.approx(uc.ft_to_m(656.168))
        assert sample.x == pytest.approx(100.0, abs=1e-3)
    finally:
        os.unlink(path)


def test_unrecognized_unit_on_dimensionless_property_warns_and_keeps_raw_value():
    path = _write("well,x,y,PORO[%],PERMX\nW-1,10,20,20.0,150\nW-2,90,80,25.0,400\n")
    try:
        dataset = read_well_csv(path)
        # PORO üçün heç bir kəmiyyət-çevirmə tanınmır -> dəyər DƏYİŞMİR,
        # amma AÇIQ bildirilir ki, çevrilmədiyi bilinsin
        assert dataset.samples[0].values["PORO"] == 20.0
        assert any("PORO" in w for w in dataset.warnings)
    finally:
        os.unlink(path)


def test_legacy_files_without_any_bracket_syntax_are_completely_unaffected():
    """Mövcud (bracket-siz) fayllar üçün YALNIZ keçiricilik xəbərdarlığı
    ƏLAVƏ OLUNUR — heç bir ƏDƏDİ DƏYƏR dəyişmir."""
    path = _write("well,x,y,PORO,PERMX,NTG\nW-1,10,20,0.20,150,0.9\nW-2,90,80,0.25,400,0.8\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.samples[0].x == 10.0 and dataset.samples[0].y == 20.0
        assert dataset.samples[0].values == {"PORO": 0.20, "PERMX": 150.0, "NTG": 0.9}
        assert dataset.samples[1].values == {"PORO": 0.25, "PERMX": 400.0, "NTG": 0.8}
    finally:
        os.unlink(path)
