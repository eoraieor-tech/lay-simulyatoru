"""Eclipse ixracı — quyu idarə sətirlərinin SÜTUNLARI (Seans 27).

TAPINTI: RATE istismarçısı `'LRAT' 2* debit 2* 1.0` kimi yazılırdı. WCONPROD
sütunları 4 ORAT · 5 WRAT · 6 GRAT · 7 LRAT · 8 RESV · 9 BHP olduğu üçün
debit 6-cı — QAZ debiti — sütununa düşürdü. Bundan başqa bizim RATE hədəfi
lay həcmidir, LRAT isə səth debitidir; uyğun rejim RESV-dir.

Qaz vuran quyu isə iki fazalı (OIL/WATER) deck-ə səssizcə SU vurucusu kimi
yazılırdı.
"""

from __future__ import annotations

import pytest

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.domain.wells import (ControlMode, Perforation, Phase, Well,
                                 WellControl, WellType)
from imex2d.io.eclipse_export import EclipseDeckWriter


def _model(producer=WellControl(ControlMode.RATE, 80.0),
           injector=WellControl(ControlMode.RATE, 90.0)):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    wells = [Well("INJ", WellType.INJECTOR, injector, [Perforation(0, 0, 0)]),
             Well("PROD", WellType.PRODUCER, producer, [Perforation(3, 3, 0)])]
    return ReservoirModelBuilder().build(geology, wells, scal=default_scal(),
                                         name="ixrac sınağı")


def _record(deck: str, keyword: str, well: str) -> list:
    """Açar sözün altındakı quyu sətri — `n*` defoltları AÇILMIŞ halda."""
    lines = deck.split(keyword + "\n", 1)[1].splitlines()
    line = next(l for l in lines if l.strip().startswith(f"'{well}'"))
    items = []
    for token in line.replace("/", " ").split():
        if token.endswith("*") and token[:-1].isdigit():
            items += ["*"] * int(token[:-1])
        else:
            items.append(token.strip("'"))
    return items


def test_rate_producer_target_lands_in_the_resv_column():
    items = _record(EclipseDeckWriter().render(_model()), "WCONPROD", "PROD")
    assert items[2] == "RESV"
    assert items[3:7] == ["*"] * 4, "ORAT/WRAT/GRAT/LRAT boş qalmalıdır"
    assert float(items[7]) == pytest.approx(80.0), "8-ci sütun RESV debitidir"
    assert float(items[8]) == pytest.approx(1.0), "9-cu sütun BHP həddidir"


def test_bhp_producer_line_is_unchanged():
    items = _record(EclipseDeckWriter().render(
        _model(producer=WellControl(ControlMode.BHP, 150.0))), "WCONPROD", "PROD")
    assert items[2] == "BHP" and float(items[8]) == pytest.approx(150.0)


def test_rate_injector_target_lands_in_the_resv_column():
    items = _record(EclipseDeckWriter().render(_model()), "WCONINJE", "INJ")
    assert items[1:4] == ["WATER", "OPEN", "RESV"]
    assert items[4] == "*", "5-ci sütun (səth debiti) boş qalmalıdır"
    assert float(items[5]) == pytest.approx(90.0)
    assert float(items[6]) == pytest.approx(1000.0)


def test_gas_injector_is_refused_instead_of_written_as_water():
    model = _model(injector=WellControl(ControlMode.RATE, 5000.0,
                                        injected_phase=Phase.GAS))
    with pytest.raises(ValueError, match="qaz"):
        EclipseDeckWriter().render(model)
