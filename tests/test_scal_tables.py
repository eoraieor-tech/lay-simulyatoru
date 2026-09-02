"""Cədvəl əsaslı SCAL və region dəstəyi (B4 + B7)."""

import os
import tempfile

import numpy as np

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.domain.scal import CoreyParameters
from imex2d.domain.scal_tables import SaturationTable, SaturationTableSet
from imex2d.io.scal_io import (ScalFormatError, read_scal_csv, read_swof,
                               write_scal_csv)
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.scal_tables_provider import (
    TableCapillaryPressureProvider, TableRelativePermeabilityProvider)


def _table(**kwargs) -> SaturationTable:
    return SaturationTable.from_corey(CoreyParameters(**kwargs), points=41)


def _set(**kwargs) -> SaturationTableSet:
    tables = SaturationTableSet()
    tables.add(1, _table(**kwargs))
    return tables


def _two_regions() -> SaturationTableSet:
    tables = SaturationTableSet()
    tables.add(1, SaturationTable.from_corey(
        CoreyParameters(swc=0.15, sor=0.20, krw_end=0.5, nw=2.0), name="qum"))
    tables.add(2, SaturationTable.from_corey(
        CoreyParameters(swc=0.30, sor=0.35, krw_end=0.2, nw=3.5),
        name="gilli"))
    return tables


def _write(text: str, suffix=".csv") -> str:
    handle, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


# ── cədvəl obyekti ────────────────────────────────────────────────────
def test_endpoints_are_derived_from_the_table():
    table = SaturationTable([0.20, 0.30, 0.60, 0.80],
                            [0.00, 0.02, 0.18, 0.35],
                            [0.80, 0.55, 0.10, 0.00])
    assert abs(table.swc - 0.20) < 1e-9
    assert abs(table.sor - 0.20) < 1e-9        # kro son dəfə 0.60-da müsbət
    assert abs(table.krw_end - 0.35) < 1e-9


def test_corey_conversion_reproduces_the_formula():
    scal = CoreyParameters()
    table = SaturationTable.from_corey(scal, points=81)
    sw = np.linspace(scal.swc + 0.01, 1.0 - scal.sor - 0.01, 25)
    assert np.max(np.abs(table.interpolate_krw(sw) - scal.krw(sw))) < 1e-3
    assert np.max(np.abs(table.interpolate_kro(sw) - scal.kro(sw))) < 1e-3


def test_validation_rejects_non_monotonic_curves():
    """Monotonluq pozulanda diskretizasiya qeyri-stabil olur."""
    rising_oil = SaturationTable([0.2, 0.5, 0.8], [0.0, 0.1, 0.3],
                                 [0.0, 0.4, 0.8])
    falling_water = SaturationTable([0.2, 0.5, 0.8], [0.3, 0.1, 0.0],
                                    [0.8, 0.4, 0.0])
    assert any("kro" in message for message in rising_oil.validate())
    assert any("krw" in message for message in falling_water.validate())


def test_validation_rejects_out_of_range_values():
    assert SaturationTable([0.2, 0.8], [0.0, 1.4], [0.8, 0.0]).validate()
    assert SaturationTable([0.2, 0.8], [0.0, -0.1], [0.8, 0.0]).validate()
    assert SaturationTable([0.8, 0.2], [0.0, 0.3], [0.8, 0.0]).validate()


def test_slope_is_the_exact_piecewise_gradient():
    """Hamar törəmə cədvəlin özü ilə uyğun gəlmir (bax A6_PLAN.md)."""
    table = SaturationTable([0.2, 0.4, 0.8], [0.0, 0.1, 0.5],
                            [0.9, 0.5, 0.0])
    assert abs(table.slope(table.krw, 0.3)[0] - 0.5) < 1e-9     # 0.1/0.2
    assert abs(table.slope(table.krw, 0.6)[0] - 1.0) < 1e-9     # 0.4/0.4
    assert abs(table.slope(table.krw, 0.9)[0]) < 1e-12          # kənarda


# ── provider ──────────────────────────────────────────────────────────
def test_table_provider_matches_the_corey_adapter():
    """Eyni əyri iki yolla — nəticə eyni olmalıdır."""
    scal = default_scal()
    table_provider = TableRelativePermeabilityProvider(
        SaturationTableSet(tables={1: SaturationTable.from_corey(scal, 81)}))
    adapter = CoreyRelativePermeabilityAdapter(scal)
    sw = np.linspace(scal.swc + 0.02, 1.0 - scal.sor - 0.02, 20)

    assert np.max(np.abs(table_provider.krw(sw) - adapter.krw(sw))) < 1e-3
    assert np.max(np.abs(table_provider.kro(sw) - adapter.kro(sw))) < 1e-3
    assert table_provider.saturation_limits() == adapter.saturation_limits()


def test_regions_give_different_curves_at_the_same_saturation():
    tables = _two_regions()
    ids = np.array([1] * 5 + [2] * 5)
    provider = TableRelativePermeabilityProvider(tables, ids)
    krw = provider.krw(np.full(10, 0.55))

    assert np.allclose(krw[:5], krw[0])
    assert np.allclose(krw[5:], krw[5])
    assert krw[0] > krw[5] * 2, krw


def test_explicit_region_array_overrides_the_stored_region_ids():
    """TAPILAN SƏHV (audit): `krw(sw, region=<massiv>)` çağırılanda
    `region` massivi əvvəllər SƏSSİZCƏ atılırdı — yalnız konstruktorda
    verilmiş `self.region_ids`-in İLK elementi bütün sorğu üçün işlədilirdi.
    Mühərrikin özü `region` arqumentini heç vaxt ötürmür (bax
    `IRelativePermeabilityProvider` çağırış yerləri), ona görə bu, indiyə
    qədər TUTULMAMIŞDI — amma interfeys `region: Optional[np.ndarray]`
    elan edir və çağıran açıq per-hüceyrə massivi ötürə bilər."""
    tables = _two_regions()
    provider = TableRelativePermeabilityProvider(tables, region_ids=np.array([1, 1, 2]))
    sw = np.full(3, 0.55)

    default = provider.krw(sw)                              # self.region_ids: [1,1,2]
    override = provider.krw(sw, region=np.array([2, 2, 2]))  # açıq per-hüceyrə override
    mixed = provider.krw(sw, region=np.array([1, 2, 1]))

    table1, table2 = tables.get(1), tables.get(2)
    expected1 = table1.interpolate_krw(sw)
    expected2 = table2.interpolate_krw(sw)

    assert np.allclose(default, [expected1[0], expected1[0], expected2[0]])
    assert np.allclose(override, expected2)                  # HAMISI region-2
    assert np.allclose(mixed, [expected1[0], expected2[0], expected1[0]])


def test_saturation_limits_take_the_narrowest_interval():
    """Mühərrik bir hədd cütü ilə kəsir — ən məhdudlaşdırıcı lazımdır.

    Əks halda gilli zonada Sw cədvəldən kənara çıxardı.
    """
    provider = TableRelativePermeabilityProvider(_two_regions())
    low, high = provider.saturation_limits()
    # hər iki hədd ƏN MƏHDUDLAŞDIRICI regiondan gəlir — burada gilli zona
    assert abs(low - 0.30) < 1e-9          # ən böyük Swc
    assert abs(high - 0.65) < 1e-9         # ən kiçik 1−Sor


def test_per_region_limits_are_still_available():
    provider = TableRelativePermeabilityProvider(_two_regions())
    assert abs(provider.saturation_limits(1)[0] - 0.15) < 1e-9
    assert abs(provider.saturation_limits(2)[0] - 0.30) < 1e-9


def test_cfl_derivative_uses_the_worst_region():
    """Zaman addımı ən sərt hüceyrə ilə məhdudlaşmalıdır, orta ilə yox."""
    tables = _two_regions()
    provider = TableRelativePermeabilityProvider(tables)
    combined = provider.max_fractional_flow_derivative(0.5, 3.0)
    per_region = [provider.max_fractional_flow_derivative(0.5, 3.0, region)
                  for region in tables.regions]
    assert abs(combined - max(per_region)) < 1e-9


def test_invalid_tables_are_rejected_at_construction():
    broken = SaturationTableSet()
    broken.add(1, SaturationTable([0.2, 0.8], [0.5, 0.0], [0.8, 0.0]))
    try:
        TableRelativePermeabilityProvider(broken)
    except ValueError:
        return
    raise AssertionError("Yararsız cədvəl qəbul edildi")


def test_capillary_provider_reads_the_pc_column():
    tables = SaturationTableSet()
    tables.add(1, SaturationTable([0.2, 0.5, 0.8], [0.0, 0.1, 0.35],
                                  [0.8, 0.3, 0.0], pc=[1.0, 0.4, 0.0]))
    provider = TableCapillaryPressureProvider(tables)
    assert provider.has_capillary_pressure()
    assert abs(float(provider.pcow(np.array([0.35]))[0]) - 0.7) < 1e-9
    assert float(provider.dpcow_dsw(np.array([0.35]))[0]) < 0


def test_capillary_provider_is_zero_without_a_pc_column():
    provider = TableCapillaryPressureProvider(_set())
    assert not provider.has_capillary_pressure()
    assert np.allclose(provider.pcow(np.array([0.3, 0.5])), 0.0)


# ── mühərriklə inteqrasiya ────────────────────────────────────────────
def test_engine_runs_with_table_provider():
    from imex2d.application.simulation_service import SimulationService

    scal = default_scal()
    model = five_spot_model(nx=9, ny=9, scal=scal)
    provider = TableRelativePermeabilityProvider(
        SaturationTableSet(tables={1: SaturationTable.from_corey(scal, 81)}))
    result = SimulationService(relperm_provider=provider).run(
        model, short_config(end_time=200.0))
    assert result.converged
    assert result.final_recovery_factor > 0


def test_table_and_corey_give_the_same_recovery():
    """Eyni əyri — eyni nəticə. Cədvəl yolu fizikanı dəyişməməlidir."""
    from imex2d.application.simulation_service import SimulationService

    scal = default_scal()
    table_result = SimulationService(
        relperm_provider=TableRelativePermeabilityProvider(
            SaturationTableSet(
                tables={1: SaturationTable.from_corey(scal, 201)}))).run(
        five_spot_model(nx=9, ny=9, scal=scal), short_config(end_time=300.0))
    corey_result = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal)).run(
        five_spot_model(nx=9, ny=9, scal=scal), short_config(end_time=300.0))

    assert abs(table_result.final_recovery_factor
               - corey_result.final_recovery_factor) < 0.5


# ── CSV / SWOF ────────────────────────────────────────────────────────
def test_reads_csv_with_regions():
    path = _write("""region,sw,krw,kro,pc
1,0.20,0.000,0.800,0.5
1,0.50,0.100,0.300,0.1
1,0.80,0.350,0.000,0.0
2,0.30,0.000,0.750,0.8
2,0.65,0.200,0.000,0.0
""")
    try:
        tables = read_scal_csv(path)
        assert tables.regions == [1, 2]
        assert tables.get(1).has_capillary
        assert abs(tables.get(2).swc - 0.30) < 1e-9
    finally:
        os.unlink(path)


def test_csv_without_region_column_uses_one_region():
    path = _write("sw,krw,kro\n0.2,0.0,0.8\n0.5,0.1,0.3\n0.8,0.35,0.0\n")
    try:
        assert read_scal_csv(path).regions == [1]
    finally:
        os.unlink(path)


def test_csv_rows_are_sorted_by_saturation():
    path = _write("sw,krw,kro\n0.8,0.35,0.0\n0.2,0.0,0.8\n0.5,0.1,0.3\n")
    try:
        table = read_scal_csv(path).get(1)
        assert list(table.sw) == [0.2, 0.5, 0.8]
    finally:
        os.unlink(path)


def test_csv_missing_columns_is_rejected():
    path = _write("sw,krw\n0.2,0.0\n0.8,0.35\n")
    try:
        read_scal_csv(path)
    except ScalFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("Natamam CSV qəbul edildi")


def test_non_monotonic_file_is_rejected_on_read():
    path = _write("sw,krw,kro\n0.2,0.3,0.8\n0.5,0.1,0.3\n0.8,0.0,0.0\n")
    try:
        read_scal_csv(path)
    except ScalFormatError as error:
        assert "krw" in str(error)
        return
    finally:
        os.unlink(path)
    raise AssertionError("Monoton olmayan cədvəl qəbul edildi")


def test_reads_eclipse_swof_with_two_regions():
    path = _write("""PROPS

SWOF
-- Sw    krw    kro    Pc
  0.20  0.000  0.800  0.5
  0.50  0.100  0.300  0.1
  0.80  0.350  0.000  0.0
/
  0.30  0.000  0.750  0.8
  0.55  0.060  0.250  0.2
  0.65  0.200  0.000  0.0
/

DENSITY
  850 1010 1.0 /
""", suffix=".DATA")
    try:
        tables = read_swof(path)
        assert tables.regions == [1, 2]
        assert tables.get(1).sw.size == 3
        assert tables.get(2).has_capillary
    finally:
        os.unlink(path)


def test_swof_without_the_keyword_is_rejected():
    path = _write("PROPS\n\nDENSITY\n 850 1010 1.0 /\n", suffix=".DATA")
    try:
        read_swof(path)
    except ScalFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("SWOF olmadan fayl qəbul edildi")


def test_csv_round_trip_preserves_the_curves():
    tables = _two_regions()
    handle, path = tempfile.mkstemp(suffix=".csv")
    os.close(handle)
    try:
        write_scal_csv(path, tables)
        restored = read_scal_csv(path)
    finally:
        os.unlink(path)

    assert restored.regions == tables.regions
    for region in tables.regions:
        original, copy = tables.get(region), restored.get(region)
        assert np.allclose(copy.sw, original.sw, atol=1e-5)
        assert np.allclose(copy.krw, original.krw, atol=1e-5)
        assert np.allclose(copy.kro, original.kro, atol=1e-5)


def test_written_deck_can_be_read_back_as_scal():
    """Eclipse ixracı ilə SCAL importu bir-birini tamamlamalıdır."""
    from imex2d.io.eclipse_export import EclipseDeckWriter

    model = five_spot_model(nx=7, ny=7, scal=default_scal())
    handle, path = tempfile.mkstemp(suffix=".DATA")
    os.close(handle)
    try:
        EclipseDeckWriter().write(model, path)
        tables = read_swof(path)
    finally:
        os.unlink(path)

    assert tables.regions == [1]
    table = tables.get(1)
    assert abs(table.swc - model.scal_parameters.swc) < 0.02
    assert abs(table.krw_end - model.scal_parameters.krw_end) < 1e-3
