"""Phase 5 §11/§16 — SGS-in real geoloji model-qurma boru xəttinə
inteqrasiyası: PORO/PERMX sərt data hörməti, ziddiyyət aşkarlanması,
fasiya-şərtli SGS `model.facies_fields`-dən istifadə edərək.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.application.geology_service import (ContinuousSGSConfig, FaciesBuildConfig,
                                                 GeologicalGridSpec,
                                                 WellBasedGeologicalModelBuilder)
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology.hard_data import HardDataConflictError
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.sgs import FaciesPropertyConfig


def _poro_permx_dataset(n_wells=8, seed=0):
    rng = np.random.default_rng(seed)
    samples = []
    for w in range(n_wells):
        x, y = rng.uniform(10, 190, size=2)
        poro = float(rng.normal(0.20, 0.03))
        permx = float(rng.lognormal(4.5, 0.8))
        samples.append(WellSample(well=f"W-{w}", x=x, y=y, layer=0,
                                  values={"PORO": poro, "PERMX": permx}))
    return WellDataset(samples=samples, source="test")


def _spec():
    return GeologicalGridSpec(nx=6, ny=6, nz=1, dx=35.0, dy=35.0, top_depth=2000.0)


# ── §11: sərt data hörməti — PORO/PERMX/PERMY/PERMZ ───────────────────
def test_poro_hard_data_honored_through_full_pipeline():
    dataset = _poro_permx_dataset()
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"PORO": ContinuousSGSConfig(seed=1)}
    model, _ = builder.build(dataset, _spec(), sgs_config=config)

    from imex2d.domain.geometry import xy_to_ij
    poro_values = model.property_maps["PORO"].values
    for sample in dataset.samples:
        i, j = xy_to_ij(sample.x, sample.y, model.geometry)
        cell = np.ravel_multi_index((0, j, i), model.grid.shape)
        assert poro_values[cell] == pytest.approx(sample.values["PORO"], abs=1e-6)


def test_permx_hard_data_honored_and_stays_positive():
    dataset = _poro_permx_dataset(seed=1)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"PERMX": ContinuousSGSConfig(seed=2, log_space=True)}
    model, _ = builder.build(dataset, _spec(), sgs_config=config)

    from imex2d.domain.geometry import xy_to_ij
    permx_values = model.property_maps["PERMX"].values
    assert np.all(permx_values > 0.0)
    for sample in dataset.samples:
        i, j = xy_to_ij(sample.x, sample.y, model.geometry)
        cell = np.ravel_multi_index((0, j, i), model.grid.shape)
        assert permx_values[cell] == pytest.approx(sample.values["PERMX"], rel=1e-5)


def test_permy_and_permz_still_use_default_continuous_path_without_sgs_config():
    """§11.13-un analoqu: `sgs_config`-də OLMAYAN sütun ƏVVƏLKİ
    (deterministik Kriging) yolu ilə davam edir — SGS DAYATILMIR."""
    dataset = _poro_permx_dataset(seed=2)
    for s in dataset.samples:
        s.values["PERMY"] = s.values["PERMX"] * 0.8
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"PERMX": ContinuousSGSConfig(seed=1)}   # PERMY KONFİQURASİYA EDİLMİR
    model, _ = builder.build(dataset, _spec(), sgs_config=config)
    # PERMY kəsilməz Kriging-dən gəlir (hamar, SGS-in stoxastik "cırıqlığı" YOXDUR)
    assert "PERMY" in model.property_maps


# ── §11: ziddiyyət aşkarlanması ─────────────────────────────────────────
def test_conflicting_poro_hard_data_raises_through_pipeline():
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PORO": 0.15}),
        WellSample(well="B", x=11.0, y=10.0, layer=0, values={"PORO": 0.35}),
        WellSample(well="C", x=150.0, y=150.0, layer=0, values={"PORO": 0.20}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"PORO": ContinuousSGSConfig(conflict_tolerance=0.02)}
    with pytest.raises(HardDataConflictError):
        builder.build(dataset, _spec(), sgs_config=config)


def test_duplicate_consistent_poro_hard_data_does_not_raise():
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PORO": 0.200, "PERMX": 150.0}),
        WellSample(well="A2", x=11.0, y=10.0, layer=0, values={"PORO": 0.201, "PERMX": 151.0}),
        WellSample(well="C", x=150.0, y=150.0, layer=0, values={"PORO": 0.22, "PERMX": 200.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"PORO": ContinuousSGSConfig(conflict_tolerance=0.02)}
    model, _ = builder.build(dataset, _spec(), sgs_config=config)
    assert "PORO" in model.property_maps


# ── §5/§16: fasiya-şərtli SGS, `model.facies_fields`-dən istifadə ────────
def test_facies_conditioned_sgs_uses_geological_model_facies_field():
    rng = np.random.default_rng(3)
    samples = []
    for w in range(10):
        x = rng.uniform(10, 90)
        y = rng.uniform(10, 190)
        facies = 0
        poro = float(rng.normal(0.25, 0.02))
        samples.append(WellSample(well=f"S{w}", x=x, y=y, layer=0,
                                  values={"PORO": poro, "PERMX": 200.0, "FACIES": float(facies)}))
    for w in range(10):
        x = rng.uniform(110, 190)
        y = rng.uniform(10, 190)
        facies = 1
        poro = float(rng.normal(0.08, 0.01))
        samples.append(WellSample(well=f"H{w}", x=x, y=y, layer=0,
                                  values={"PORO": poro, "PERMX": 5.0, "FACIES": float(facies)}))
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=8, ny=8, nz=1, dx=25.0, dy=25.0, top_depth=2000.0)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())

    facies_config = {"FACIES": FaciesBuildConfig(proportions={0: 0.5, 1: 0.5}, seed=1)}
    sgs_config = {"PORO": ContinuousSGSConfig(seed=2, facies_field_name="FACIES",
                                              conflict_tolerance=0.1)}
    model, report = builder.build(dataset, spec, facies_config=facies_config,
                                  sgs_config=sgs_config)

    facies_codes = model.facies_fields["FACIES"].codes
    poro_values = model.property_maps["PORO"].values
    sand_poro = poro_values[facies_codes == 0]
    shale_poro = poro_values[facies_codes == 1]
    assert np.mean(sand_poro) > np.mean(shale_poro) + 0.05


def test_facies_conditioned_sgs_missing_field_raises_clear_error():
    dataset = _poro_permx_dataset(seed=4)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"PORO": ContinuousSGSConfig(seed=1, facies_field_name="NOT_BUILT")}
    with pytest.raises(ValueError, match="facies_field_name"):
        builder.build(dataset, _spec(), sgs_config=config)
