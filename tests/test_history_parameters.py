"""Uyğunlaşdırma parametrləri və model modifikatoru (C5/2)."""

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.history.mismatch import MismatchCalculator
from imex2d.history.parameters import (ModelModifier, ParameterDefinition,
                                       ParameterKind, ParameterSet,
                                       standard_parameters)
from test_history_mismatch import _observations_from
from test_implicit_newton import _rate_controlled


def _model(nx=9):
    return _rate_controlled(nx=nx, scal=default_scal())


def _set(model=None):
    model = model or _model()
    return model, ParameterSet(standard_parameters(model))


# ── tərif ─────────────────────────────────────────────────────────────
def test_definition_rejects_inverted_bounds():
    try:
        ParameterDefinition("X", lambda m, v: None, minimum=2.0, maximum=1.0)
    except ValueError:
        return
    raise AssertionError("Tərs hədlər qəbul edildi")


def test_log_scale_rejects_non_positive_minimum():
    try:
        ParameterDefinition("X", lambda m, v: None, minimum=0.0, maximum=10.0,
                            log_scale=True)
    except ValueError:
        return
    raise AssertionError("Log miqyasda sıfır minimum qəbul edildi")


def test_initial_value_is_clipped_into_the_bounds():
    definition = ParameterDefinition("X", lambda m, v: None,
                                     minimum=0.5, maximum=2.0, initial=9.0)
    assert definition.initial == 2.0


def test_unit_mapping_round_trips():
    for log_scale, minimum, maximum in ((False, 0.5, 2.0), (True, 0.1, 10.0)):
        definition = ParameterDefinition("X", lambda m, v: None,
                                         minimum=minimum, maximum=maximum,
                                         log_scale=log_scale)
        for value in np.linspace(minimum, maximum, 7):
            assert abs(definition.from_unit(definition.to_unit(value))
                       - value) < 1e-9


def test_unit_mapping_is_bounded():
    definition = ParameterDefinition("X", lambda m, v: None,
                                     minimum=1.0, maximum=3.0)
    assert definition.to_unit(-5.0) == 0.0
    assert definition.to_unit(99.0) == 1.0
    assert definition.from_unit(-1.0) == 1.0
    assert definition.from_unit(2.0) == 3.0


def test_log_scale_makes_reciprocal_factors_symmetric():
    """0.5 və 2.0 fiziki cəhətdən simmetrik dəyişikliklərdir.

    Xətti miqyasda 1.0-dan məsafələri fərqlidir (0.5 vs 1.0), ona görə
    optimallaşdırıcı azaltmağa üstünlük verərdi. Log fəzada hər ikisi
    mərkəzdən eyni məsafədədir.
    """
    definition = ParameterDefinition("X", lambda m, v: None,
                                     minimum=0.1, maximum=10.0, log_scale=True)
    low = definition.to_unit(0.5)
    high = definition.to_unit(2.0)
    assert abs((0.5 - low) - (high - 0.5)) < 1e-9


# ── standart dəst ─────────────────────────────────────────────────────
def test_standard_set_covers_the_main_uncertainties():
    _, parameters = _set()
    names = parameters.names
    for expected in ("PERM_MULT", "PORO_MULT", "SOR", "SWC", "KRW_END"):
        assert expected in names, expected


def test_vertical_ratio_appears_only_with_permz():
    model = _model()
    assert "KV_KH" in ParameterSet(standard_parameters(model)).names
    model.rock.permz = None
    assert "KV_KH" not in ParameterSet(standard_parameters(model)).names


def test_viscosity_parameter_disappears_when_pvt_is_present():
    from imex2d.simulation.pvt.correlations import build_pvt_table

    model = _model()
    assert "MU_OIL" in ParameterSet(standard_parameters(model)).names
    model.pvt_table = build_pvt_table(bubble_point_bar=150.0)
    assert "MU_OIL" not in ParameterSet(standard_parameters(model)).names


def test_contact_parameter_is_bounded_by_the_model_depths():
    import dataclasses

    model = _model()
    surface = np.add.outer(np.zeros(model.grid.ny),
                           np.arange(model.grid.nx) * 5.0) + 2000.0
    model.geometry = dataclasses.replace(model.geometry,
                                         top_depth_map=surface)
    model.initial_conditions.oil_water_contact = 2020.0
    definition = next(item for item in standard_parameters(model)
                      if item.name == "OWC")
    depths = model.geometry.cell_depths()
    assert abs(definition.minimum - depths.min()) < 1e-9
    assert abs(definition.maximum - depths.max()) < 1e-9


def test_flat_model_still_gets_a_usable_contact_range():
    """Düz layda bütün dərinliklər eynidir — hədlər genişləndirilir."""
    model = _model()
    model.initial_conditions.oil_water_contact = 5.0
    definition = next(item for item in standard_parameters(model)
                      if item.name == "OWC")
    assert definition.maximum > definition.minimum


def test_initial_values_reproduce_the_base_model():
    model, parameters = _set()
    for definition, value in zip(parameters.definitions,
                                 parameters.initial_values):
        if definition.kind is ParameterKind.MULTIPLIER:
            assert abs(value - 1.0) < 1e-9, definition.name


# ── tətbiq ────────────────────────────────────────────────────────────
def test_multiplier_scales_the_field_and_keeps_heterogeneity():
    model = _model()
    model.rock.permx.values[:] = np.linspace(50.0, 500.0, model.ncell)
    model.rock.permy.values[:] = model.rock.permx.values.copy()
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)

    original = model.rock.permx.values.copy()
    values = parameters.initial_values.copy()
    values[parameters.names.index("PERM_MULT")] = 3.0
    modified = modifier.apply(values)

    assert np.allclose(modified.rock.permx.values, original * 3.0)
    # heterogenliyin forması qorunur
    assert np.allclose(modified.rock.permx.values / original.mean() / 3.0,
                       original / original.mean())


def test_base_model_is_never_modified():
    """Ən vacib xassə: optimallaşdırma yüzlərlə dəfə tətbiq edir."""
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)

    original_perm = model.rock.permx.values.copy()
    original_sor = model.scal_parameters.sor

    values = parameters.initial_values.copy()
    values[parameters.names.index("PERM_MULT")] = 4.0
    values[parameters.names.index("SOR")] = 0.4
    for _ in range(5):
        modifier.apply(values)

    assert np.allclose(model.rock.permx.values, original_perm)
    assert model.scal_parameters.sor == original_sor


def test_repeated_application_is_deterministic():
    """Çarpanlar üst-üstə yığılmamalıdır: 3.0 → 3.0, 9.0 deyil."""
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)
    values = parameters.initial_values.copy()
    values[parameters.names.index("PERM_MULT")] = 3.0

    first = modifier.apply(values).rock.permx.values.copy()
    second = modifier.apply(values).rock.permx.values
    assert np.allclose(first, second)


def test_modified_model_does_not_share_arrays_with_the_base():
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    modified = ModelModifier(model, parameters).apply(
        parameters.initial_values)
    modified.rock.permx.values[0] = 99999.0
    assert model.rock.permx.values[0] != 99999.0


def test_absolute_parameters_are_set_not_scaled():
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    values = parameters.initial_values.copy()
    values[parameters.names.index("SOR")] = 0.33
    modified = ModelModifier(model, parameters).apply(values)
    assert abs(modified.scal_parameters.sor - 0.33) < 1e-9


def test_vertical_ratio_rebuilds_permz_from_horizontal():
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    values = parameters.initial_values.copy()
    values[parameters.names.index("KV_KH")] = 0.02
    modified = ModelModifier(model, parameters).apply(values)

    horizontal = np.sqrt(modified.rock.permx.values
                         * modified.rock.permy.values)
    assert np.allclose(modified.rock.permz.values, horizontal * 0.02)


def test_values_outside_the_bounds_are_clipped():
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    values = parameters.initial_values.copy()
    values[parameters.names.index("SOR")] = 5.0
    modified = ModelModifier(model, parameters).apply(values)
    assert modified.scal_parameters.sor <= 0.45 + 1e-12


def test_wrong_value_count_is_rejected():
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    try:
        ModelModifier(model, parameters).apply(np.ones(2))
    except ValueError:
        return
    raise AssertionError("Yanlış sayda dəyər qəbul edildi")


def test_modified_model_still_passes_validation_and_runs():
    scal = default_scal()
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    values = parameters.initial_values.copy()
    values[parameters.names.index("PERM_MULT")] = 2.0
    values[parameters.names.index("SOR")] = 0.3

    modified = ModelModifier(model, parameters).apply(values)
    assert modified.validate() == []
    result = make_service(scal).run(modified, short_config(end_time=100.0))
    assert result.converged


# ── uyğunsuzluğa təsir ────────────────────────────────────────────────
def test_parameters_measurably_change_the_mismatch():
    """Parametrlər nəticəyə təsir etmirsə, optimallaşdırma mənasızdır.

    Müddət su gəlişindən sonranı əhatə etməlidir: SCAL parametrləri
    (Sor, krw) yalnız iki faza birlikdə hərəkət edəndə özünü göstərir.
    """
    scal = default_scal()
    end_time = 1500.0
    model = _model(nx=11)
    truth = make_service(scal).run(model, short_config(end_time=end_time,
                                                       snapshots=5))
    assert truth.breakthrough_time is not None, "Ssenaridə su gəlişi yoxdur"

    observations = _observations_from(truth)
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)
    calculator = MismatchCalculator()
    base_score = calculator.evaluate(truth, observations).total

    for name, value in (("PERM_MULT", 0.25), ("PORO_MULT", 1.3),
                        ("SOR", 0.40), ("KRW_END", 0.9)):
        values = parameters.initial_values.copy()
        values[parameters.names.index(name)] = value
        result = make_service(scal).run(modifier.apply(values),
                                        short_config(end_time=end_time))
        score = calculator.evaluate(result, observations).total
        assert score > base_score + 0.01, (name, base_score, score)


def test_scal_changes_reach_the_simulation_engine():
    """Provider modeldən qurulmalıdır, konstruktorda saxlanılmamalı.

    Əks halda optimallaşdırıcı Sor-u dəyişir, model yenilənir, lakin
    nisbi keçiricilik adapteri köhnə Sor ilə qalır — nəticə heç
    dəyişmir və axtarış mənasız olur.
    """
    scal = default_scal()
    model = _model(nx=11)
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)

    base = make_service(scal).run(model, short_config(end_time=1500.0))
    values = parameters.initial_values.copy()
    values[parameters.names.index("SOR")] = 0.40
    changed = make_service(scal).run(modifier.apply(values),
                                     short_config(end_time=1500.0))

    assert abs(changed.final_recovery_factor
               - base.final_recovery_factor) > 1.0
    assert changed.breakthrough_time < base.breakthrough_time


def test_changing_swc_keeps_the_model_valid():
    """Swc dəyişəndə ilkin Sw avtomatik uyğunlaşdırılmalıdır.

    Əks halda model yoxlamadan keçmir və optimallaşdırma dayanır.
    """
    scal = default_scal()
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)

    for swc in (0.06, 0.35):
        values = parameters.initial_values.copy()
        values[parameters.names.index("SWC")] = swc
        modified = modifier.apply(values)
        assert modified.validate() == [], (swc, modified.validate())
        assert make_service(scal).run(
            modified, short_config(end_time=100.0)).converged


def test_unit_space_covers_the_whole_range():
    model = _model()
    parameters = ParameterSet(standard_parameters(model))
    low = parameters.from_unit(np.zeros(len(parameters)))
    high = parameters.from_unit(np.ones(len(parameters)))
    for definition, minimum, maximum in zip(parameters.definitions, low, high):
        assert abs(minimum - definition.minimum) < 1e-9
        assert abs(maximum - definition.maximum) < 1e-9
