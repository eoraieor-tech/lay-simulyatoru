"""Phase 4.1 — xassə növü reyestri (CONTINUOUS vs CATEGORICAL)."""

from __future__ import annotations

from imex2d.geology.property_types import (PropertyType, classify_property,
                                            is_categorical)


def test_known_continuous_properties():
    for name in ("PORO", "PERMX", "PERMY", "PERMZ", "NTG", "SW", "VSH", "PRESSURE"):
        assert classify_property(name) is PropertyType.CONTINUOUS
        assert not is_categorical(name)


def test_known_categorical_properties():
    for name in ("FACIES", "LITHOLOGY", "ROCKTYPE"):
        assert classify_property(name) is PropertyType.CATEGORICAL
        assert is_categorical(name)


def test_case_insensitive():
    assert is_categorical("facies")
    assert is_categorical("Facies")


def test_unknown_property_defaults_to_continuous_not_silently_categorical():
    assert classify_property("SOME_NEW_PROPERTY") is PropertyType.CONTINUOUS


def test_overrides_extend_registry_without_mutating_default():
    overrides = {"CUSTOM_FACIES": PropertyType.CATEGORICAL}
    assert is_categorical("CUSTOM_FACIES", overrides)
    assert not is_categorical("CUSTOM_FACIES")   # defolt reyestr TOXUNULMAYIB


def test_overrides_can_reclassify_a_default_entry():
    overrides = {"NTG": PropertyType.CATEGORICAL}
    assert classify_property("NTG", overrides) is PropertyType.CATEGORICAL
    assert classify_property("NTG") is PropertyType.CONTINUOUS
