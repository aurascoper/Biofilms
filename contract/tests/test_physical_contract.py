import pytest

from physical_contract import (ALLOWED_BOUNDARY_CONDITIONS,
                               ALLOWED_SOURCE_ANGULAR, ALLOWED_SOURCE_SPATIAL,
                               EVIDENCE_POLICIES, EXECUTION_CLASSES,
                               MaterialSpec, closed_composition_problems,
                               git_provenance, render_material_toml)

WATER = {"H": 0.111894, "O": 0.888106}


def test_loader_and_contract_share_one_vocabulary():
    # The whole reason this package exists: if the coupling loader kept private
    # copies, a boundary type added on one side would silently not exist on the
    # other. Importing the loader here is fine — coupling may depend on the
    # contract; it is calibration that may not depend on coupling.
    from biofilm_openmc import config

    assert config._ALLOWED_BC is ALLOWED_BOUNDARY_CONDITIONS
    assert config._ALLOWED_SPATIAL is ALLOWED_SOURCE_SPATIAL
    assert config._ALLOWED_ANGULAR is ALLOWED_SOURCE_ANGULAR


def test_evidence_policy_is_not_a_system_provenance_value():
    # `system_provenance` answers "what kind of system"; `evidence_policy`
    # answers "may this run use unmeasured values". An engineered composite can
    # hold measured, certified, derived and declared components at once, so
    # collapsing the two axes would make that system inexpressible.
    provenance = {"published_replica", "certified_component",
                  "engineered_composite", "declared"}
    assert not (EVIDENCE_POLICIES & provenance)
    assert "synthetic_validation" in EXECUTION_CLASSES


def test_open_composition_is_refused():
    assert not closed_composition_problems(WATER)
    problems = closed_composition_problems({"H": 0.2, "O": 0.5})
    assert len(problems) == 1 and "not 1" in problems[0]
    assert "negative" in closed_composition_problems({"H": -0.1, "O": 1.1})[0]
    assert closed_composition_problems({}) == ["no elemental composition"]


def test_render_round_trips_through_tomllib():
    import tomllib

    spec = MaterialSpec(name="liquid water", density_g_cm3=1.0, elements=WATER)
    data = tomllib.loads(render_material_toml(spec, "medium", header="declared"))
    assert data["materials"]["medium"]["density_g_cm3"] == 1.0
    assert data["materials"]["medium"]["elements"] == WATER


def test_render_refuses_an_unusable_material():
    # Rendering is where a composition stops being a dict and starts being a
    # file the transport loader will read, so it is the last place to refuse.
    with pytest.raises(ValueError, match="not 1"):
        render_material_toml(MaterialSpec("x", 1.0, {"H": 0.5}), "medium")
    with pytest.raises(ValueError, match="not positive"):
        render_material_toml(MaterialSpec("x", 0.0, WATER), "medium")


def test_git_provenance_marks_a_dirty_tree_in_the_commit_string():
    got = git_provenance()
    assert set(got) == {"git_commit", "git_dirty"}
    if got["git_commit"] is not None:
        # The marker rides inside the commit string precisely so it cannot be
        # dropped by a consumer that only reads one field.
        assert got["git_commit"].endswith("-dirty") == got["git_dirty"]


def test_git_provenance_on_a_non_repository_is_null_not_a_guess(tmp_path):
    assert git_provenance(str(tmp_path)) == {"git_commit": None, "git_dirty": None}
