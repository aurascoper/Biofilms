"""The apparatus half of a transport config, and what `synthetic` may skip.

The load-bearing assertion in this file is `test_synthetic_skips_only_evidence`:
if `evidence_policy = "synthetic"` ever relaxes a structural or coherence check,
it has become a universal escape hatch around the contracts that exist to stop
incoherent configurations from being built.
"""

from __future__ import annotations

import dataclasses
import tomllib

import pytest

from biofilm_calibration.integration import VoxelBinding
from biofilm_calibration.materials.export import MaterialClass
from biofilm_calibration.reference_system import (ReferenceSystemSpec,
                                                  emit_problems,
                                                  emit_transport_config,
                                                  evidence_problems,
                                                  structural_problems)

BINDING = VoxelBinding(
    entity_kind="biomass_parcel",
    lineage_semantics="computational",
    material_class="baseline_biomass",
    material_represents="hydrated_bulk_biofilm",
    material_model_kind="hydrated_effective_medium",
    density_basis="wet_bulk",
    lattice_pitch_um=12000.0,          # 1.2 cm — dosimetry scale, NOT a biofilm pitch
    density_g_cm3=1.05,
)

WATER = {"H": 0.111894, "O": 0.888106}


def _material(name, density, elements, evidence="declared"):
    return MaterialClass(name=name, density_g_cm3=density, elements=elements,
                         evidence_basis=evidence, source_id="SYNTHETIC_DECL",
                         system_conditions="invented for a software path test")


SPEC = ReferenceSystemSpec(
    reference_system_id="synthetic_biofilm_e2e",
    system_provenance="declared",
    evidence_policy="synthetic",
    execution_class="synthetic_validation",
    target_calibration=False,
    lattice_n=20,
    origin_cm=(0.0, 0.0, 0.0),
    cylinder_radius_cm=11.0,
    cylinder_length_cm=24.0,
    membrane_thickness_cm=0.5,
    axial_bc="vacuum",
    radial_outer_bc="vacuum",
    source_spatial="point_origin",
    source_angular="isotropic",
    source_position_cm=(12.3, 12.3, 12.3),
    spectrum_energies_eV=(661655.0,),
    spectrum_probabilities=(1.0,),
    medium=_material("water", 1.0, WATER),
    membrane=_material("invented membrane", 0.94, {"C": 0.856, "H": 0.144}),
    biomass_elements=WATER,
    batches=80,
    particles=50000,
    seed=1,
    mesh_coarsening_factor=2,
    notes="synthetic reference system; calibrates nothing",
)

GATES = {"spatial_verdict": "PROVISIONAL",
         "material_openmc_verdict": "PROVISIONAL"}


def test_the_synthetic_spec_is_structurally_sound():
    assert structural_problems(BINDING, SPEC) == []


# --- invariant 1: no silent defaults -------------------------------------

_NOT_REQUIRED = {"photons_per_second", "notes"}


@pytest.mark.parametrize("field", [f.name for f in dataclasses.fields(SPEC)
                                   if f.name not in _NOT_REQUIRED])
def test_every_required_field_is_named_when_it_is_blank(field):
    """Deleting any required field must produce a refusal that NAMES it. A
    refusal that says 'something is missing' sends you hunting; a refusal that
    says which field is a fix."""
    blanked = dataclasses.replace(SPEC, **{field: None})
    problems = " ".join(structural_problems(BINDING, blanked))
    assert field in problems


def test_dosimetry_requires_a_source_rate_and_transport_does_not():
    assert structural_problems(BINDING, SPEC, stage="transport") == []
    problems = structural_problems(BINDING, SPEC, stage="dosimetry")
    assert any("photons_per_second is unset" in p for p in problems)

    with_rate = dataclasses.replace(SPEC, photons_per_second=1.0e12)
    assert structural_problems(BINDING, with_rate, stage="dosimetry") == []


# --- amendment A3: what synthetic may and may not skip -------------------

def test_synthetic_skips_only_evidence():
    """`synthetic` bypasses the evidence gate and NOTHING else."""
    # The gates are PROVISIONAL and the materials are `declared`, so a
    # measured_only policy is blocked on evidence...
    measured = dataclasses.replace(SPEC, evidence_policy="measured_only")
    blocked = emit_problems(BINDING, measured, **GATES)
    assert any("spatial gate is PROVISIONAL" in p for p in blocked)
    assert any("cannot be exported as a calibrated material" in p for p in blocked)

    # ...while the synthetic policy passes the identical inputs.
    assert emit_problems(BINDING, SPEC, **GATES) == []

    # But a geometry that cannot be built is refused under BOTH policies.
    for policy in ("synthetic", "measured_only"):
        impossible = dataclasses.replace(SPEC, evidence_policy=policy,
                                         cylinder_radius_cm=40.0)
        assert any("exceeds the lattice half-width" in p
                   for p in emit_problems(BINDING, impossible, **GATES))


def test_synthetic_may_not_claim_to_calibrate_the_target():
    claiming = dataclasses.replace(SPEC, target_calibration=True)
    assert any("cannot calibrate the target system" in p
               for p in structural_problems(BINDING, claiming))


def test_a_source_on_a_lattice_plane_is_refused_even_when_synthetic():
    # 12.0 cm = 10 lattice planes at a 1.2 cm pitch: the derived centre.
    on_plane = dataclasses.replace(SPEC, source_position_cm=(12.0, 12.0, 12.0))
    problems = structural_problems(BINDING, on_plane)
    assert sum("lies on a lattice plane" in p for p in problems) == 3


def test_an_open_composition_is_refused_even_when_synthetic():
    open_comp = dataclasses.replace(SPEC, biomass_elements={"H": 0.2, "O": 0.5})
    assert any("baseline_biomass" in p and "not 1" in p
               for p in structural_problems(BINDING, open_comp))


def test_membrane_reaching_outside_the_tallied_cube_is_refused():
    # radius 11 + thickness 1.5 = 12.5 > half-width 12.0
    fat = dataclasses.replace(SPEC, membrane_thickness_cm=1.5)
    assert any("reaches outside the tallied lattice cube" in p
               for p in structural_problems(BINDING, fat))


def test_a_mesh_factor_that_does_not_divide_the_lattice_is_refused():
    bad = dataclasses.replace(SPEC, mesh_coarsening_factor=3)   # 20 % 3 != 0
    assert any("does not divide the lattice size" in p
               for p in structural_problems(BINDING, bad))


def test_an_unnormalized_spectrum_is_refused():
    # The loader does not normalize either — the total emission rate belongs in
    # photons_per_second, not in the PMF.
    raw_yields = dataclasses.replace(SPEC, spectrum_probabilities=(0.8501,))
    assert any("not 1" in p for p in structural_problems(BINDING, raw_yields))


# --- invariant 2: the emitted config is complete and loadable ------------

def test_emission_refuses_and_names_every_blocker():
    blank = VoxelBinding(**{**BINDING.__dict__, "lattice_pitch_um": None,
                            "density_g_cm3": None})
    with pytest.raises(ValueError, match="refusing to emit") as exc:
        emit_transport_config(blank, SPEC, **GATES)
    message = exc.value.args[0]
    assert "lattice_pitch_um is unset" in message
    assert "density_g_cm3 is unset" in message


def test_emitted_config_is_complete_toml_not_a_fragment():
    """The previous emitter produced a pitch and one density. That is not a
    config: the loader also requires the membrane thickness, the cylinder, the
    boundaries, the source and all three material classes."""
    data = tomllib.loads(emit_transport_config(BINDING, SPEC, **GATES))

    assert data["provenance"] == {
        "reference_system_id": "synthetic_biofilm_e2e",
        "system_provenance": "declared",
        "evidence_policy": "synthetic",
        "execution_class": "synthetic_validation",
        "target_calibration": False}
    assert data["geometry"]["voxel_pitch_cm"] == pytest.approx(1.2)
    assert data["geometry"]["membrane_thickness_cm"] == 0.5
    assert data["source"]["position_cm"] == [12.3, 12.3, 12.3]
    assert set(data["materials"]) == {"medium", "baseline_biomass", "membrane"}
    assert data["materials"]["baseline_biomass"]["density_g_cm3"] == 1.05
    assert data["transport"]["mesh"]["coarsening_factor"] == 2
    # transport reads no source activity, so none is written
    assert "photons_per_second" not in data["source"]


def test_dosimetry_emission_adds_the_source_rate():
    with_rate = dataclasses.replace(SPEC, photons_per_second=1.0e12)
    data = tomllib.loads(emit_transport_config(BINDING, with_rate,
                                               stage="dosimetry", **GATES))
    assert data["source"]["photons_per_second"] == 1.0e12


def test_emission_is_byte_identical_across_calls():
    """A config that changes between two identical emissions cannot be a
    reproducibility anchor."""
    assert (emit_transport_config(BINDING, SPEC, **GATES)
            == emit_transport_config(BINDING, SPEC, **GATES))


def test_a_synthetic_config_says_so_in_its_own_text():
    """The label has to survive the file being copied out of context."""
    text = emit_transport_config(BINDING, SPEC, **GATES)
    assert "SYNTHETIC" in text
    assert "calibrates nothing" in text
