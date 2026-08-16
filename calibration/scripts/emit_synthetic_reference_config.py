#!/usr/bin/env python
"""Emit the SYNTHETIC reference system's transport configs.

    reference_system_id = synthetic_biofilm_e2e
    system_provenance   = declared
    evidence_policy     = synthetic
    execution_class     = synthetic_validation
    target_calibration  = false

WHAT THIS IS FOR. Reference D's two central numbers — the CPM lattice pitch and
the wet bulk biofilm density — are unmeasured, so the target biofilm dose path
cannot run. But the SOFTWARE path downstream of those numbers has never run
either, on anything. This system supplies invented values so the whole chain
executes now, while a failure is free, instead of first meeting the pipeline on
the day expensive paired measurements arrive.

WHAT IT IS NOT. Every physical value below is invented. None of it enters
data/parameter_provenance.csv or data/calibration/voxel_binding.csv, and a
successful run licenses exactly one claim —

    "the one-way physical biofilm-dose architecture executed end to end under a
     synthetic, non-calibrating reference system"

— and never "the target biofilm was calibrated or dosed".

WHY THE PITCH IS 1.2 cm AND NOT A BIOFILM PITCH. Two constraints fix it.
`check_lattice_congruence` needs the CSG cylinder inside the cube of side
n*pitch. And at a real biofilm pitch (10 um, n=20) the domain is 200 um across,
where a 661.7 keV photon leaves without interacting, heating is numerically
zero, and every downstream invariant becomes untestable. So the pitch is
declared at DOSIMETRY scale for the same reason A0 is a 15 cm phantom, and is
labelled here as what it is: a number chosen to make the software path
exercisable, not an estimate of anything biological.

Usage:

    python calibration/scripts/emit_synthetic_reference_config.py --check-refusal
    python calibration/scripts/emit_synthetic_reference_config.py \
        --out config/reference_synthetic_biofilm_e2e.toml
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from biofilm_calibration.integration import VoxelBinding
from biofilm_calibration.materials.export import MaterialClass
from biofilm_calibration.reference_system import (ReferenceSystemSpec,
                                                  emit_transport_config)

REFERENCE_SYSTEM = "synthetic_biofilm_e2e"
SOURCE_ID = "SYNTHETIC_E2E_DECL"

# The snapshot this config is run against: export_checkpoint.jl's demo export
# hard-codes CPMParams(N = 20). The geometry checks are meaningless without it.
LATTICE_N = 20
PITCH_UM = 12000.0                      # 1.2 cm — see the module docstring
PITCH_CM = PITCH_UM * 1e-4
SIDE_CM = LATTICE_N * PITCH_CM          # 24.0 cm, roughly 2 attenuation lengths

# INVENTED. Liquid-water stoichiometry is used for the medium because it is the
# one composition in this file that is not arbitrary; the biomass reuses it so
# that a dose difference between biomass and medium comes only from density,
# which makes the ordering invariant readable.
WATER = {"H": 0.111894, "O": 0.888106}


def synthetic_binding() -> VoxelBinding:
    """The same coherent binding the repository declares, with the two blanks
    filled by invented numbers. The declared row in voxel_binding.csv keeps its
    blanks — this is a separate object, not an edit to the ledger."""
    return VoxelBinding(
        entity_kind="biomass_parcel",
        lineage_semantics="computational",
        material_class="baseline_biomass",
        material_represents="hydrated_bulk_biofilm",
        material_model_kind="hydrated_effective_medium",
        density_basis="wet_bulk",
        lattice_pitch_um=PITCH_UM,      # INVENTED
        density_g_cm3=1.05,             # INVENTED
        notes="synthetic: both numbers invented to exercise the software path",
    )


def _material(name, density, elements):
    # evidence_basis="synthetic" is REFUSED by materials/export.py for a
    # calibrated export. That refusal is correct and is not being weakened:
    # the emitter skips the evidence group only because evidence_policy says
    # synthetic, and every structural check still runs.
    return MaterialClass(name=name, density_g_cm3=density, elements=elements,
                         evidence_basis="synthetic", source_id=SOURCE_ID,
                         system_conditions="invented; no growth medium, no "
                                           "temperature, no hydration state")


def synthetic_spec(**overrides) -> ReferenceSystemSpec:
    spec = ReferenceSystemSpec(
        reference_system_id=REFERENCE_SYSTEM,
        system_provenance="declared",
        evidence_policy="synthetic",
        execution_class="synthetic_validation",
        target_calibration=False,
        lattice_n=LATTICE_N,
        origin_cm=(0.0, 0.0, 0.0),
        # 11.0 <= side/2 = 12.0, and 11.0 + 0.5 membrane = 11.5 still inside
        # the tallied cube, so no membrane dose falls outside the mesh.
        cylinder_radius_cm=11.0,
        cylinder_length_cm=24.0,
        membrane_thickness_cm=0.5,
        axial_bc="vacuum",
        radial_outer_bc="vacuum",
        source_spatial="point_origin",
        source_angular="isotropic",
        # Half a pitch off every lattice plane. The derived centre would be
        # (12.0, 12.0, 12.0) = 10 planes at a 1.2 cm pitch, i.e. on a surface
        # in all three dimensions at once.
        source_position_cm=(12.6, 12.6, 12.6),
        # Cs-137, the one line A0 retains. Evaluated EMISSION energy.
        spectrum_energies_eV=(661655.0,),
        spectrum_probabilities=(1.0,),
        medium=_material("invented medium (water stoichiometry)", 1.0, WATER),
        membrane=_material("invented membrane", 0.94, {"C": 0.856, "H": 0.144}),
        biomass_elements=WATER,
        # A0's decided setting, reused. This is NOT a biofilm resolution study:
        # that study needs a measured pitch and is a separate threshold.
        batches=80,
        particles=50000,
        seed=1,
        mesh_coarsening_factor=2,       # divides 20 -> a 10^3 tally mesh
        notes="synthetic reference system; every physical value is invented; "
              "calibrates nothing",
    )
    return dataclasses.replace(spec, **overrides) if overrides else spec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path,
                    help="write the transport-stage config here; the dosimetry "
                         "variant goes beside it as *.dosimetry.toml")
    ap.add_argument("--check-refusal", action="store_true",
                    help="blank the two invented numbers and show that emission "
                         "refuses, naming every missing field")
    ap.add_argument("--photons-per-second", type=float, default=1.0e12,
                    help="INVENTED source rate for the dosimetry stage")
    args = ap.parse_args(argv)

    if args.check_refusal:
        blank = dataclasses.replace(synthetic_binding(),
                                    lattice_pitch_um=None, density_g_cm3=None)
        try:
            emit_transport_config(blank, synthetic_spec())
        except ValueError as exc:
            print(exc)
            return 0
        print("EMISSION SUCCEEDED WITH BLANK NUMBERS — the refusal is broken.",
              file=sys.stderr)
        return 1

    if not args.out:
        ap.error("--out is required unless --check-refusal is given")

    binding = synthetic_binding()
    transport = emit_transport_config(binding, synthetic_spec(),
                                      stage="transport")
    dosimetry = emit_transport_config(
        binding, synthetic_spec(photons_per_second=args.photons_per_second),
        stage="dosimetry")

    args.out.write_text(transport, encoding="utf-8")
    dose_path = args.out.with_suffix(".dosimetry.toml")
    dose_path.write_text(dosimetry, encoding="utf-8")
    print(f"wrote {args.out}\nwrote {dose_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
