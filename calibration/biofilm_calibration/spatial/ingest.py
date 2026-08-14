"""Format-agnostic ingestion of a biomass field.

The package stays numpy-only. Readers for ND2, CZI, LIF or TIFF live behind an
optional `dataset` extra that the core never imports, so the harness remains
testable without a multi-gigabyte file and without an imaging stack. What this
module defines is the CONTRACT every reader must satisfy, which is the part
that matters scientifically:

    an array, its physical voxel calibration, and a declared statement of what
    the mask encloses — validated together, or refused.

WHY THE THREE TRAVEL TOGETHER. A biomass array on its own is dimensionless
pixels. With a voxel size it becomes a volume. Only with a segmentation basis
does it become a volume OF something, and it is that third piece which decides
whether the volume may be divided into a measured mass. Conversion pipelines
routinely drop the voxel calibration and never carried the basis at all, so
both are required here rather than looked up later.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..acquisition import SEGMENTATION_BASIS
from .field import FieldError, check_field


@dataclass(frozen=True)
class BiomassField:
    """A validated biomass field with everything needed to interpret it."""
    values: np.ndarray                  # 3-D, in [0, 1]
    voxel_size_um: tuple                # (x, y, z)
    segmentation_basis: str
    sample_id: str
    source_id: str
    reader: str = "in_memory"
    notes: str = ""

    @property
    def voxel_volume_um3(self) -> float:
        return float(np.prod(np.asarray(self.voxel_size_um, dtype=float)))

    @property
    def shape(self) -> tuple:
        return tuple(self.values.shape)

    def extent_um(self) -> tuple:
        return tuple(float(n * v) for n, v in
                     zip(self.values.shape, self.voxel_size_um))

    def hydrated_volume_cm3(self):
        """V_hydrated = integral of B dV, tagged with what it encloses.

        Returns a `materials.basis_conversion.HydratedVolume`, so the volume
        arrives at the density calculation already carrying its basis and
        cannot be divided into a mismatched mass.
        """
        from ..materials.basis_conversion import HydratedVolume

        um3 = float(self.values.sum()) * self.voxel_volume_um3
        return HydratedVolume(um3 * 1e-12, self.segmentation_basis)  # um3 -> cm3


def ingest(values, voxel_size_um, segmentation_basis: str, *,
           sample_id: str, source_id: str, reader: str = "in_memory",
           notes: str = "") -> BiomassField:
    """Validate an array, its calibration and its basis together."""
    arr = check_field(values)

    v = np.asarray(voxel_size_um, dtype=float)
    if v.shape == ():
        v = np.repeat(v, 3)
    if v.shape != (3,):
        raise FieldError(
            f"voxel_size_um must be three values (x, y, z), got {voxel_size_um!r} "
            "— anisotropic z-spacing is the norm in confocal stacks and must not "
            "be collapsed to a scalar by accident")
    if np.any(v <= 0) or not np.isfinite(v).all():
        raise FieldError(f"voxel_size_um must be positive and finite, got {list(v)}")

    if segmentation_basis not in SEGMENTATION_BASIS:
        raise FieldError(
            f"unknown segmentation basis {segmentation_basis!r}; expected one "
            f"of {sorted(SEGMENTATION_BASIS)}. What the mask encloses decides "
            "whether its volume may be divided into a measured mass.")
    if not sample_id.strip():
        raise FieldError("sample_id is required — an unattributed field cannot "
                         "be paired with a gravimetry sibling")
    if not source_id.strip():
        raise FieldError("source_id is required")

    return BiomassField(values=arr, voxel_size_um=tuple(float(x) for x in v),
                        segmentation_basis=segmentation_basis,
                        sample_id=sample_id, source_id=source_id,
                        reader=reader, notes=notes)


def available_readers() -> dict:
    """Which optional readers are importable, without importing them into the
    core. Absence is reported, never raised: the harness works without any."""
    out = {}
    for name, module in (("nd2", "nd2"), ("tiff", "tifffile")):
        try:
            __import__(module)
            out[name] = True
        except ImportError:
            out[name] = False
    return out
