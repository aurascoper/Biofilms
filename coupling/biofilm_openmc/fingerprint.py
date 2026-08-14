"""Three separate state hashes and the EXACT transport-rerun policy.

- transport_state_hash: everything transport physically depends on —
  voxel material classes, compositions/densities, geometry, source,
  transport settings, nuclear-data identity. Any change -> rerun.
- label_state_hash: cell/species/lineage/generation arrays. Changing ONLY
  labels reuses the cached field and re-attributes.
- (restart_state_hash lives Julia-side with the restart checkpoint.)

No approximate voxel-churn threshold (review amendment 7): one changed
high-Z voxel near the source can matter more than thousands of unchanged
low-Z ones. Threshold-based skipping is a later, *measured* approximation.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np

from .config import Config
from .materials import voxel_class_array
from .snapshot import Snapshot


def _sha(*chunks: bytes) -> str:
    h = hashlib.sha256()
    for c in chunks:
        h.update(c)
    return h.hexdigest()


def transport_state_hash(snapshot: Snapshot, config: Config,
                         nuclear_data_id: str = "unset") -> str:
    classes = voxel_class_array(snapshot, config)
    manifest = {
        "materials": {name: {"density_g_cm3": mc.density_g_cm3,
                             "elements": list(mc.elements)}
                      for name, mc in sorted(config.materials.items())},
        "geometry": [config.voxel_pitch_cm, list(config.origin_cm),
                     config.cylinder_radius_cm, config.cylinder_length_cm,
                     config.membrane_thickness_cm,
                     config.axial_bc, config.radial_outer_bc],
        "source": [config.photons_per_second,
                   list(config.spectrum_energies_eV),
                   list(config.spectrum_probabilities),
                   config.source_spatial, config.source_angular],
        "transport": [config.batches, config.particles, config.seed],
        "nuclear_data": nuclear_data_id,
        "mesh": list(snapshot.cell_id.shape),
    }
    class_bytes = "\x00".join(classes.ravel(order="C")).encode()
    return _sha(class_bytes, json.dumps(manifest, sort_keys=True).encode())


def label_state_hash(snapshot: Snapshot) -> str:
    return _sha(*(np.ascontiguousarray(a).tobytes() for a in (
        snapshot.cell_id, snapshot.species_id,
        snapshot.lineage_id, snapshot.generation)))


def needs_rerun(previous_transport_hash: str | None, snapshot: Snapshot,
                config: Config, nuclear_data_id: str = "unset") -> bool:
    """EXACT policy: rerun on any transport-hash change; reuse only on an
    exact match (then only label attribution needs recomputing)."""
    if previous_transport_hash is None:
        return True
    return transport_state_hash(snapshot, config, nuclear_data_id) != \
        previous_transport_hash
