"""CPM biofilm <-> OpenMC photon-transport coupling (offline-first).

Layering: config -> snapshot -> materials -> model -> dose -> lineage.
`openmc` is imported only inside model/transport functions, so every
schema/postprocessing path runs in a bare numpy+h5py environment.
"""

from .config import (BIOFILM_CYLINDER, MODEL_KINDS, STAGES, WATER_PHANTOM,
                     Config, ConfigError, CPMFeedbackConfig, DoseRateConfig,
                     MembraneFeedbackConfig, TransportConfig, load_config,
                     load_cpm_feedback_config, load_dose_rate_config,
                     load_membrane_feedback_config, load_transport_config,
                     required_keys, required_material_classes)
from .snapshot import Snapshot, load_snapshot, to_openmc_lattice_order, from_openmc_lattice_order
