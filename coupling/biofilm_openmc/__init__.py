"""CPM biofilm <-> OpenMC photon-transport coupling (offline-first).

Layering: config -> snapshot -> materials -> model -> dose -> lineage.
`openmc` is imported only inside model/transport functions, so every
schema/postprocessing path runs in a bare numpy+h5py environment.
"""

from .config import Config, ConfigError, load_config
from .snapshot import Snapshot, load_snapshot, to_openmc_lattice_order, from_openmc_lattice_order
