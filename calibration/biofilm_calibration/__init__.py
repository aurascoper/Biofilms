"""Calibration contracts and analysis harnesses for the CPM biofilm model.

Two independent workstreams, deliberately not sharing state:

    spatial/    what a CPM lattice site and a cell ID physically represent
    materials/  what material occupies an "occupied" voxel, in what units

Neither fits parameters. There is no measured biological data in this
repository, so these establish what would have to be measured, what the current
model can represent, and where the two meet — and refuse to invent values.
"""

from .schema import (EVIDENCE_BASIS, STATUS, Column, SchemaError, TableSchema,
                     read_table, write_template)

__all__ = ["Column", "TableSchema", "SchemaError", "read_table",
           "write_template", "STATUS", "EVIDENCE_BASIS"]
