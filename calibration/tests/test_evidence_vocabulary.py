"""The calibration schema reads the contract's evidence vocabulary as ONE object.

The set itself and the two-relation controls are tested in
contract/tests/test_evidence_basis.py. This file exists because that tier of CI
does not install this package: an `is` check on an object it cannot import
would be a skip, and a skip is uncovered surface.
"""

from physical_contract import PARAMETER_EVIDENCE_BASIS

from biofilm_calibration.materials.export import EXPORTABLE_EVIDENCE
from biofilm_calibration.schema import EVIDENCE_BASIS, STATUS


def test_schema_reads_the_contract_object_not_a_copy():
    assert EVIDENCE_BASIS is PARAMETER_EVIDENCE_BASIS
    assert frozenset(set(PARAMETER_EVIDENCE_BASIS)) is not PARAMETER_EVIDENCE_BASIS  # a copy would pass ==


def test_export_policy_is_a_subset_and_unresolved_is_a_status():
    assert EXPORTABLE_EVIDENCE <= PARAMETER_EVIDENCE_BASIS
    assert "unresolved" in STATUS and "unresolved" not in EVIDENCE_BASIS
