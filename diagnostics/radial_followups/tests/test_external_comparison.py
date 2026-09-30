import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location('compare_external', Path(__file__).parents[1] / 'compare_external.py')
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_supplied_artifact_mutants_are_rejected(tmp_path):
    path = Path(__file__).parents[3] / 'diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv'
    rows = comparison.read(path)
    key = next(k for k in rows if k[1] == 'full' and k[2] == 100)
    changed = {k: dict(v) for k, v in rows.items()}
    changed[key]['mean_r_over_R'] = str(float(changed[key]['mean_r_over_R']) + 0.01)
    assert changed[key] != rows[key]
    with pytest.raises(ValueError, match='printed fields differ'):
        comparison.compare_rows(changed, rows)
    del changed[key]
    with pytest.raises(ValueError, match='coverage differs'):
        comparison.compare_rows(changed, rows)
    lines = path.read_text().splitlines()
    row = next(x for x in lines if x.startswith('42,full,0,1,'))
    duplicate = tmp_path / 'duplicate.csv'
    duplicate.write_text('\n'.join(lines + [row]) + '\n')
    with pytest.raises(ValueError, match='duplicate row'):
        comparison.read(duplicate)
