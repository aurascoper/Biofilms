"""Compare supplied cloud CSVs with the recorded local shards at printed precision.

This does not authenticate where the supplied files were computed, infer full-state
identity, or replace either study's prediction. No attached program is executed.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import statistics

FACTORS = (1, 10, 100, 1000, 10000, 100000, 1000000)
SEEDS = range(42, 298)
COLUMNS = ('seed', 'arm', 'mcs', 'species', 'species_name', 'ncells', 'volume',
           'mean_r', 'mean_r_over_R', 'mean_melanin')


def read(path):
    reader = csv.DictReader(x for x in path.read_text().splitlines() if x and not x.startswith('#'))
    if tuple(reader.fieldnames or ()) != COLUMNS:
        raise ValueError(f'wrong schema: {path}')
    rows = list(reader)
    keys = [(int(r['seed']), r['arm'], int(r['mcs']), int(r['species'])) for r in rows]
    if len(set(keys)) != len(rows):
        raise ValueError(f'duplicate row: {path}')
    return dict(zip(keys, rows))


def compare_rows(actual, expected):
    if set(actual) != set(expected):
        raise ValueError('row-key coverage differs')
    differences = []
    for key, row in actual.items():
        for column in COLUMNS:
            if row[column] != expected[key][column]:
                differences.append((key, column, row[column], expected[key][column]))
    if differences:
        raise ValueError(f'{len(differences)} printed fields differ; first: {differences[0]}')
    return len(actual), len(actual) * len(COLUMNS)


def compare(downloads, local):
    result = {'scope': 'All supplied rows and printed fields; six-decimal floating outputs, not unrounded states.',
              'external_runtime': '1.12.7 (supplied metadata)', 'local_runtime': '1.12.6',
              'studies': {}, 'inputs': {}}
    for study in ('beta', 'alpha', 'horizon'):
        path = downloads / f'followup_{study}_n40_p6_seeds42-297.csv'
        actual = read(path)
        arms = ('full', 'null') if study == 'horizon' else ('full', 'null', *(f'{study}_x{f}' for f in FACTORS[1:]))
        times = (0, 100, 200, 300, 400) if study == 'horizon' else (0, 100)
        keys = set(itertools.product(SEEDS, arms, times, range(1, 8)))
        if set(actual) != keys:
            raise ValueError(f'incomplete supplied {study} study')
        expected = {}
        for seed, arm in itertools.product(SEEDS, arms):
            if study == 'horizon':
                shard = local / 'horizon' / arm / f'seed_{seed}.csv'
            elif arm == 'null':
                shard = local / 'baseline/null' / f'seed_{seed}.csv'
            else:
                factor = 1 if arm == 'full' else int(arm.split('_x')[1])
                shard = local / study / f'factor_{factor}' / f'seed_{seed}.csv'
            for (_, _, mcs, species), row in read(shard).items():
                row = {**row, 'arm': arm}
                expected[(seed, arm, mcs, species)] = row
        common_times = (0, 100, 400) if study == 'horizon' else times
        common_keys = set(itertools.product(SEEDS, arms, common_times, range(1, 8)))
        if set(expected) != common_keys:
            raise ValueError(f'incomplete local {study} study')
        comparable = {k: r for k, r in actual.items() if k[2] in common_times}
        n, fields = compare_rows(comparable, expected)
        result['studies'][study] = {'rows': n, 'fields_equal': fields, 'exact_at_printed_precision': True,
                                   'compared_mcs': common_times, 'external_only_rows': len(actual)-n,
                                   'external_only_mcs': [t for t in times if t not in common_times]}
        result['inputs'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        if study == 'alpha':
            summaries = []
            for factor, species in itertools.product(FACTORS, (1, 3)):
                arm = 'full' if factor == 1 else f'alpha_x{factor}'
                rows = [actual[(seed, arm, 100, species)] for seed in SEEDS]
                shifts = [float(r['mean_r_over_R']) - float(actual[(int(r['seed']), 'null', 100, species)]['mean_r_over_R']) for r in rows]
                melanin = [float(r['mean_melanin']) for r in rows]
                summaries.append({'factor': factor, 'species': species, 'mean_shift': statistics.mean(shifts),
                                  'positive_seeds': sum(x > 0 for x in shifts),
                                  'mean_melanin': statistics.mean(melanin), 'min_melanin': min(melanin), 'max_melanin': max(melanin)})
            result['alpha_summary'] = summaries
    result['status'] = 'passed'
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('downloads', type=Path)
    parser.add_argument('local_archive', type=Path)
    parser.add_argument('receipt', type=Path)
    args = parser.parse_args()
    record = compare(args.downloads, args.local_archive)
    args.receipt.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record['studies'], indent=2))
