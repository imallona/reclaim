#!/usr/bin/env python3
"""Merge TEcount per-sample cntTables into a single feature x sample TSV.

TEcount writes one cntTable per --project run with two columns:
  feature_id\\tcount

We invoke TEcount once per BAM (one per simulated cell for SmartSeq2),
then this helper joins the per-sample tables on feature_id and writes a
merged TSV with header `feature_id\\t<s1>\\t<s2>...`.

Asserts:
  - len(--inputs) == len(--samples)
  - every cntTable has at least one feature row
  - merged feature universe is non-empty
"""
import argparse
import csv
from collections import defaultdict


def parse_cnttable(path):
    """Return {feature_id: count}. Skips header rows that start with
    'gene/TE' or contain non-numeric counts."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            feat, val = parts[0], parts[1]
            if not feat or feat.lower().startswith(('gene/te', 'gene_id',
                                                    'feature_id')):
                continue
            try:
                count = int(float(val))
            except ValueError:
                continue
            out[feat] = count
    return out


def merge(per_sample):
    """Build a sorted feature universe and a feature -> {sample: count}
    mapping from a {sample: {feature: count}} dict."""
    assert per_sample, 'no per-sample tables to merge'
    feature_universe = set()
    for d in per_sample.values():
        feature_universe.update(d.keys())
    assert feature_universe, 'merged feature universe is empty'
    out = defaultdict(dict)
    for sample, d in per_sample.items():
        for feat, count in d.items():
            if count > 0:
                out[feat][sample] = count
    return sorted(feature_universe), {k: dict(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', nargs='+', required=True,
                    help='Per-sample cntTable paths in --samples order')
    ap.add_argument('--samples', nargs='+', required=True,
                    help='Sample IDs parallel to --inputs')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    assert len(args.inputs) == len(args.samples), (
        '--inputs and --samples must have the same length')
    per_sample = {}
    for path, sample in zip(args.inputs, args.samples):
        per_sample[sample] = parse_cnttable(path)
        assert per_sample[sample], (
            f'cntTable {path!r} has no feature rows; check TEcount output')

    features, indexed = merge(per_sample)
    with open(args.output, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['feature_id'] + list(args.samples))
        for f in features:
            row = [f] + [indexed.get(f, {}).get(s, 0) for s in args.samples]
            w.writerow(row)


if __name__ == '__main__':
    main()
