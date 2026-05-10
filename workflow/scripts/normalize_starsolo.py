#!/usr/bin/env python3
"""
Convert a STARsolo Solo.out/Gene/raw directory to the standard
feature x cell TSV used by evaluate.py.

Supports multi-level aggregation for repeat feature sets via --locus-map
and --granularity:
  gene_id   - output as-is (STARsolo native gene_id level)
  family_id - aggregate gene_ids to RepeatMasker family using locus map
  class_id  - aggregate gene_ids to RepeatMasker class using locus map

The locus map TSV has columns (no header): transcript_id, gene_id, family_id, class_id.
gene_id -> group is derived by scanning column 1 -> target column (unique pairs).
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
import scipy.io


def read_lines(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return [line.rstrip('\n') for line in fh]


def load_starsolo_raw(raw_dir, mtx_name='matrix.mtx'):
    barcodes = read_lines(os.path.join(raw_dir, 'barcodes.tsv'))
    feature_lines = read_lines(os.path.join(raw_dir, 'features.tsv'))
    feature_ids = [line.split('\t')[0] for line in feature_lines]
    mat = scipy.io.mmread(os.path.join(raw_dir, mtx_name)).tocsc()
    return barcodes, feature_ids, mat


def build_gene_to_group(locus_map_path, granularity):
    """Build dict: gene_id -> group_id from locus map (columns 1 and 2 or 3)."""
    col = {'gene_id': 1, 'family_id': 2, 'class_id': 3}[granularity]
    mapping = {}
    with open(locus_map_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene_id = parts[1]
            group = parts[col]
            mapping[gene_id] = group
    return mapping


def _format_value(v):
    """Render a count as integer when whole (the common case) and as a
    short float when fractional. sc_count_features.py emits 1/NH weights
    in --multimapper multi mode, so silent int-truncation here would
    invalidate tagcount-multi evaluation against ground truth."""
    if v == 0:
        return '0'
    iv = int(v)
    if iv == v:
        return str(iv)
    return f'{v:.6f}'


def write_count_matrix(barcodes, feature_ids, sparse_mat, output_path, gene_to_group=None):
    if gene_to_group is None:
        # gene_id passthrough
        with open(output_path, 'w') as fh:
            fh.write('feature_id\t' + '\t'.join(barcodes) + '\n')
            for feat_idx, feat_id in enumerate(feature_ids):
                row = sparse_mat.getrow(feat_idx)
                if row.nnz == 0:
                    continue
                values = row.toarray().flatten()
                fh.write(feat_id + '\t' + '\t'.join(_format_value(v) for v in values) + '\n')
    else:
        # aggregate by group
        n_cells = len(barcodes)
        group_counts = defaultdict(lambda: [0.0] * n_cells)
        skipped = 0
        for feat_idx, feat_id in enumerate(feature_ids):
            group = gene_to_group.get(feat_id)
            if group is None:
                skipped += 1
                continue
            row = sparse_mat.getrow(feat_idx)
            if row.nnz == 0:
                continue
            values = row.toarray().flatten()
            for ci, v in enumerate(values):
                group_counts[group][ci] += float(v)
        if skipped:
            print('{} gene_ids had no group mapping (skipped)'.format(skipped), file=sys.stderr)
        with open(output_path, 'w') as fh:
            fh.write('feature_id\t' + '\t'.join(barcodes) + '\n')
            for group in sorted(group_counts):
                vals = group_counts[group]
                if all(v == 0 for v in vals):
                    continue
                fh.write(group + '\t' + '\t'.join(_format_value(v) for v in vals) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-dir', required=True, help='Solo.out/Gene/raw directory')
    ap.add_argument('--output', required=True, help='Output TSV path')
    ap.add_argument('--locus-map', default=None,
                    help='TSV: transcript_id, gene_id, family_id, class_id (no header)')
    ap.add_argument('--multimapper-mode', default='unique',
                    choices=['unique', 'multi'])
    ap.add_argument('--granularity', default='gene_id',
                    choices=['gene_id', 'family_id', 'class_id'])
    args = ap.parse_args()

    em_mtx = os.path.join(args.raw_dir, 'UniqueAndMult-EM.mtx')
    mtx_name = 'UniqueAndMult-EM.mtx' if args.multimapper_mode == 'multi' and os.path.exists(em_mtx) else 'matrix.mtx'
    print('Using matrix: {}'.format(mtx_name), file=sys.stderr)
    barcodes, feature_ids, mat = load_starsolo_raw(args.raw_dir, mtx_name)
    print('{} features, {} barcodes'.format(len(feature_ids), len(barcodes)), file=sys.stderr)

    gene_to_group = None
    if args.locus_map and args.granularity != 'gene_id':
        print('Building gene_id -> {} map from {}'.format(args.granularity, args.locus_map),
              file=sys.stderr)
        gene_to_group = build_gene_to_group(args.locus_map, args.granularity)
        print('  {} gene entries loaded'.format(len(gene_to_group)), file=sys.stderr)

    write_count_matrix(barcodes, feature_ids, mat, args.output, gene_to_group)
    print('Written to {}'.format(args.output), file=sys.stderr)


if __name__ == '__main__':
    main()
