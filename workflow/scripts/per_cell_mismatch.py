#!/usr/bin/env python3
"""
Per-cell mismatch summary from a STARsolo chromium BAM.

For each cell barcode listed in the perturbation table (filtered to the given
gem_group), sums aligned bases and NM tags from primary, non-supplementary
alignments. Writes one row per kept cell:

    sample_id, gem_group, cell_barcode,
    n_reads, n_aligned_bases, sum_nm, mismatch_rate

mismatch_rate is sum_nm / n_aligned_bases. NM (per the SAM spec) includes
substitutions plus insertions and deletions, so this rate is mismatch+indel,
not pure substitution.
"""

import argparse
import csv
import sys
from collections import defaultdict

import pysam


def load_cells(cells_tsv, gem_group):
    keep = set()
    with open(cells_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if int(row['gem_group']) != gem_group:
                continue
            cb = row['cell_barcode'].split('-')[0]
            keep.add(cb)
    return keep


def walk(bam_path, keep_cb):
    counts = defaultdict(lambda: [0, 0, 0])
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            try:
                cb = r.get_tag('CB')
            except KeyError:
                continue
            if cb not in keep_cb:
                continue
            try:
                nm = r.get_tag('NM')
            except KeyError:
                continue
            row = counts[cb]
            row[0] += 1
            row[1] += r.query_alignment_length
            row[2] += nm
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--cells', required=True,
                    help='cells_to_perturbation.tsv')
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--gem-group', type=int, required=True)
    ap.add_argument('--min-reads', type=int, default=50,
                    help='skip cells with fewer reads than this (default 50)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    keep = load_cells(args.cells, args.gem_group)
    if not keep:
        sys.exit(f'no cells for gem_group {args.gem_group} in {args.cells}')

    counts = walk(args.bam, keep)

    cols = ['sample_id', 'gem_group', 'cell_barcode',
            'n_reads', 'n_aligned_bases', 'sum_nm', 'mismatch_rate']
    with open(args.out, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for cb, (n, bases, nm) in sorted(counts.items()):
            if n < args.min_reads or bases == 0:
                continue
            rate = nm / bases
            out.write(f'{args.sample_id}\t{args.gem_group}\t{cb}\t'
                      f'{n}\t{bases}\t{nm}\t{rate:.6f}\n')


if __name__ == '__main__':
    main()
