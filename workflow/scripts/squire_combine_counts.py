#!/usr/bin/env python3
"""Combine SQuIRE Count per-sample outputs into a single TSV consumed by
harmonize_external_counts.py --tool squire.

SQuIRE Count writes per-sample squire_count_<sample>.txt files with
columns:
  TE_ID, TE_chr, TE_start, TE_stop, TE_name, milliDiv, strand,
  uniq_counts, tot_counts, tot_reads, score, fpkm

Our combined output retains the locus identity and assigns subfamily,
family, and class via the TE_name field, then adds one count column per
sample (uniq_counts is preferred to match REclaim's STAR unique mode).

This script is intentionally minimal so the pipeline can run on machines
where SQuIRE itself does not install cleanly: as long as the .txt files
are present, the combine step succeeds.

Asserts:
  - At least one squire_count_*.txt file exists in --indir.
  - Every per-sample table has the required columns.
  - The final TSV has the columns the harmoniser expects:
    TE_ID, TE_name, family, class, subfamily, <sample columns>.
"""
import argparse
import csv
import os
from collections import defaultdict


REQUIRED_COLS = ('TE_ID', 'TE_name', 'uniq_counts')


def parse_te_name(te_name):
    """SQuIRE TE_name is 'subfamily:family:class' (e.g. AluY:Alu:SINE)."""
    parts = te_name.split(':')
    if len(parts) < 3:
        return te_name, '', ''
    return parts[0], parts[1], parts[2]


def parse_squire_count_file(path, sample):
    out = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        missing = [c for c in REQUIRED_COLS if c not in reader.fieldnames]
        assert not missing, (
            f'{path!r} missing required columns {missing}')
        for row in reader:
            out.append({
                'TE_ID': row['TE_ID'],
                'TE_name': row['TE_name'],
                'sample': sample,
                'count': float(row['uniq_counts'] or 0),
            })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', required=True,
                    help='Directory holding squire_count_<sample>.txt files')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    files = sorted(
        os.path.join(args.indir, f) for f in os.listdir(args.indir)
        if f.startswith('squire_count_') and f.endswith('.txt'))
    assert files, f'no squire_count_*.txt in {args.indir!r}'

    per_te = defaultdict(lambda: {'TE_name': '', 'samples': defaultdict(float)})
    sample_set = []
    for path in files:
        sample = os.path.basename(path)[len('squire_count_'):-len('.txt')]
        sample_set.append(sample)
        for r in parse_squire_count_file(path, sample):
            per_te[r['TE_ID']]['TE_name'] = r['TE_name']
            per_te[r['TE_ID']]['samples'][sample] += r['count']
    assert per_te, 'no rows merged from squire_count tables'

    with open(args.output, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        header = ['TE_ID', 'TE_name', 'family', 'class', 'subfamily'] + sample_set
        w.writerow(header)
        for te_id in sorted(per_te.keys()):
            te_name = per_te[te_id]['TE_name']
            sub, fam, cls = parse_te_name(te_name)
            row = [te_id, te_name, fam, cls, sub]
            for s in sample_set:
                row.append(per_te[te_id]['samples'].get(s, 0))
            w.writerow(row)


if __name__ == '__main__':
    main()
