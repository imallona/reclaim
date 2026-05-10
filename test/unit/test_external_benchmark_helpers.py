"""Unit tests for the small external-benchmark helpers.

Covers:
  - merge_tecount_tables.parse_cnttable: header skip, integer cast, missing
    samples handled by caller.
  - merge_tecount_tables.merge: feature universe and per-sample lookup.
  - merge_tecount_tables.main CLI: end-to-end integer rendering.
  - squire_combine_counts.parse_te_name: 'subfamily:family:class' parsing.
  - squire_combine_counts CLI: combines two per-sample squire_count_*.txt
    files into a single TSV with the columns harmonize_external_counts
    expects.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = REPO_ROOT / 'workflow' / 'scripts'
sys.path.insert(0, str(SCRIPT_DIR))

import merge_tecount_tables as mtt  # noqa: E402
import squire_combine_counts as scc  # noqa: E402


def write(path, lines):
    with open(path, 'w') as fh:
        for line in lines:
            fh.write('\t'.join(str(x) for x in line) + '\n')


def test_parse_cnttable_skips_header_and_blanks(tmp_path):
    p = tmp_path / 'cnt.txt'
    write(p, [
        ('gene/TE', 'sample.bam'),
        ('AluY:Alu:SINE', 12),
        ('L1HS:L1:LINE', 0),
        ('', ''),
        ('MIR:MIR:SINE', 'NA'),
        ('SVA:SVA:Retroposon', 8),
    ])
    out = mtt.parse_cnttable(str(p))
    assert out == {'AluY:Alu:SINE': 12, 'L1HS:L1:LINE': 0,
                   'SVA:SVA:Retroposon': 8}


def test_merge_assembles_feature_universe():
    per_sample = {
        's1': {'AluY:Alu:SINE': 10, 'L1HS:L1:LINE': 0},
        's2': {'AluY:Alu:SINE': 8, 'MIR:MIR:SINE': 3},
    }
    feats, indexed = mtt.merge(per_sample)
    assert feats == ['AluY:Alu:SINE', 'L1HS:L1:LINE', 'MIR:MIR:SINE']
    assert indexed['AluY:Alu:SINE'] == {'s1': 10, 's2': 8}
    assert 'L1HS:L1:LINE' not in indexed  # both samples zero
    assert indexed['MIR:MIR:SINE'] == {'s2': 3}


def test_merge_rejects_empty_input():
    with pytest.raises(AssertionError):
        mtt.merge({})


def test_merge_cli_end_to_end(tmp_path):
    p1 = tmp_path / 's1.cntTable'
    p2 = tmp_path / 's2.cntTable'
    write(p1, [('gene/TE', 's1'), ('AluY:Alu:SINE', 10), ('L1HS:L1:LINE', 4)])
    write(p2, [('gene/TE', 's2'), ('AluY:Alu:SINE', 8), ('MIR:MIR:SINE', 3)])
    out = tmp_path / 'merged.tsv'
    script = SCRIPT_DIR / 'merge_tecount_tables.py'
    res = subprocess.run([
        sys.executable, str(script),
        '--inputs', str(p1), str(p2),
        '--samples', 's1', 's2',
        '--output', str(out),
    ], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    with open(out) as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh]
    assert rows[0] == ['feature_id', 's1', 's2']
    by_feat = {r[0]: r[1:] for r in rows[1:]}
    assert by_feat['AluY:Alu:SINE'] == ['10', '8']
    assert by_feat['L1HS:L1:LINE'] == ['4', '0']
    assert by_feat['MIR:MIR:SINE'] == ['0', '3']


def test_merge_cli_rejects_length_mismatch(tmp_path):
    p1 = tmp_path / 's1.cntTable'
    write(p1, [('gene/TE', 's1'), ('AluY:Alu:SINE', 5)])
    out = tmp_path / 'merged.tsv'
    script = SCRIPT_DIR / 'merge_tecount_tables.py'
    res = subprocess.run([
        sys.executable, str(script),
        '--inputs', str(p1),
        '--samples', 's1', 's2_extra',
        '--output', str(out),
    ], capture_output=True, text=True)
    assert res.returncode != 0
    assert 'same length' in res.stderr


def test_parse_te_name_full():
    assert scc.parse_te_name('AluY:Alu:SINE') == ('AluY', 'Alu', 'SINE')


def test_parse_te_name_underspecified():
    assert scc.parse_te_name('AluY') == ('AluY', '', '')


def test_squire_combine_cli_end_to_end(tmp_path):
    sample_a = tmp_path / 'squire_count_s1.txt'
    sample_b = tmp_path / 'squire_count_s2.txt'
    cols = ['TE_ID', 'TE_chr', 'TE_start', 'TE_stop', 'TE_name', 'milliDiv',
            'strand', 'uniq_counts', 'tot_counts', 'tot_reads', 'score', 'fpkm']
    rows_a = [
        cols,
        ['locus_aluy_1', 'chr1', '100', '300', 'AluY:Alu:SINE',
         '0.05', '+', '5', '7', '7', '1.0', '10.0'],
        ['locus_l1hs_1', 'chr1', '500', '6000', 'L1HS:L1:LINE',
         '0.03', '+', '2', '3', '3', '1.0', '5.0'],
    ]
    rows_b = [
        cols,
        ['locus_aluy_1', 'chr1', '100', '300', 'AluY:Alu:SINE',
         '0.05', '+', '4', '5', '5', '1.0', '8.0'],
    ]
    write(sample_a, rows_a)
    write(sample_b, rows_b)
    out = tmp_path / 'combined.tsv'
    script = SCRIPT_DIR / 'squire_combine_counts.py'
    res = subprocess.run([
        sys.executable, str(script),
        '--indir', str(tmp_path),
        '--output', str(out),
    ], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    with open(out) as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh]
    header = rows[0]
    assert header[:5] == ['TE_ID', 'TE_name', 'family', 'class', 'subfamily']
    assert sorted(header[5:]) == ['s1', 's2']
    by_te = {r[0]: dict(zip(header, r)) for r in rows[1:]}
    aluy = by_te['locus_aluy_1']
    assert aluy['family'] == 'Alu'
    assert aluy['class'] == 'SINE'
    assert aluy['subfamily'] == 'AluY'
    assert float(aluy['s1']) == 5.0
    assert float(aluy['s2']) == 4.0
    assert 'locus_l1hs_1' in by_te
    assert float(by_te['locus_l1hs_1']['s2']) == 0.0


def test_squire_combine_cli_rejects_empty_dir(tmp_path):
    out = tmp_path / 'combined.tsv'
    script = SCRIPT_DIR / 'squire_combine_counts.py'
    res = subprocess.run([
        sys.executable, str(script),
        '--indir', str(tmp_path),
        '--output', str(out),
    ], capture_output=True, text=True)
    assert res.returncode != 0
    assert 'no squire_count' in res.stderr
