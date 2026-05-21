"""Unit tests for workflow/scripts/harmonize_external_counts.py.

Covers:
  - Pure helper map_features_to_reclaim: identity match, suffix-stripping
    match (TEtranscripts subfamily:family:class), no-match.
  - parse_locus_map: 4-col TSV roundtrip and conflict warning.
  - aggregate_to_granularity: gene_id passthrough, family_id rollup,
    class_id rollup, dropped-feature accounting.
  - parse_tetranscripts_table: header parsing and subfamily extraction.
  - parse_squire_count: subfamily aggregation across loci.
  - check_overlap: passes when overlap above threshold, asserts when below.
  - write_count_matrix: shape and ordering.

scTE parsing (parse_scte_h5ad) is NOT exercised here because anndata is not
a hard dependency of the test suite; it is integration-tested when the scte
conda env is available.
"""
import os
import sys
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = REPO_ROOT / 'workflow' / 'scripts'
sys.path.insert(0, str(SCRIPT_DIR))

import harmonize_external_counts as hec  # noqa: E402


def write_tsv(path, rows):
    with open(path, 'w') as fh:
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')


def test_map_features_identity():
    out = hec.map_features_to_reclaim(
        ['L1HS', 'AluY', 'novel'],
        ['L1HS', 'AluY', 'L1PA2'])
    assert out == {'L1HS': 'L1HS', 'AluY': 'AluY', 'novel': None}


def test_map_features_suffix_strip():
    out = hec.map_features_to_reclaim(
        ['L1HS:L1:LINE', 'AluY:Alu:SINE', 'rogue:foo:bar'],
        ['L1HS', 'AluY'],
        strip_suffix=':')
    assert out == {
        'L1HS:L1:LINE': 'L1HS',
        'AluY:Alu:SINE': 'AluY',
        'rogue:foo:bar': None}


def test_map_features_rejects_non_iter():
    with pytest.raises(AssertionError):
        hec.map_features_to_reclaim('AluY', ['AluY'])


def test_parse_locus_map_roundtrip(tmp_path):
    p = tmp_path / 'lm.tsv'
    write_tsv(p, [
        ('transcript_id', 'gene_id', 'family_id', 'class_id'),
        ('AluY_chr1_100', 'AluY', 'Alu', 'SINE'),
        ('AluY_chr2_200', 'AluY', 'Alu', 'SINE'),
        ('L1HS_chr3_300', 'L1HS', 'L1', 'LINE'),
    ])
    g2c, f2c, g2f = hec.parse_locus_map(str(p))
    assert g2c['AluY'] == 'SINE'
    assert g2c['L1HS'] == 'LINE'
    assert f2c['Alu'] == 'SINE'
    assert g2f['L1HS'] == 'L1'


def test_parse_locus_map_warns_on_conflict(tmp_path, capsys):
    p = tmp_path / 'lm_conflict.tsv'
    write_tsv(p, [
        ('transcript_id', 'gene_id', 'family_id', 'class_id'),
        ('X1', 'AluY', 'Alu', 'SINE'),
        ('X2', 'AluY', 'Alu', 'LINE'),  # same gene, different class
    ])
    hec.parse_locus_map(str(p))
    captured = capsys.readouterr()
    assert 'multiple class_ids' in captured.err


def test_aggregate_gene_passthrough():
    counts = {'AluY': {'s1': 10, 's2': 20}, 'L1HS': {'s1': 5}}
    out, dropped = hec.aggregate_to_granularity(counts, {}, {}, 'gene_id')
    assert out == counts
    assert dropped == 0


def test_aggregate_to_family():
    counts = {
        'AluY':   {'s1': 10, 's2': 20},
        'AluYa5': {'s1': 5,  's2': 10},
        'L1HS':   {'s1': 100},
    }
    g2f = {'AluY': 'Alu', 'AluYa5': 'Alu', 'L1HS': 'L1'}
    g2c = {'AluY': 'SINE', 'AluYa5': 'SINE', 'L1HS': 'LINE'}
    out, dropped = hec.aggregate_to_granularity(counts, g2f, g2c, 'family_id')
    assert out['Alu'] == {'s1': 15, 's2': 30}
    assert out['L1'] == {'s1': 100}
    assert dropped == 0


def test_aggregate_to_class_drops_unmapped():
    counts = {'AluY': {'s1': 10}, 'rogue': {'s1': 1}}
    g2f = {'AluY': 'Alu'}
    g2c = {'AluY': 'SINE'}
    out, dropped = hec.aggregate_to_granularity(counts, g2f, g2c, 'class_id')
    assert out == {'SINE': {'s1': 10}}
    assert dropped == 1


def test_parse_tetranscripts_table(tmp_path):
    p = tmp_path / 'te.tsv'
    write_tsv(p, [
        ('feature_id', 's1', 's2'),
        ('AluY:Alu:SINE', 12, 8),
        ('L1HS:L1:LINE', 0, 4),
        ('MIR:MIR:SINE', 3, 0),
    ])
    out = hec.parse_tetranscripts_table(str(p), ['s1', 's2'])
    assert out == {
        'AluY': {'s1': 12, 's2': 8},
        'L1HS': {'s2': 4},
        'MIR':  {'s1': 3},
    }


def test_parse_tetranscripts_missing_sample_raises(tmp_path):
    p = tmp_path / 'te.tsv'
    write_tsv(p, [
        ('feature_id', 's1'),
        ('AluY:Alu:SINE', 1),
    ])
    with pytest.raises(AssertionError):
        hec.parse_tetranscripts_table(str(p), ['s1', 's2_missing'])


def test_parse_squire_count_aggregates_across_loci(tmp_path):
    p = tmp_path / 'sq.tsv'
    write_tsv(p, [
        ('TE_ID', 'TE_name', 'family', 'class', 'subfamily', 's1', 's2'),
        ('locus1', 'AluY_locus1', 'Alu', 'SINE', 'AluY', 5, 2),
        ('locus2', 'AluY_locus2', 'Alu', 'SINE', 'AluY', 3, 4),
        ('locus3', 'L1HS_locus3', 'L1',  'LINE', 'L1HS', 7, 0),
    ])
    out = hec.parse_squire_count(str(p), ['s1', 's2'])
    assert out['AluY'] == {'s1': 8, 's2': 6}
    assert out['L1HS'] == {'s1': 7}


def test_check_overlap_passes_above_threshold(capsys):
    hec.check_overlap(['AluY', 'L1HS'], ['AluY', 'L1HS', 'MIR'],
                      threshold=0.5, label='gene_id')
    captured = capsys.readouterr()
    assert '2/2' in captured.err


def test_check_overlap_fails_below_threshold():
    with pytest.raises(AssertionError):
        hec.check_overlap(['rogue1', 'rogue2'], ['AluY', 'L1HS'],
                          threshold=0.5, label='gene_id')


def test_write_count_matrix_shape_and_ordering(tmp_path):
    counts = {'L1HS': {'s1': 5}, 'AluY': {'s1': 10, 's2': 20}}
    out_path = tmp_path / 'out.tsv'
    hec.write_count_matrix(counts, ['s1', 's2'], str(out_path))
    with open(out_path) as fh:
        lines = [line.rstrip('\n').split('\t') for line in fh]
    assert lines[0] == ['feature_id', 's1', 's2']
    feat_rows = [r[0] for r in lines[1:]]
    assert feat_rows == sorted(feat_rows)
    aluy_row = next(r for r in lines[1:] if r[0] == 'AluY')
    assert aluy_row == ['AluY', '10', '20']


def test_write_count_matrix_empty_raises(tmp_path):
    with pytest.raises(AssertionError):
        hec.write_count_matrix({}, ['s1'], str(tmp_path / 'x.tsv'))


def test_cli_end_to_end_tetranscripts(tmp_path):
    """End-to-end CLI invocation: TEtranscripts table + locus_map ->
    harmonised TSV at family_id granularity."""
    lm = tmp_path / 'lm.tsv'
    write_tsv(lm, [
        ('transcript_id', 'gene_id', 'family_id', 'class_id'),
        ('AluY_chr1_100', 'AluY', 'Alu', 'SINE'),
        ('L1HS_chr3_300', 'L1HS', 'L1', 'LINE'),
    ])
    te = tmp_path / 'te.tsv'
    write_tsv(te, [
        ('feature_id', 's1', 's2'),
        ('AluY:Alu:SINE', 12, 8),
        ('L1HS:L1:LINE', 0, 4),
    ])
    out = tmp_path / 'harmonised.tsv'
    script = REPO_ROOT / 'workflow' / 'scripts' / 'harmonize_external_counts.py'
    res = subprocess.run([
        sys.executable, str(script),
        '--tool', 'tetranscripts',
        '--input', str(te),
        '--locus-map', str(lm),
        '--samples', 's1,s2',
        '--granularity', 'family_id',
        '--min-overlap-fraction', '0.0',
        '--output', str(out),
    ], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    with open(out) as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh]
    assert rows[0] == ['feature_id', 's1', 's2']
    fam_rows = {r[0]: r[1:] for r in rows[1:]}
    assert fam_rows['Alu'] == ['12', '8']
    assert fam_rows['L1'] == ['0', '4']


def test_parse_barcode_map_roundtrip(tmp_path):
    p = tmp_path / 'bc.tsv'
    write_tsv(p, [
        ('barcode', 'cell_id'),
        ('ACGTACGTACGTACGT', 'cell_001'),
        ('TTTTGGGGAAAACCCC', 'cell_002'),
    ])
    out = hec.parse_barcode_map(str(p))
    assert out == {
        'ACGTACGTACGTACGT': 'cell_001',
        'TTTTGGGGAAAACCCC': 'cell_002',
    }


def test_parse_barcode_map_skips_header_only(tmp_path):
    p = tmp_path / 'empty_bc.tsv'
    write_tsv(p, [('barcode', 'cell_id')])
    assert hec.parse_barcode_map(str(p)) == {}
