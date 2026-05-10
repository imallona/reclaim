"""
Unit tests for workflow/scripts/parse_star_log.py.

The parser scans STAR Log.final.out files and emits a tidy TSV with mapping
rate, mismatch rate, deletion/insertion rates, and read counts. The tests
here cover:

  - key:value extraction and column mapping
  - percent-suffix stripping
  - sample_id / multimapper_mode / feature_set inference from the log path
  - missing keys produce empty cells (no crash, no garbage)
  - end-to-end main() emits a TSV with the correct schema
"""

import csv
import os
import sys

import pytest

import parse_star_log as psl


SAMPLE_LOG = """\
                                 Started job on |\tSat May 02 06:22:00 2026
                             Started mapping on |\tSat May 02 06:22:00 2026
                                    Finished on |\tSat May 02 07:18:42 2026
       Mapping speed, Million of reads per hour |\t12.34

                          Number of input reads |\t11719321
                      Average input read length |\t90
                                    UNIQUE READS:
                   Uniquely mapped reads number |\t9384732
                        Uniquely mapped reads % |\t80.07%
                          Average mapped length |\t89.21
                       Number of splices: Total |\t214120
             Number of splices: Annotated (sjdb) |\t214120
                       Number of splices: GT/AG |\t211000
                       Number of splices: GC/AG |\t2900
                       Number of splices: AT/AC |\t150
                  Number of splices: Non-canonical |\t70
                          Mismatch rate per base, % |\t0.32%
                                Deletion rate per base |\t0.01%
                              Deletion average length |\t1.31
                              Insertion rate per base |\t0.01%
                              Insertion average length |\t1.50

                              MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |\t1320132
             % of reads mapped to multiple loci |\t11.27%
        Number of reads mapped to too many loci |\t12345
             % of reads mapped to too many loci |\t0.11%

                                  UNMAPPED READS:
       Number of reads unmapped: too many mismatches |\t0
            % of reads unmapped: too many mismatches |\t0.00%
            Number of reads unmapped: too short |\t1002112
                 % of reads unmapped: too short |\t8.55%
                Number of reads unmapped: other |\t0
                     % of reads unmapped: other |\t0.00%

                                  CHIMERIC READS:
                       Number of chimeric reads |\t0
                            % of chimeric reads |\t0.00%
"""


@pytest.fixture
def log_at_path(tmp_path):
    """Drop a synthetic Log.final.out into the conventional layout
       <starsolo_base>/<sample>/<mode>_<feature>/Log.final.out so the parser
       can recover sample_id, multimapper_mode, feature_set from the path."""
    leaf = tmp_path / 'starsolo' / 'KD6_15_essential' / 'unique_genes'
    leaf.mkdir(parents=True)
    log_path = leaf / 'Log.final.out'
    log_path.write_text(SAMPLE_LOG)
    return str(log_path)


def test_parse_one_extracts_all_known_keys(log_at_path):
    row = psl.parse_one(log_at_path)
    assert row['n_input_reads'] == '11719321'
    assert row['avg_input_length'] == '90'
    assert row['avg_mapped_length'] == '89.21'
    assert row['uniquely_mapped_pct'] == '80.07'
    assert row['multi_mapped_pct'] == '11.27'
    assert row['too_many_loci_pct'] == '0.11'
    assert row['unmapped_too_short_pct'] == '8.55'
    assert row['unmapped_other_pct'] == '0.00'
    assert row['mismatch_rate_pct'] == '0.32'
    assert row['deletion_rate_pct'] == '0.01'
    assert row['insertion_rate_pct'] == '0.01'
    assert row['deletion_avg_length'] == '1.31'
    assert row['insertion_avg_length'] == '1.50'


def test_parse_one_strips_pct_suffix(log_at_path):
    row = psl.parse_one(log_at_path)
    for k in ('uniquely_mapped_pct', 'multi_mapped_pct', 'mismatch_rate_pct'):
        assert '%' not in row[k]


def test_derive_path_fields(log_at_path):
    sample_id, mode, fs = psl.derive_path_fields(log_at_path)
    assert sample_id == 'KD6_15_essential'
    assert mode == 'unique'
    assert fs == 'genes'


def test_derive_path_fields_unparseable_returns_empties(tmp_path):
    p = tmp_path / 'unrelated' / 'flat' / 'Log.final.out'
    p.parent.mkdir(parents=True)
    p.write_text('')
    sample_id, mode, fs = psl.derive_path_fields(str(p))
    # grandparent name maps to sample_id even when 'flat' has no underscore;
    # mode/fs fall back to empty strings.
    assert sample_id == 'unrelated'
    assert mode == ''
    assert fs == ''


def test_parse_missing_keys_become_empty(tmp_path):
    """If the log lacks some known keys, those columns should be empty
       rather than raising or substituting placeholder text."""
    p = tmp_path / 'minimal' / 'multi_genic_repeats' / 'Log.final.out'
    p.parent.mkdir(parents=True)
    p.write_text(
        '                          Number of input reads |\t100\n'
        '                          Mismatch rate per base, % |\t0.50%\n'
    )
    row = psl.parse_one(str(p))
    assert row['n_input_reads'] == '100'
    assert row['mismatch_rate_pct'] == '0.50'
    # Keys absent from the log are blank, not missing from the dict.
    assert row['uniquely_mapped_pct'] == ''
    assert row['avg_mapped_length'] == ''


def test_main_writes_well_formed_tsv(log_at_path, tmp_path, monkeypatch):
    out_path = tmp_path / 'out.tsv'
    monkeypatch.setattr(sys, 'argv', [
        'parse_star_log',
        '--log', log_at_path,
        '--out', str(out_path),
    ])
    psl.main()

    assert out_path.exists()
    with open(out_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = list(reader)
    assert len(rows) == 1
    row = rows[0]
    expected_columns = {
        'sample_id', 'multimapper_mode', 'feature_set', 'log_path',
        'n_input_reads', 'avg_input_length', 'avg_mapped_length',
        'uniquely_mapped_pct', 'multi_mapped_pct', 'too_many_loci_pct',
        'unmapped_too_short_pct', 'unmapped_other_pct',
        'mismatch_rate_pct', 'deletion_rate_pct', 'insertion_rate_pct',
        'deletion_avg_length', 'insertion_avg_length',
    }
    assert set(row.keys()) == expected_columns
    assert row['sample_id'] == 'KD6_15_essential'
    assert row['multimapper_mode'] == 'unique'
    assert row['feature_set'] == 'genes'
    assert row['mismatch_rate_pct'] == '0.32'


def test_main_handles_multiple_logs(tmp_path, log_at_path, monkeypatch):
    """When given multiple --log args, each becomes one row."""
    second = tmp_path / 'starsolo' / 'KD6_15_essential' / 'multi_intergenic_repeats'
    second.mkdir(parents=True)
    (second / 'Log.final.out').write_text(SAMPLE_LOG)
    out_path = tmp_path / 'out.tsv'
    monkeypatch.setattr(sys, 'argv', [
        'parse_star_log',
        '--log', log_at_path,
        '--log', str(second / 'Log.final.out'),
        '--out', str(out_path),
    ])
    psl.main()
    with open(out_path) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    assert len(rows) == 2
    modes = sorted(r['multimapper_mode'] for r in rows)
    assert modes == ['multi', 'unique']
    fs = sorted(r['feature_set'] for r in rows)
    assert fs == ['genes', 'intergenic_repeats']
