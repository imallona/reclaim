"""
Unit tests for workflow/scripts/sc_count_features.py.

This script is the alternative repeat-counting path: a side counter that
re-counts a STARsolo-tagged BAM against a feature GTF using the cell
barcode whitelist from STARsolo's gene-counting pass. It will eventually
replace the per-feature_set STARsolo runs, so its correctness is
load-bearing for the paper. The tests below exercise:

  GTF parsing
    - gene_id, family_id, class_id extraction
    - level selection: gene > transcript > exon (intron-inclusive)
    - exon-only repeat GTFs fall back cleanly

  Read assignment
    - exonic read counts toward the gene
    - intronic read counts when GTF has gene-level entries (GeneFull-style)
    - intronic read does NOT count when GTF has only exon entries
    - read overlapping two distinct gene_ids → discarded (ambiguous)
    - read off any feature → discarded
    - unmapped/secondary/supplementary records → ignored

  Cell barcode and UMI logic
    - missing CB tag → skipped, no contribution
    - missing UB tag → skipped, no contribution
    - off-whitelist CB → skipped
    - duplicate UMIs for the same (CB, gene_id) collapse to one count
    - same UMI appearing for two different (CB, gene_id) pairs counts in
      both (they are independent dedup spaces)

  Multimapper handling
    - --multimapper unique: NH > 1 reads dropped
    - --multimapper multi : NH > 1 reads kept (counted with their UB)

  Granularity rollup
    - family_id and class_id sums match per-gene counts
    - rollup is post-UMI-dedup (sum of len(set) values)

  Output integrity
    - Matrix Market header is well-formed
    - non-zero entries match the expected dictionary
    - barcodes.tsv preserves whitelist order
    - features.tsv lists only features with at least one count
"""

import os
import sys

import pytest

pysam = pytest.importorskip('pysam')
pytest.importorskip('intervaltree')

import sc_count_features as scf  # noqa: E402  (after path setup in conftest)


# Helper to write a small GTF with gene/exon entries.
def _write_gtf(path, records):
    """Write GTF lines from a list of dicts.

    Each record: {chrom, type, start, end, strand, gene_id, family_id?, class_id?}
    Coordinates are 1-based inclusive (GTF convention).
    """
    with open(path, 'w') as fh:
        for r in records:
            attrs = [f'gene_id "{r["gene_id"]}"',
                     f'transcript_id "{r["gene_id"]}_t1"']
            if 'family_id' in r:
                attrs.append(f'family_id "{r["family_id"]}"')
            if 'class_id' in r:
                attrs.append(f'class_id "{r["class_id"]}"')
            fh.write('\t'.join([
                r['chrom'], 'test', r['type'],
                str(r['start']), str(r['end']),
                '.', r.get('strand', '+'), '.',
                '; '.join(attrs) + ';',
            ]) + '\n')


# Helper to write a BAM with CB/UB-tagged reads on a single chromosome.
def _write_bam(path, references, reads):
    """Write a sorted+indexed BAM.

    references: list of (name, length) tuples
    reads: list of dicts {chrom, pos (0-based), seq, cb?, ub?, nh?, flag?}
    """
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'LN': L, 'SN': N} for N, L in references]}
    chrom_index = {N: i for i, (N, _) in enumerate(references)}
    with pysam.AlignmentFile(path, 'wb', header=header) as bf:
        # Snakemake-style fast path: build records in pos order per chrom.
        sorted_reads = sorted(reads,
                              key=lambda r: (chrom_index[r['chrom']], r['pos']))
        for i, r in enumerate(sorted_reads):
            a = pysam.AlignedSegment()
            a.query_name = f'r{i}'
            a.query_sequence = r['seq']
            a.flag = r.get('flag', 0)
            a.reference_id = chrom_index[r['chrom']]
            a.reference_start = r['pos']
            a.cigar = [(0, len(r['seq']))]  # M only
            a.mapping_quality = 60
            tags = []
            if 'cb' in r:
                tags.append(('CB', r['cb']))
            if 'ub' in r:
                tags.append(('UB', r['ub']))
            tags.append(('NH', r.get('nh', 1)))
            a.tags = tags
            bf.write(a)
    pysam.sort('-o', path + '.sorted.bam', path)
    os.replace(path + '.sorted.bam', path)
    pysam.index(path)


@pytest.fixture
def whitelist_path(tmp_path):
    p = tmp_path / 'barcodes.tsv'
    p.write_text('AAAA\nBBBB\nCCCC\n')
    return str(p)


@pytest.fixture
def repeat_gtf(tmp_path):
    p = tmp_path / 'repeats.gtf'
    _write_gtf(str(p), [
        # Two non-overlapping repeats on chr1.
        {'chrom': 'chr1', 'type': 'exon', 'start': 101, 'end': 200,
         'gene_id': 'L1HS_dup1', 'family_id': 'L1', 'class_id': 'LINE'},
        {'chrom': 'chr1', 'type': 'exon', 'start': 301, 'end': 400,
         'gene_id': 'AluY_dup1', 'family_id': 'Alu', 'class_id': 'SINE'},
        # Repeat on chr2.
        {'chrom': 'chr2', 'type': 'exon', 'start': 501, 'end': 600,
         'gene_id': 'MIR_dup1', 'family_id': 'MIR', 'class_id': 'SINE'},
    ])
    return str(p)


@pytest.fixture
def gene_gtf(tmp_path):
    """A gene-level GTF with introns: gene[100-1000], exons[100-200, 800-1000].
    The intron 200-800 is covered by the gene record but not by the exons."""
    p = tmp_path / 'genes.gtf'
    _write_gtf(str(p), [
        {'chrom': 'chr1', 'type': 'gene', 'start': 101, 'end': 1000,
         'gene_id': 'GeneA'},
        {'chrom': 'chr1', 'type': 'exon', 'start': 101, 'end': 200,
         'gene_id': 'GeneA'},
        {'chrom': 'chr1', 'type': 'exon', 'start': 801, 'end': 1000,
         'gene_id': 'GeneA'},
        # Second gene further along.
        {'chrom': 'chr1', 'type': 'gene', 'start': 1501, 'end': 2000,
         'gene_id': 'GeneB'},
        {'chrom': 'chr1', 'type': 'exon', 'start': 1501, 'end': 2000,
         'gene_id': 'GeneB'},
    ])
    return str(p)


@pytest.fixture
def overlap_gtf(tmp_path):
    """Two genes whose exon ranges overlap on chr1 — a read landing on the
    overlap should be ambiguous."""
    p = tmp_path / 'overlap.gtf'
    _write_gtf(str(p), [
        {'chrom': 'chr1', 'type': 'exon', 'start': 101, 'end': 200,
         'gene_id': 'GeneL'},
        {'chrom': 'chr1', 'type': 'exon', 'start': 151, 'end': 250,
         'gene_id': 'GeneR'},
    ])
    return str(p)


def _read_mtx(path):
    """Parse a Matrix Market file into a dict[(row, col)] = value plus the
    header dimensions."""
    with open(path) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    if not lines[0].startswith('%%MatrixMarket'):
        raise AssertionError('missing MatrixMarket header')
    body = [ln for ln in lines[1:] if not ln.startswith('%')]
    n_rows, n_cols, n_entries = (int(x) for x in body[0].split())
    entries = {}
    for ln in body[1:]:
        r, c, v = ln.split()
        entries[(int(r), int(c))] = int(v)
    assert len(entries) == n_entries
    return n_rows, n_cols, entries


def _read_features(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            out.append(ln.rstrip('\n').split('\t')[0])
    return out


# GTF parsing


def test_gtf_extracts_gene_family_class(repeat_gtf):
    intervals, fam, cls, level = scf.parse_gtf(repeat_gtf)
    assert level == 'exon'
    assert sorted(intervals) == ['chr1', 'chr2']
    assert fam == {'L1HS_dup1': 'L1', 'AluY_dup1': 'Alu', 'MIR_dup1': 'MIR'}
    assert cls == {'L1HS_dup1': 'LINE', 'AluY_dup1': 'SINE', 'MIR_dup1': 'SINE'}


def test_gtf_prefers_gene_over_exon(gene_gtf):
    intervals, fam, cls, level = scf.parse_gtf(gene_gtf)
    assert level == 'gene'
    assert ('chr1' in intervals) and len(intervals['chr1']) == 2
    # Indexed gene span includes the intron.
    spans = sorted(intervals['chr1'])
    assert spans[0] == (100, 1000, 'GeneA')
    assert fam == {} and cls == {}


def test_gtf_falls_back_to_exon_when_no_gene(repeat_gtf):
    _, _, _, level = scf.parse_gtf(repeat_gtf)
    assert level == 'exon'


# Read assignment behaviour (using the in-process worker init path).


def _run(bam_path, gtf_path, whitelist_path, multimapper='unique', threads=1,
         umi_dedup='exact', cb_tag='CB'):
    """Drive the script in-process and return:

      cbs                - whitelist barcodes in order
      raw_counts         - {(cb_idx, gene_id): {umi: weight}} (pre-dedup)
      collapsed_counts   - {(cb_idx, gene_id): scalar} after the chosen
                            UMI-dedup method (default exact, returning
                            the count of distinct UMIs so tests can
                            compare against integer expectations)
      fam, cls           - gene_id -> family_id and -> class_id maps from
                            the GTF
    """
    cbs = scf.load_whitelist(whitelist_path)
    cb_index = {cb: i for i, cb in enumerate(cbs)}
    intervals, fam, cls, _ = scf.parse_gtf(gtf_path)

    bam_header = pysam.AlignmentFile(bam_path, 'rb')
    contigs = list(bam_header.references)
    bam_header.close()

    scf._worker_init(bam_path, intervals, cb_index, multimapper, umi_dedup, cb_tag)
    results = [scf._process_contig(c) for c in contigs]
    per_contig = [p for _, p, _ in results]
    raw_counts = scf.merge_counts(per_contig)
    collapsed_counts = scf.collapse(raw_counts, umi_dedup)
    return cbs, raw_counts, collapsed_counts, fam, cls


def test_exonic_read_counts(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 120, 'seq': 'A' * 30, 'cb': 'AAAA', 'ub': 'U1'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    assert counts == {(0, 'L1HS_dup1'): 1}


def test_intronic_read_counts_with_gene_level(tmp_path, gene_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 3000)], [
        # Read sits in the intron 201-800 — outside any exon, inside the gene.
        {'chrom': 'chr1', 'pos': 400, 'seq': 'A' * 30, 'cb': 'AAAA', 'ub': 'U1'},
    ])
    _, _, counts, _, _ = _run(bam, gene_gtf, whitelist_path)
    assert counts == {(0, 'GeneA'): 1}


def test_intronic_read_skipped_with_exon_only_gtf(
        tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        # Between the two repeats — not inside any exon.
        {'chrom': 'chr1', 'pos': 250, 'seq': 'A' * 20, 'cb': 'AAAA', 'ub': 'U1'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    assert counts == {}


def test_ambiguous_overlap_discarded(tmp_path, overlap_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000)], [
        # Read at 170 spans both GeneL[100-200] and GeneR[150-250].
        {'chrom': 'chr1', 'pos': 165, 'seq': 'A' * 20, 'cb': 'AAAA', 'ub': 'U1'},
    ])
    _, _, counts, _, _ = _run(bam, overlap_gtf, whitelist_path)
    assert counts == {}


def test_unmapped_and_secondary_skipped(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'flag': 4,  # unmapped
         'cb': 'AAAA', 'ub': 'U1'},
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'flag': 256,  # secondary
         'cb': 'AAAA', 'ub': 'U2'},
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'flag': 2048,  # supplementary
         'cb': 'AAAA', 'ub': 'U3'},
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,  # primary, kept
         'cb': 'AAAA', 'ub': 'U4'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    assert counts == {(0, 'L1HS_dup1'): 1}


# CB / UB / whitelist behaviour.


def test_missing_cb_or_ub_skipped(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'ub': 'U1'},   # no CB
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'cb': 'AAAA'},  # no UB
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U_kept'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    assert counts == {(0, 'L1HS_dup1'): 1}


def test_off_whitelist_cb_skipped(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'ZZZZ', 'ub': 'U1'},   # not in whitelist
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'BBBB', 'ub': 'U2'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    assert counts == {(1, 'L1HS_dup1'): 1}


def test_umi_dedup_within_cb_gene(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SAME'},
        {'chrom': 'chr1', 'pos': 115, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SAME'},
        {'chrom': 'chr1', 'pos': 120, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'OTHER'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    # Two distinct UMIs counted once each.
    assert counts == {(0, 'L1HS_dup1'): 2}


def test_umi_independent_across_genes(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SHARED'},   # → L1HS
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SHARED'},   # → AluY (different gene)
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path)
    assert counts == {(0, 'L1HS_dup1'): 1,
                      (0, 'AluY_dup1'): 1}


# Multimapper handling.


def test_multimapper_unique_skips_nh_gt_1(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U1', 'nh': 3},
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U2', 'nh': 1},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              multimapper='unique')
    assert counts == {(0, 'L1HS_dup1'): 1}


def test_multimapper_multi_keeps_nh_gt_1(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U1', 'nh': 3},
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U2', 'nh': 1},
    ])
    # With --multimapper multi, NH=3 contributes weight 1/3 and NH=1
    # contributes weight 1.0; sum (under exact dedup) = 4/3.
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              multimapper='multi')
    assert (0, 'L1HS_dup1') in counts
    assert counts[(0, 'L1HS_dup1')] == pytest.approx(1 / 3 + 1.0)


def test_multimapper_multi_distributes_secondary_records(
        tmp_path, repeat_gtf, whitelist_path):
    """STAR emits multimappers as one primary + (NH-1) secondary records,
    all sharing UB. multi mode must include the secondary records, else
    the molecule's contribution collapses to 1/NH at the primary's gene
    instead of summing to ~1.0 across the alternative loci."""
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        # primary alignment on L1HS_dup1
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SHARED', 'nh': 2, 'flag': 0},
        # secondary alignment of the same molecule on AluY_dup1
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SHARED', 'nh': 2, 'flag': 256},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              multimapper='multi')
    # Both genes should get weight 1/NH = 0.5; sum across them is 1.0.
    assert counts[(0, 'L1HS_dup1')] == pytest.approx(0.5)
    assert counts[(0, 'AluY_dup1')] == pytest.approx(0.5)


def test_multimapper_unique_skips_secondary_records(
        tmp_path, repeat_gtf, whitelist_path):
    """unique mode must drop secondary records; only the primary
    contributes, and its NH > 1 means it's also dropped, so the molecule
    is excluded entirely (which is the documented unique-only semantic)."""
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SHARED', 'nh': 2, 'flag': 0},
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'SHARED', 'nh': 2, 'flag': 256},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              multimapper='unique')
    assert counts == {}


# Granularity rollup.


def test_family_class_rollup_post_dedup(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        # Same CB, two distinct UMIs on L1HS_dup1 (family L1, class LINE).
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U1'},
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U2'},
        # One UMI on AluY_dup1 (family Alu, class SINE).
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U3'},
        # One UMI on MIR_dup1 (family MIR, class SINE).
        {'chrom': 'chr2', 'pos': 520, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U4'},
    ])
    _, _, counts, fam, cls = _run(bam, repeat_gtf, whitelist_path)
    fam_counts = scf.rollup(counts, fam)
    cls_counts = scf.rollup(counts, cls)
    assert fam_counts[(0, 'L1')] == 2
    assert fam_counts[(0, 'Alu')] == 1
    assert fam_counts[(0, 'MIR')] == 1
    # SINE = AluY + MIR = 2; LINE = L1HS = 2.
    assert cls_counts[(0, 'SINE')] == 2
    assert cls_counts[(0, 'LINE')] == 2


# End-to-end: matrix file integrity.


def test_end_to_end_outputs(tmp_path, repeat_gtf, whitelist_path, monkeypatch):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'U1'},
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20,
         'cb': 'BBBB', 'ub': 'U1'},
    ])
    out = tmp_path / 'out'
    monkeypatch.setattr(sys, 'argv', [
        'sc_count_features',
        '--bam', bam,
        '--gtf', repeat_gtf,
        '--whitelist', whitelist_path,
        '--out-dir', str(out),
        '--multimapper', 'unique',
        '--threads', '1',
    ])
    scf.main()

    gene_dir = out / 'Solo.out' / 'Gene' / 'raw'
    assert gene_dir.is_dir()
    rows, cols, entries = _read_mtx(str(gene_dir / 'matrix.mtx'))
    assert cols == 3  # whitelist size
    feats = _read_features(str(gene_dir / 'features.tsv'))
    assert sorted(feats) == ['AluY_dup1', 'L1HS_dup1']
    barcodes = (gene_dir / 'barcodes.tsv').read_text().splitlines()
    assert barcodes == ['AAAA', 'BBBB', 'CCCC']
    # Two non-zero entries: (L1HS_dup1, AAAA, 1) and (AluY_dup1, BBBB, 1)
    # (feature index by sorted gene_id).
    feat_index = {f: i + 1 for i, f in enumerate(sorted(feats))}
    assert entries[(feat_index['L1HS_dup1'], 1)] == 1
    assert entries[(feat_index['AluY_dup1'], 2)] == 1

    fam_dir = out / 'Solo.out' / 'Gene_family_id' / 'raw'
    assert fam_dir.is_dir()
    fam_rows, fam_cols, fam_entries = _read_mtx(str(fam_dir / 'matrix.mtx'))
    fam_feats = _read_features(str(fam_dir / 'features.tsv'))
    fam_index = {f: i + 1 for i, f in enumerate(sorted(fam_feats))}
    assert fam_entries[(fam_index['L1'], 1)] == 1
    assert fam_entries[(fam_index['Alu'], 2)] == 1

    cls_dir = out / 'Solo.out' / 'Gene_class_id' / 'raw'
    assert cls_dir.is_dir()
    _, _, cls_entries = _read_mtx(str(cls_dir / 'matrix.mtx'))
    cls_feats = _read_features(str(cls_dir / 'features.tsv'))
    cls_index = {f: i + 1 for i, f in enumerate(sorted(cls_feats))}
    # LINE for AAAA, SINE for BBBB.
    assert cls_entries[(cls_index['LINE'], 1)] == 1
    assert cls_entries[(cls_index['SINE'], 2)] == 1


# --umi-dedup none semantics: SmartSeq2 path with no real UMIs.
# Reads are keyed by a synthetic (qname, ref_start), so every record
# contributes its weight regardless of whether a UB tag is present or
# whether multiple reads share the same UB.


def test_dedup_none_counts_every_read(tmp_path, repeat_gtf, whitelist_path):
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        # Three reads at the L1HS locus, all with the SAME fixed UB. Under
        # exact/1mm_all this would collapse to 1; under none it stays 3.
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'GATTACAGCTA'},
        {'chrom': 'chr1', 'pos': 130, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'GATTACAGCTA'},
        {'chrom': 'chr1', 'pos': 150, 'seq': 'A' * 20,
         'cb': 'AAAA', 'ub': 'GATTACAGCTA'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              umi_dedup='none')
    assert counts == {(0, 'L1HS_dup1'): 3}


def test_dedup_none_works_without_ub_tag(tmp_path, repeat_gtf, whitelist_path):
    """In none mode, missing UB tags must NOT drop the read; the synthetic
    key is built from the alignment record itself, not from UB."""
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'cb': 'AAAA'},
        {'chrom': 'chr1', 'pos': 130, 'seq': 'A' * 20, 'cb': 'AAAA'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              umi_dedup='none')
    assert counts == {(0, 'L1HS_dup1'): 2}


def test_dedup_none_multimapper_distributes_across_loci(
        tmp_path, repeat_gtf, whitelist_path):
    """A multi-mapper (NH=2) primary + secondary should still distribute
    weight 1/NH across loci under --multimapper multi + --umi-dedup none.
    The synthetic key (qname, ref_start) keeps primary and secondary
    distinct because they sit at different ref_start."""
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20,
         'cb': 'AAAA', 'nh': 2, 'flag': 0},
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20,
         'cb': 'AAAA', 'nh': 2, 'flag': 256},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              multimapper='multi', umi_dedup='none')
    assert counts[(0, 'L1HS_dup1')] == pytest.approx(0.5)
    assert counts[(0, 'AluY_dup1')] == pytest.approx(0.5)


def test_dedup_none_off_whitelist_still_skipped(
        tmp_path, repeat_gtf, whitelist_path):
    """none mode does not bypass whitelist filtering: a read whose CB is
    not in the whitelist still contributes nothing."""
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'cb': 'ZZZZ'},
        {'chrom': 'chr1', 'pos': 130, 'seq': 'A' * 20, 'cb': 'AAAA'},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              umi_dedup='none')
    assert counts == {(0, 'L1HS_dup1'): 1}


def test_dedup_none_endtoend_outputs(tmp_path, repeat_gtf, whitelist_path,
                                      monkeypatch):
    """End-to-end: run main() with --umi-dedup none on a SS2-like BAM that
    has CB but no UB, and verify the matrix sums match read counts."""
    bam = str(tmp_path / 'r.bam')
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        # Cell AAAA: 2 reads on L1HS, 1 on AluY.
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'cb': 'AAAA'},
        {'chrom': 'chr1', 'pos': 130, 'seq': 'A' * 20, 'cb': 'AAAA'},
        {'chrom': 'chr1', 'pos': 320, 'seq': 'A' * 20, 'cb': 'AAAA'},
        # Cell BBBB: 1 read on MIR.
        {'chrom': 'chr2', 'pos': 520, 'seq': 'A' * 20, 'cb': 'BBBB'},
    ])
    out = tmp_path / 'out'
    monkeypatch.setattr(sys, 'argv', [
        'sc_count_features',
        '--bam', bam,
        '--gtf', repeat_gtf,
        '--whitelist', whitelist_path,
        '--out-dir', str(out),
        '--multimapper', 'unique',
        '--umi-dedup', 'none',
        '--threads', '1',
    ])
    scf.main()

    gene_dir = out / 'Solo.out' / 'Gene' / 'raw'
    rows, cols, entries = _read_mtx(str(gene_dir / 'matrix.mtx'))
    feats = _read_features(str(gene_dir / 'features.tsv'))
    feat_index = {f: i + 1 for i, f in enumerate(sorted(feats))}
    # AAAA = column 1, BBBB = column 2 (whitelist order).
    assert entries[(feat_index['L1HS_dup1'], 1)] == 2
    assert entries[(feat_index['AluY_dup1'], 1)] == 1
    assert entries[(feat_index['MIR_dup1'],  2)] == 1


def test_paired_bam_skips_read2_to_avoid_fragment_double_count(
        tmp_path, repeat_gtf, whitelist_path):
    """Paired SmartSeq2 BAMs deposit R1 and R2 of the same fragment at
    different ref_start. Without an is_read2 filter the synthetic UB key
    keeps them distinct and the fragment is counted twice. Skipping R2
    when is_paired is set collapses the pair to one count."""
    bam = str(tmp_path / 'r.bam')
    # Two fragments on L1HS, each as a R1 + R2 pair. Both reads of each
    # fragment carry is_paired (flag bit 1) and is_read1 (flag bit 64) /
    # is_read2 (flag bit 128). We only need the flag bits set; the BAM
    # writer doesn't validate mate fields.
    F1, F2 = 0x1 | 0x40, 0x1 | 0x80   # paired + read1, paired + read2
    _write_bam(bam, [('chr1', 1000), ('chr2', 1000)], [
        # Fragment 1, R1 and R2.
        {'chrom': 'chr1', 'pos': 110, 'seq': 'A' * 20, 'cb': 'AAAA',
         'flag': F1},
        {'chrom': 'chr1', 'pos': 130, 'seq': 'A' * 20, 'cb': 'AAAA',
         'flag': F2},
        # Fragment 2, R1 and R2.
        {'chrom': 'chr1', 'pos': 150, 'seq': 'A' * 20, 'cb': 'AAAA',
         'flag': F1},
        {'chrom': 'chr1', 'pos': 170, 'seq': 'A' * 20, 'cb': 'AAAA',
         'flag': F2},
    ])
    _, _, counts, _, _ = _run(bam, repeat_gtf, whitelist_path,
                              umi_dedup='none')
    # Two fragments, R2 skipped, so 2 counts not 4.
    assert counts == {(0, 'L1HS_dup1'): 2}
