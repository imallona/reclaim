#!/usr/bin/env python3
"""
Single-cell feature counter for STARsolo-tagged BAMs.

Counts reads from a STARsolo BAM against a feature GTF, with cell-barcode
whitelist filtering and UMI deduplication. Output mirrors STARsolo's
Solo.out/Gene/raw/ layout (matrix.mtx + features.tsv + barcodes.tsv).

Why this script:
  - STARsolo has no native re-count-from-BAM mode. To count a different
    feature_set against a gene-aligned BAM, you either re-align (expensive)
    or use a side counter.
  - Cell whitelist matters: STARsolo's "raw" matrix has one column per CB
    seen during alignment (millions, mostly noise). Restricting to the
    genes-pass barcodes.tsv keeps the matrix small and meaningful.
  - featureCounts is the obvious alternative but loads alignments into
    memory; this script streams the BAM and keeps memory bounded by the
    UMI dedup state (typically <1 GB per library).

Indexing strategy (intron-inclusive, like STARsolo --soloFeatures GeneFull):
  - If the GTF has `gene` records, those are indexed (one interval per
    gene spanning the entire transcript including introns; reads in
    introns count toward the gene).
  - Else if `transcript` records exist, those are indexed.
  - Else `exon` records (the repeat GTFs use exon-only with one exon per
    repeat locus; the entire repeat span is then covered by definition).
  This is more inclusive than STARsolo's default --soloFeatures Gene
  (which only counts exonic reads) and matches GeneFull semantics.

Assignment policy:
  - A read counts toward gene_id G iff its alignment overlaps exactly one
    GTF gene's indexed interval. Multi-feature overlap → discarded.
  - Reads without CB or UB tags are skipped.
  - Reads whose CB is not in the whitelist are skipped.

Multimapper handling (--multimapper):
  - unique : skip reads with NH > 1. Matches STARsolo's `Unique` setting
             exactly; pysam-unique and STARsolo-Unique counts on the same
             BAM and same dedup setting must agree.
  - multi  : include NH > 1 reads with weight 1/NH per alignment record.
             A multimapper read with NH=3 contributes 1/3 to each of its
             three feature assignments. This is non-iterative and
             linear-cost, unlike STARsolo's --soloMultiMappers EM which
             iterates a redistribution step. Numerical differences from
             EM are systematic: EM concentrates weight on genes with
             stronger unique-mapper support, while 1/NH spreads it
             uniformly across the candidate loci. For repeat-family
             quantification, where multimappers usually come from
             uniformly-mappable young subfamilies, the 1/NH uniform spread
             is the more conservative and arguably more honest choice.

             Known limitation: when a multimapper has more than one
             alternative falling in the same gene (e.g., two paralog
             positions of the same family), the implementation stores
             the max weight per (CB, gene, UMI) instead of summing
             positions. The molecule then contributes only 1/NH to that
             gene rather than k/NH (k = # alternatives in that gene).
             Sum across genes is therefore <= 1 for multimappers with
             intra-gene collisions, vs STARsolo's EM which would credit
             the gene proportionally to its alternative density. The
             simulation benchmark quantifies this bias against ground
             truth.

UMI deduplication (--umi-dedup):
  - 1mm_all    (default): collapse UMIs within Hamming distance 1 of each
                          other via union-find on the per-(CB, gene_id)
                          UMI graph. Matches STARsolo's default
                          --soloUMIdedup 1MM_All.
  - exact               : count distinct UB strings. Faster but
                          systematically over-counts unique molecules
                          because sequencing errors create spurious UMIs.
                          Use only for testing/debugging.

Granularity outputs:
  When the GTF carries family_id and/or class_id attributes (the repeat
  GTFs do; Ensembl genes do not), the script also writes per-granularity
  matrices side-by-side. Family/class counts are the post-UMI-dedup sum
  of their member gene_ids — equivalent to the rollup performed by
  granularity_rollup.R but produced in one pass.

Parallelism:
  - One worker per BAM contig. Each worker opens its own pysam handle
    (pysam handles are not fork-safe), processes one contig, and returns
    its local dict {(cb_idx, gene_id) -> set(UMIs)}. The main process
    merges per-key sets.

Output layout:
  <out_dir>/Solo.out/Gene/raw/{matrix.mtx|UniqueAndMult-EM.mtx,
                               features.tsv, barcodes.tsv}      gene_id
  <out_dir>/Solo.out/Gene_family_id/raw/...                     family_id
  <out_dir>/Solo.out/Gene_class_id/raw/...                      class_id
"""

import argparse
import gzip
import multiprocessing as mp
import os
import sys
from collections import defaultdict


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')


def parse_attribute(attrs, key):
    """Pull the first matching key from a GTF attributes column."""
    needle = key + ' "'
    i = attrs.find(needle)
    if i < 0:
        return None
    j = attrs.find('"', i + len(needle))
    if j < 0:
        return None
    return attrs[i + len(needle):j]


def parse_gtf(gtf_path):
    """Return (intervals_by_chrom, gene_to_family, gene_to_class, level).

    intervals_by_chrom: {chrom: list[(start, end, gene_id)]}, indexed at the
        broadest level present in the GTF: gene > transcript > exon.
    gene_to_family, gene_to_class: {gene_id: label} (empty when absent).
    level: which feature_type was chosen ('gene', 'transcript', or 'exon').

    Indexing at the gene level includes introns; a read in any intron of
    gene G counts toward G as long as no other gene's interval overlaps.
    """
    by_type = {'gene': defaultdict(list),
               'transcript': defaultdict(list),
               'exon': defaultdict(list)}
    gene_to_family = {}
    gene_to_class = {}

    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in by_type:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GTF 1-based inclusive -> 0-based half-open
            end = int(parts[4])
            if end <= start:
                continue
            attrs = parts[8]
            gene_id = parse_attribute(attrs, 'gene_id')
            if gene_id is None:
                continue
            by_type[feature_type][chrom].append((start, end, gene_id))
            fam = parse_attribute(attrs, 'family_id')
            cls = parse_attribute(attrs, 'class_id')
            if fam is not None:
                gene_to_family[gene_id] = fam
            if cls is not None:
                gene_to_class[gene_id] = cls

    for level in ('gene', 'transcript', 'exon'):
        if by_type[level]:
            chosen = by_type[level]
            n_intervals = sum(len(v) for v in chosen.values())
            print(
                f'GTF indexed at level={level}: {n_intervals} intervals across '
                f'{len(chosen)} chromosomes; family_id known for '
                f'{len(gene_to_family)} genes, class_id for '
                f'{len(gene_to_class)}', file=sys.stderr)
            return chosen, gene_to_family, gene_to_class, level

    sys.exit(f'no gene/transcript/exon records found in {gtf_path}')


def load_whitelist(path):
    cbs = []
    with open_text(path) as fh:
        for line in fh:
            cb = line.rstrip('\n').split('\t', 1)[0]
            if cb:
                cbs.append(cb)
    if len(set(cbs)) != len(cbs):
        sys.exit('whitelist contains duplicate barcodes')
    print(f'whitelist size: {len(cbs)}', file=sys.stderr)
    return cbs


_WORKER_BAM = None
_WORKER_TREE = None
_WORKER_CB_INDEX = None
_WORKER_MULTIMAPPER = None
_WORKER_UMI_DEDUP = None
_WORKER_CB_TAG = None


def _worker_init(bam_path, intervals_by_chrom, cb_index, multimapper,
                  umi_dedup, cb_tag):
    import pysam
    from intervaltree import IntervalTree
    global _WORKER_BAM, _WORKER_TREE, _WORKER_CB_INDEX
    global _WORKER_MULTIMAPPER, _WORKER_UMI_DEDUP, _WORKER_CB_TAG
    _WORKER_BAM = pysam.AlignmentFile(bam_path, 'rb')
    _WORKER_TREE = {}
    for chrom, recs in intervals_by_chrom.items():
        tree = IntervalTree()
        for s, e, gid in recs:
            tree.addi(s, e, gid)
        _WORKER_TREE[chrom] = tree
    _WORKER_CB_INDEX = cb_index
    _WORKER_MULTIMAPPER = multimapper
    _WORKER_UMI_DEDUP = umi_dedup
    _WORKER_CB_TAG = cb_tag


_NO_OVERLAP = object()
_AMBIGUOUS = object()


def _gene_id_for_read(chrom, blocks):
    """Resolve the gene_id for an alignment given its aligned blocks.

    `blocks` is the list of (start, end) tuples returned by
    pysam.AlignedSegment.get_blocks(): one entry per CIGAR M/=/X run, with
    skipped intronic regions (CIGAR N) NOT included. We union overlaps
    across blocks so a spliced read counts toward the gene if every aligned
    block sits inside (or touches) that gene's interval — and any overlap
    with a second gene's interval makes the read ambiguous.

    Returns:
      gene_id (str)   -- exactly one gene_id overlaps any block
      _NO_OVERLAP     -- chrom not in the tree, or no block touches any feature
      _AMBIGUOUS      -- two or more distinct gene_ids overlap aligned blocks
    """
    tree = _WORKER_TREE.get(chrom)
    if tree is None:
        return _NO_OVERLAP
    gene_ids = set()
    for s, e in blocks:
        for iv in tree.overlap(s, e):
            gene_ids.add(iv.data)
            if len(gene_ids) > 1:
                return _AMBIGUOUS
    if not gene_ids:
        return _NO_OVERLAP
    return next(iter(gene_ids))


def _process_contig(chrom):
    """Walk reads on one contig and return per-key UMI weights.

    Per-(cb_idx, gene_id) state is dict[umi -> total weight], where weight
    is 1.0 for unique-mappers and 1.0/NH for multi-mappers (when
    --multimapper multi). When the same UMI is observed in multiple
    records for the same (cb, gene_id), we keep the maximum weight: a
    UMI seen as a unique-mapper trumps the same UMI seen as a multi-mapper
    on a sibling alignment record.
    """
    counts = defaultdict(dict)
    n_total = n_kept = n_no_cb = n_no_ub = n_off_wl = 0
    n_no_feature = n_ambiguous = n_multi_skipped = 0
    bf = _WORKER_BAM
    cb_idx = _WORKER_CB_INDEX
    mm = _WORKER_MULTIMAPPER
    has_chrom = chrom in _WORKER_TREE
    for r in bf.fetch(contig=chrom):
        n_total += 1
        if r.is_unmapped or r.is_supplementary:
            continue
        # Paired BAMs (typical SmartSeq2): count R1 only so a fragment is
        # not double-counted at both ends. Single-end BAMs have is_paired
        # False, so the filter is a no-op there.
        if r.is_paired and r.is_read2:
            continue
        # STAR emits multimappers as one primary + (NH-1) secondary records,
        # all sharing the same UB. unique mode keeps only the primary;
        # multi mode keeps every locus so the per-record 1/NH weights sum
        # toward 1.0 across the molecule's reference loci.
        if r.is_secondary and mm == 'unique':
            continue
        try:
            cb = r.get_tag(_WORKER_CB_TAG)
        except KeyError:
            n_no_cb += 1
            continue
        cb_i = cb_idx.get(cb)
        if cb_i is None:
            n_off_wl += 1
            continue
        if _WORKER_UMI_DEDUP == 'none':
            # No real UMI: use a per-record synthetic key so each BAM record
            # is distinct in the (cb, gid) bucket. (qname, ref_start) is
            # unique per record because different reads have distinct qnames,
            # while primary + secondary records of the same multi-mapper
            # share qname but sit at different ref_start.
            ub = (r.query_name, r.reference_start)
        else:
            try:
                ub = r.get_tag('UB')
            except KeyError:
                n_no_ub += 1
                continue
        nh = r.get_tag('NH') if r.has_tag('NH') else 1
        if nh > 1 and mm == 'unique':
            n_multi_skipped += 1
            continue
        weight = 1.0 if nh == 1 else 1.0 / nh
        if not has_chrom:
            n_no_feature += 1
            continue
        # get_blocks() excludes CIGAR-N gaps (introns) so a spliced read
        # only counts as overlapping features its exonic blocks actually
        # cover. Falls back to the full span when get_blocks() returns
        # nothing (e.g. soft-clipped-only alignments, which shouldn't
        # arise but defended against).
        blocks = r.get_blocks() or [(r.reference_start, r.reference_end)]
        gid = _gene_id_for_read(chrom, blocks)
        if gid is _NO_OVERLAP:
            n_no_feature += 1
            continue
        if gid is _AMBIGUOUS:
            n_ambiguous += 1
            continue
        umi_map = counts[(cb_i, gid)]
        prev = umi_map.get(ub)
        if prev is None or weight > prev:
            umi_map[ub] = weight
        n_kept += 1
    stats = (n_total, n_kept, n_no_cb, n_no_ub, n_off_wl,
             n_multi_skipped, n_no_feature, n_ambiguous)
    return chrom, dict(counts), stats


def merge_counts(per_contig):
    """Merge per-contig {(cb_i, gene_id) -> {umi -> weight}} dicts. Distinct
    gene_ids per chrom usually means no key collisions; UMI dicts are merged
    by max-weight per UMI when collisions do occur."""
    merged = {}
    for partial in per_contig:
        for k, v in partial.items():
            existing = merged.get(k)
            if existing is None:
                merged[k] = v
            else:
                for umi, w in v.items():
                    if umi not in existing or w > existing[umi]:
                        existing[umi] = w
    return merged


def _hamming_at_most_one(a, b):
    """True if a and b differ at <= 1 position. Equal-length only."""
    if a == b:
        return True
    if len(a) != len(b):
        return False
    diff = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            diff += 1
            if diff > 1:
                return False
    return True


def dedup_exact(umi_weights):
    """Sum of weights across distinct UMI strings."""
    return sum(umi_weights.values())


def dedup_1mm_all(umi_weights):
    """Collapse UMIs within Hamming distance 1 via union-find. The deduped
    count is the sum of max-weight per connected component, mirroring
    STARsolo's 1MM_All semantics where one count is attributed to each
    component (the parent UMI's weight is conserved).

    For unique-mode (all weights = 1.0), this returns the integer number
    of components, matching STARsolo's default.
    """
    umis = list(umi_weights.keys())
    n = len(umis)
    if n == 0:
        return 0.0
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if _hamming_at_most_one(umis[i], umis[j]):
                pi, pj = find(i), find(j)
                if pi != pj:
                    parent[pi] = pj

    component_max = {}
    for i, umi in enumerate(umis):
        root = find(i)
        w = umi_weights[umi]
        if component_max.get(root, -1.0) < w:
            component_max[root] = w
    return sum(component_max.values())


def dedup_none(umi_weights):
    """No UMI deduplication: sum every weight in the bucket. Pair with
    --umi-dedup none, which makes _process_contig key the bucket by a
    per-record synthetic UB so every BAM record contributes its weight.

    Use case: SmartSeq2 BAMs without real UMIs. Paired-end BAMs are
    handled upstream by skipping is_read2 records so a fragment is counted
    once. The pipeline does not attempt duplicate removal: with no UMI it
    cannot distinguish technical PCR duplicates from biologically real
    reads at the same coordinates.
    """
    return sum(umi_weights.values())


DEDUP_FUNCS = {
    'exact': dedup_exact,
    '1mm_all': dedup_1mm_all,
    'none': dedup_none,
}


def collapse(counts, method):
    """Apply UMI dedup to every (cb, gene_id) bucket and return scalar counts."""
    fn = DEDUP_FUNCS[method]
    return {k: fn(v) for k, v in counts.items()}


def write_matrix(counts, n_cells, out_dir, mtx_name):
    """Write Matrix Market matrix.mtx + features.tsv. barcodes.tsv is
    written by the caller (shared across granularities). counts values
    are scalar; integer values are written as int, fractional as float
    so STARsolo-style readers parse them transparently."""
    os.makedirs(out_dir, exist_ok=True)
    feature_set = sorted({gid for _, gid in counts.keys()})
    feature_index = {gid: i for i, gid in enumerate(feature_set)}
    with open(os.path.join(out_dir, 'features.tsv'), 'w') as fh:
        for gid in feature_set:
            fh.write(f'{gid}\t{gid}\tGene Expression\n')
    entries = []
    for (cb_i, gid), n in counts.items():
        if n == 0:
            continue
        entries.append((feature_index[gid] + 1, cb_i + 1, n))
    entries.sort()
    with open(os.path.join(out_dir, mtx_name), 'w') as fh:
        fh.write('%%MatrixMarket matrix coordinate real general\n')
        fh.write('%\n')
        fh.write(f'{len(feature_set)} {n_cells} {len(entries)}\n')
        for fi, ci, n in entries:
            if isinstance(n, float) and not n.is_integer():
                fh.write(f'{fi} {ci} {n:.6f}\n')
            else:
                fh.write(f'{fi} {ci} {int(n)}\n')
    print(f'wrote {os.path.join(out_dir, mtx_name)}', file=sys.stderr)


def rollup(counts, gene_to_group):
    """Sum already-deduped scalar counts up to coarser groups. Genes with
    no group label are dropped."""
    rolled = defaultdict(float)
    for (cb_i, gid), n in counts.items():
        group = gene_to_group.get(gid)
        if group is None:
            continue
        rolled[(cb_i, group)] += n
    return rolled


def write_outputs(counts, cbs, out_root, multimapper, umi_dedup,
                   gene_to_family, gene_to_class):
    """counts is {(cb_i, gene_id) -> scalar} (UMI dedup already applied).
    Family/class rollups sum the gene-level scalars."""
    n_cells = len(cbs)
    mtx_name = ('UniqueAndMult-EM.mtx' if multimapper == 'multi'
                else 'matrix.mtx')

    def write_barcodes(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as fh:
            for cb in cbs:
                fh.write(cb + '\n')

    gene_dir = os.path.join(out_root, 'Solo.out', 'Gene', 'raw')
    write_barcodes(os.path.join(gene_dir, 'barcodes.tsv'))
    write_matrix(counts, n_cells, gene_dir, mtx_name)

    if gene_to_family:
        fam_dir = os.path.join(out_root, 'Solo.out', 'Gene_family_id', 'raw')
        write_barcodes(os.path.join(fam_dir, 'barcodes.tsv'))
        write_matrix(rollup(counts, gene_to_family), n_cells, fam_dir, mtx_name)

    if gene_to_class:
        cls_dir = os.path.join(out_root, 'Solo.out', 'Gene_class_id', 'raw')
        write_barcodes(os.path.join(cls_dir, 'barcodes.tsv'))
        write_matrix(rollup(counts, gene_to_class), n_cells, cls_dir, mtx_name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--whitelist', required=True)
    ap.add_argument('--out-dir', required=True,
                    help='Output root; Solo.out/Gene[_<granularity>]/raw/ '
                         'subdirs are created beneath it')
    ap.add_argument('--multimapper', choices=('unique', 'multi'),
                    default='unique',
                    help='unique (default) skips NH > 1 reads; multi '
                         'includes them with weight 1/NH (approximate)')
    ap.add_argument('--umi-dedup', choices=tuple(DEDUP_FUNCS.keys()),
                    default='1mm_all',
                    help='1mm_all (default, matches STARsolo) collapses '
                         'UMIs within Hamming-1; exact counts distinct UB '
                         'strings (debug only, over-counts)')
    ap.add_argument('--cb-tag', default='CB',
                    help='BAM tag carrying the cell identifier. CB '
                         '(default) for Chromium and bowtie2 chromium '
                         'output. RG for STARsolo SmartSeq output, where '
                         'CB is forbidden by STAR and the cell-id sits in '
                         'the read-group tag.')
    ap.add_argument('--threads', type=int, default=1)
    args = ap.parse_args()

    cbs = load_whitelist(args.whitelist)
    cb_index = {cb: i for i, cb in enumerate(cbs)}
    intervals_by_chrom, gene_to_family, gene_to_class, level = parse_gtf(args.gtf)

    import pysam
    bam_header = pysam.AlignmentFile(args.bam, 'rb')
    contigs = list(bam_header.references)
    bam_header.close()

    init_args = (args.bam, intervals_by_chrom, cb_index, args.multimapper, args.umi_dedup, args.cb_tag)
    if args.threads <= 1:
        _worker_init(*init_args)
        results = [_process_contig(c) for c in contigs]
    else:
        ctx = mp.get_context('fork')
        with ctx.Pool(processes=args.threads, initializer=_worker_init,
                      initargs=init_args) as pool:
            results = pool.map(_process_contig, contigs)

    totals = [0] * 8
    per_contig = []
    for chrom, partial, stats in results:
        per_contig.append(partial)
        for i, v in enumerate(stats):
            totals[i] += v
    raw_counts = merge_counts(per_contig)
    counts = collapse(raw_counts, args.umi_dedup)

    n_total, n_kept, n_no_cb, n_no_ub, n_off_wl, n_multi, n_nf, n_amb = totals
    print(
        f'GTF indexing level     : {level}\n'
        f'multimapper            : {args.multimapper}\n'
        f'umi dedup              : {args.umi_dedup}\n'
        f'reads scanned: {n_total}\n'
        f'  kept (counted)         : {n_kept}\n'
        f'  no CB tag              : {n_no_cb}\n'
        f'  no UB tag              : {n_no_ub}\n'
        f'  CB not in whitelist    : {n_off_wl}\n'
        f'  multimapper skipped    : {n_multi}\n'
        f'  no overlapping feature : {n_nf}\n'
        f'  multi-feature ambiguous: {n_amb}',
        file=sys.stderr)

    write_outputs(counts, cbs, args.out_dir, args.multimapper,
                   args.umi_dedup, gene_to_family, gene_to_class)


if __name__ == '__main__':
    main()
