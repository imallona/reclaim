#!/usr/bin/env python3
"""
Harmonize per-tool external-benchmark counts to the REclaim feature x cell TSV
format consumed by evaluate.py.

Per-tool input formats:
  - tetranscripts: a per-sample TSV with one column per sample and one row
    per "subfamily:family:class" feature_id. Stripped to the subfamily
    component (matches REclaim gene_id granularity).
  - scte: an h5ad cell x repeat AnnData. Repeat features are family-level
    by default (matches REclaim family_id granularity).
  - squire: a per-locus TSV with TE_ID, TE_name, family, class, subfamily,
    counts. Aggregated to the requested granularity.

The harmonized output is a feature_id x cell_id TSV, identical in shape to
counts/{aligner}_{feature_set}_{granularity}_{multimapper_mode}.tsv. evaluate.py
already understands that shape and rolls up to {gene_id, family_id, class_id}
internally as needed; this script's only job is to put the external tool's
output into a TSV with the correct rownames and colnames.

Asserts at boundaries:
  - Input file readable and non-empty.
  - At least one column matches the user-supplied sample list.
  - Feature-ID overlap with the locus_map exceeds --min-overlap-fraction.
  - Output is non-empty.

Drops are logged to stderr with counts so the user can see how many features
fell off due to no mapping.

Pure helper map_features_to_reclaim is exported for unit testing.
"""
import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict

TOOLS = ('tetranscripts', 'scte', 'squire')
GRANULARITIES = ('gene_id', 'family_id', 'class_id')


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path)


def parse_locus_map(locus_map_path):
    """Return three dicts: gene_id->class_id, family_id->class_id, gene_id->family_id.

    The locus_map is the canonical 4-column TSV produced by reference.snmk:
    transcript_id, gene_id, family_id, class_id. We deduplicate in this
    function and warn on the (unexpected) case where one gene maps to
    multiple classes.
    """
    if not locus_map_path or not os.path.exists(locus_map_path):
        return {}, {}, {}
    gene_to_class = {}
    family_to_class = {}
    gene_to_family = {}
    conflicts = []
    with open_text(locus_map_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            _, gene_id, family_id, class_id = parts[:4]
            if not gene_id or gene_id == 'gene_id':  # header skip
                continue
            if gene_id in gene_to_class and gene_to_class[gene_id] != class_id:
                conflicts.append(gene_id)
            gene_to_class[gene_id] = class_id
            family_to_class[family_id] = class_id
            gene_to_family[gene_id] = family_id
    if conflicts:
        sys.stderr.write(
            f'  warning: {len(set(conflicts))} gene_id(s) map to multiple class_ids; '
            f'first wins (e.g. {sorted(set(conflicts))[:3]})\n')
    return gene_to_class, family_to_class, gene_to_family


def map_features_to_reclaim(external_ids, reclaim_ids, strip_suffix=None):
    """Pure helper: map an external tool's feature IDs to REclaim's IDs.

    Strategy: try identity first; for IDs that don't match, try stripping a
    user-supplied suffix delimiter (e.g. ':' for TEtranscripts'
    'AluY:Alu:SINE'). Returns a dict {external_id: reclaim_id_or_None}.

    The function is intentionally simple so it stays testable. Tool-specific
    parsing (which delimiter, which field) lives in the caller.
    """
    assert isinstance(external_ids, (list, tuple, set))
    reclaim_set = set(reclaim_ids)
    out = {}
    for eid in external_ids:
        if eid in reclaim_set:
            out[eid] = eid
            continue
        if strip_suffix is not None and strip_suffix in eid:
            head = eid.split(strip_suffix, 1)[0]
            out[eid] = head if head in reclaim_set else None
        else:
            out[eid] = None
    return out


def aggregate_to_granularity(counts_per_feature, gene_to_family, gene_to_class,
                             granularity):
    """Roll a {gene_id: {sample: count}} dict up to family_id or class_id.

    For granularity='gene_id' this is a no-op (returns the input). For
    family_id, sums on gene_to_family. For class_id, sums on gene_to_class.
    Features without a known mapping are dropped and counted.
    """
    assert granularity in GRANULARITIES
    if granularity == 'gene_id':
        return counts_per_feature, 0
    mapping = gene_to_family if granularity == 'family_id' else gene_to_class
    out = defaultdict(lambda: defaultdict(float))
    dropped = 0
    for gene_id, per_sample in counts_per_feature.items():
        rolled_key = mapping.get(gene_id)
        if rolled_key is None:
            dropped += 1
            continue
        for sample, count in per_sample.items():
            out[rolled_key][sample] += count
    return {k: dict(v) for k, v in out.items()}, dropped


def parse_tetranscripts_table(path, sample_columns):
    """TEtranscripts writes a header row with sample columns and rows
    feature_id\tcount1\tcount2\\... where feature_id is 'subfamily:family:class'.
    Returns {gene_id: {sample: count}} keyed on the subfamily component.
    """
    counts = defaultdict(dict)
    with open_text(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        feature_col = header[0]
        del feature_col  # informational
        # Map declared sample order to the file's columns. We require all
        # requested samples to be present.
        col_idx = {}
        for s in sample_columns:
            if s not in header:
                raise AssertionError(
                    f'sample {s!r} not found in TEtranscripts header {header!r}')
            col_idx[s] = header.index(s)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts or not parts[0]:
                continue
            raw_id = parts[0]
            subfamily = raw_id.split(':', 1)[0] if ':' in raw_id else raw_id
            for s, ix in col_idx.items():
                try:
                    v = float(parts[ix])
                except (ValueError, IndexError):
                    v = 0.0
                if v > 0:
                    counts[subfamily][s] = v
    return dict(counts)


def parse_scte_h5ad(path, sample_columns):
    """scTE emits a cell x feature AnnData. We expect one h5ad per
    'sample' (here: one per simulated cell), and concatenate them. Returns
    {family_id: {cell_id: count}}.

    For the simulation benchmark the cell ids inside each h5ad's obs index
    are passed straight through; we assume the user supplies them as
    sample_columns when the technology is single cell.

    Imports anndata lazily so the rest of the script is usable without it.
    """
    try:
        import anndata as ad
    except ImportError as e:
        raise SystemExit(f'scte path requires anndata; install or use --tool other: {e}')
    counts = defaultdict(dict)
    if os.path.isdir(path):
        h5_files = sorted(
            os.path.join(path, f) for f in os.listdir(path)
            if f.endswith('.h5ad'))
    else:
        h5_files = [path]
    assert h5_files, f'no .h5ad files at {path}'
    sample_set = set(sample_columns)
    for h5 in h5_files:
        a = ad.read_h5ad(h5)
        # Repeats columns: scTE prepends 'TE_' or stores under .var with
        # 'is_TE' True; conservative path is to use all columns and let the
        # downstream locus_map mapping pick what is recognisable.
        var_names = list(a.var_names)
        for cell_idx, cell_id in enumerate(a.obs_names):
            if sample_set and cell_id not in sample_set:
                continue
            row = a.X[cell_idx]
            try:
                row_arr = row.toarray().ravel()
            except AttributeError:
                row_arr = row
            for fi, val in enumerate(row_arr):
                v = float(val)
                if v > 0:
                    counts[var_names[fi]][cell_id] = v
    return dict(counts)


def parse_squire_count(path, sample_columns):
    """SQuIRE per-locus output: TE_ID, TE_name, family, class, subfamily,
    counts (one per sample). For benchmark use we aggregate locus rows up
    to the subfamily column so granularity matches REclaim gene_id.
    Returns {subfamily: {sample: count}}.
    """
    counts = defaultdict(lambda: defaultdict(float))
    with open_text(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            sub_col = header.index('subfamily')
        except ValueError:
            raise AssertionError(f'subfamily column missing in {path!r}')
        col_idx = {}
        for s in sample_columns:
            if s not in header:
                raise AssertionError(f'sample {s!r} not in SQuIRE header')
            col_idx[s] = header.index(s)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts or not parts[0]:
                continue
            sub = parts[sub_col]
            for s, ix in col_idx.items():
                try:
                    v = float(parts[ix])
                except (ValueError, IndexError):
                    v = 0.0
                if v > 0:
                    counts[sub][s] += v
    return {k: dict(v) for k, v in counts.items()}


def _fmt_count(v):
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def write_count_matrix(counts, sample_columns, out_path):
    """Writes a feature_id x cell_id TSV. Zero-count cells are written as 0
    (the REclaim count tables include the 0s; evaluate.py only reads non-
    zero entries from observed but also walks all rows for the feature
    universe). Whole-number counts are rendered as integers to match the
    rest of the pipeline's count TSVs.
    """
    features = sorted(counts.keys())
    assert features, 'cannot write an empty count matrix'
    with open(out_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['feature_id'] + list(sample_columns))
        for f in features:
            row = [f] + [_fmt_count(counts[f].get(s, 0)) for s in sample_columns]
            w.writerow(row)


def check_overlap(features_external, features_reclaim, threshold, label):
    if not features_reclaim:
        sys.stderr.write(
            f'  warning: locus_map empty; cannot check overlap for {label}\n')
        return
    overlap = len(set(features_external) & set(features_reclaim))
    frac = overlap / len(set(features_external)) if features_external else 0
    sys.stderr.write(
        f'  {label}: {overlap}/{len(set(features_external))} '
        f'({frac:.1%}) external features matched a REclaim {label} key\n')
    assert frac >= threshold, (
        f'feature overlap {frac:.1%} below required {threshold:.0%} '
        f'for {label}; check tool output and locus_map alignment')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tool', choices=TOOLS, required=True)
    ap.add_argument('--input', required=True,
                    help='Path to tool output (file or directory; see per-tool parser).')
    ap.add_argument('--locus-map', required=True,
                    help='4-col REclaim locus_map TSV (transcript_id, gene_id, '
                         'family_id, class_id).')
    ap.add_argument('--samples', required=True,
                    help='Comma-separated list of cell/sample IDs to keep.')
    ap.add_argument('--granularity', choices=GRANULARITIES, default='gene_id',
                    help='Output granularity. Tool outputs at one native level '
                         'and we roll up here; for class_id this also relies on '
                         'the locus_map.')
    ap.add_argument('--min-overlap-fraction', type=float, default=0.05,
                    help='Asserts that this fraction of the tool feature IDs '
                         'maps to a known REclaim feature. Default 0.05.')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    samples = [s for s in args.samples.split(',') if s]
    assert samples, '--samples must include at least one sample id'
    assert os.path.exists(args.input), f'input not found: {args.input}'

    gene_to_class, family_to_class, gene_to_family = parse_locus_map(args.locus_map)

    if args.tool == 'tetranscripts':
        # TEtranscripts feature_id: 'subfamily:family:class' (matches REclaim gene_id)
        per_feature = parse_tetranscripts_table(args.input, samples)
        check_overlap(list(per_feature.keys()), list(gene_to_class.keys()),
                      args.min_overlap_fraction, 'gene_id (subfamily)')
        rolled, dropped = aggregate_to_granularity(
            per_feature, gene_to_family, gene_to_class, args.granularity)
    elif args.tool == 'scte':
        # scTE feature_id is family-level by default
        per_feature = parse_scte_h5ad(args.input, samples)
        check_overlap(list(per_feature.keys()), list(family_to_class.keys()),
                      args.min_overlap_fraction, 'family_id')
        if args.granularity == 'gene_id':
            sys.stderr.write(
                '  warning: scTE outputs at family granularity; gene_id '
                'rollup is not possible (returns family-level rows '
                'as-is and evaluate.py will treat them as gene_ids)\n')
            rolled, dropped = per_feature, 0
        elif args.granularity == 'family_id':
            rolled, dropped = per_feature, 0
        else:  # class_id
            rolled = defaultdict(lambda: defaultdict(float))
            dropped = 0
            for fam, per_sample in per_feature.items():
                cls = family_to_class.get(fam)
                if cls is None:
                    dropped += 1
                    continue
                for s, c in per_sample.items():
                    rolled[cls][s] += c
            rolled = {k: dict(v) for k, v in rolled.items()}
    elif args.tool == 'squire':
        per_feature = parse_squire_count(args.input, samples)
        check_overlap(list(per_feature.keys()), list(gene_to_class.keys()),
                      args.min_overlap_fraction, 'gene_id (subfamily)')
        rolled, dropped = aggregate_to_granularity(
            per_feature, gene_to_family, gene_to_class, args.granularity)
    else:
        raise SystemExit(f'unsupported tool {args.tool!r}')

    if dropped:
        sys.stderr.write(
            f'  rolled to {args.granularity}: dropped {dropped} feature(s) '
            f'with no mapping\n')

    assert rolled, 'rolled count matrix empty after harmonization'
    write_count_matrix(rolled, samples, args.output)
    sys.stderr.write(
        f'  wrote {len(rolled)} features x {len(samples)} samples to '
        f'{args.output}\n')


if __name__ == '__main__':
    main()
