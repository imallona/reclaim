#!/usr/bin/env python
"""
Phase 1 selection for the repogle_2022 use case.

Reads the figshare K562 essential bulk and singlecell h5ads (Replogle et al.
2022) and the KD6 raw-files manifest. Picks the top N perturbations by
Replogle's published per-perturbation TE_ratio, samples a control pool, and
emits three slim TSVs that downstream phases consume:

  selected_perturbations.tsv  - top perturbations + their replogle_te_ratio
  cells_to_perturbation.tsv   - per-cell join key (gem_group, cell_barcode,
                                 gene_transcript, sgID_AB, umi_count)
  selected_libraries.tsv      - per-gemgroup library metadata (mRNA fastqs,
                                 sgRNA fastqs, BAM filename) parsed from the
                                 KD6 manifest

Inputs are the raw files; only the slim TSVs are tracked in the repo.
"""

import argparse
import csv
import h5py
import numpy as np
import sys


def decode_obj_array(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def load_bulk(path, min_cells):
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        gene_transcript = decode_obj_array(obs["gene_transcript"][:])
        te_ratio = obs["TE_ratio"][:]
        core_control = obs["core_control"][:]
        num_cells = obs["num_cells_unfiltered"][:]
    keep = (~core_control) & (num_cells >= min_cells)
    return gene_transcript[keep], te_ratio[keep], num_cells[keep]


def pick_top_n(gene_transcript, te_ratio, num_cells, top_n):
    order = np.argsort(-te_ratio)[:top_n]
    return [
        {
            "gene_transcript": gene_transcript[i],
            "replogle_te_ratio": float(te_ratio[i]),
            "n_cells_replogle": int(num_cells[i]),
        }
        for i in order
    ]


def load_singlecell_obs(path):
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        cats = obs["__categories"]
        codes_gene_transcript = obs["gene_transcript"][:]
        codes_gene = obs["gene"][:]
        codes_sgID = obs["sgID_AB"][:]
        levels_gene_transcript = decode_obj_array(cats["gene_transcript"][:])
        levels_gene = decode_obj_array(cats["gene"][:])
        levels_sgID = decode_obj_array(cats["sgID_AB"][:])
        cell_barcode = decode_obj_array(obs["cell_barcode"][:])
        gem_group = obs["gem_group"][:]
        umi = obs["UMI_count"][:]
    return {
        "cell_barcode": cell_barcode,
        "gem_group": gem_group,
        "gene": levels_gene[codes_gene],
        "sgID_AB": levels_sgID[codes_sgID],
        "gene_transcript": levels_gene_transcript[codes_gene_transcript],
        "umi_count": umi,
    }


def select_cells(sc, top_perts, n_controls, gemgroup_subset, seed):
    gt = sc["gene_transcript"]
    target_set = set(p["gene_transcript"] for p in top_perts)
    is_top = np.isin(gt, list(target_set))
    is_control = np.array([
        "non-targeting" in x.lower() for x in gt
    ])

    if gemgroup_subset:
        gg_keep = np.isin(sc["gem_group"], gemgroup_subset)
        is_top = is_top & gg_keep
        is_control = is_control & gg_keep

    rng = np.random.default_rng(seed)
    control_idx = np.where(is_control)[0]
    n_pick = min(n_controls, len(control_idx))
    sampled = rng.choice(control_idx, size=n_pick, replace=False) if n_pick else np.array([], dtype=int)

    keep = np.zeros(len(gt), dtype=bool)
    keep[is_top] = True
    keep[sampled] = True
    return keep, int(is_top.sum()), int(n_pick)


def parse_manifest(path, gemgroup_subset):
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gg = int(row["gemgroup"])
            if gemgroup_subset and gg not in gemgroup_subset:
                continue
            mrna = []
            sgrna = []
            for k, v in row.items():
                if v is None or k in {"library", "gemgroup", "bam", "bam_index",
                                       "barcodes", "features", "matrix"}:
                    continue
                if not v.endswith(".fastq.gz"):
                    continue
                if "_mRNA_" in v:
                    mrna.append(v)
                elif "_sgRNA_" in v:
                    sgrna.append(v)
            rows.append({
                "gem_group": gg,
                "library": row["library"],
                "bam": row.get("bam", ""),
                "barcodes_tsv": row.get("barcodes", ""),
                "features_tsv": row.get("features", ""),
                "matrix_mtx": row.get("matrix", ""),
                "n_mrna_fastqs": len(mrna),
                "n_sgrna_fastqs": len(sgrna),
                "mrna_fastqs": ";".join(sorted(mrna)),
                "sgrna_fastqs": ";".join(sorted(sgrna)),
            })
    rows.sort(key=lambda r: r["gem_group"])
    return rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def parse_gemgroup_subset(s):
    if not s:
        return []
    return sorted({int(x) for x in s.split(",") if x.strip()})


def read_gene_list(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sym = line.split("\t")[0].strip()
            if sym.lower() == "gene":
                continue
            genes.append(sym)
    return set(genes)


def pick_by_gene_list(sc, gt_bulk, te_ratio, gene_list_path):
    wanted = read_gene_list(gene_list_path)
    gt = sc["gene_transcript"]
    gn = sc["gene"]
    is_control = np.array(["non-targeting" in x.lower() for x in gt])
    mask = np.isin(gn, list(wanted)) & ~is_control
    sel_gt, counts = np.unique(gt[mask], return_counts=True)
    te_map = dict(zip(gt_bulk, te_ratio))
    return [
        {
            "gene_transcript": g,
            "replogle_te_ratio": float(te_map.get(g, float("nan"))),
            "n_cells_replogle": int(c),
        }
        for g, c in zip(sel_gt, counts)
    ]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bulk", required=True)
    p.add_argument("--singlecell", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--top-n", type=int, default=0,
                   help="top N perturbations by TE_ratio; ignored when "
                        "--gene-list is given")
    p.add_argument("--gene-list", default="",
                   help="file with one gene symbol per line (or a tsv whose "
                        "first column is the symbol); selects those genes' "
                        "perturbations instead of top-N by TE_ratio")
    p.add_argument("--n-controls", type=int, required=True)
    p.add_argument("--min-cells", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--gemgroups", default="",
                   help="comma-separated gemgroup ids; empty = all")
    p.add_argument("--out-perturbations", required=True)
    p.add_argument("--out-cells", required=True)
    p.add_argument("--out-libraries", required=True)
    args = p.parse_args()

    gemgroup_subset = parse_gemgroup_subset(args.gemgroups)

    print("Loading bulk h5ad", file=sys.stderr)
    gt_bulk, te_ratio, num_cells = load_bulk(args.bulk, args.min_cells)

    print("Loading singlecell h5ad obs", file=sys.stderr)
    sc = load_singlecell_obs(args.singlecell)
    print(f"  {len(sc['cell_barcode'])} cells, "
          f"{len(np.unique(sc['gem_group']))} gemgroups", file=sys.stderr)

    if args.gene_list:
        top_perts = pick_by_gene_list(sc, gt_bulk, te_ratio, args.gene_list)
        print(f"  gene-list selection: {len(top_perts)} perturbations",
              file=sys.stderr)
    else:
        top_perts = pick_top_n(gt_bulk, te_ratio, num_cells, args.top_n)
        print(f"  picked top {len(top_perts)} perturbations "
              f"(min_cells>={args.min_cells})", file=sys.stderr)

    keep, n_top, n_ctrl = select_cells(
        sc, top_perts, args.n_controls, gemgroup_subset, args.seed)
    print(f"  selected {n_top} perturbed + {n_ctrl} control cells "
          f"(total {keep.sum()})", file=sys.stderr)

    print("Parsing KD6 manifest", file=sys.stderr)
    library_rows = parse_manifest(args.manifest, gemgroup_subset)
    print(f"  {len(library_rows)} libraries", file=sys.stderr)

    write_tsv(args.out_perturbations, top_perts,
              ["gene_transcript", "replogle_te_ratio", "n_cells_replogle"])

    cell_rows = []
    for i in np.where(keep)[0]:
        cell_rows.append({
            "cell_barcode": sc["cell_barcode"][i],
            "gem_group": int(sc["gem_group"][i]),
            "gene": sc["gene"][i],
            "sgID_AB": sc["sgID_AB"][i],
            "gene_transcript": sc["gene_transcript"][i],
            "umi_count": int(sc["umi_count"][i]),
        })
    cell_rows.sort(key=lambda r: (r["gem_group"], r["cell_barcode"]))
    write_tsv(args.out_cells, cell_rows,
              ["cell_barcode", "gem_group", "gene", "sgID_AB",
               "gene_transcript", "umi_count"])

    write_tsv(args.out_libraries, library_rows,
              ["gem_group", "library", "bam", "barcodes_tsv", "features_tsv",
               "matrix_mtx", "n_mrna_fastqs", "n_sgrna_fastqs",
               "mrna_fastqs", "sgrna_fastqs"])

    print("Done", file=sys.stderr)


if __name__ == "__main__":
    main()
