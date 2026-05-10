# repogle_2022

Reanalysis of the K562 essential Perturb-seq screen from Replogle et al. 2022, Cell.

## Motivation

Replogle et al. report a per-perturbation `TE_ratio` (fraction of UMIs in repetitive elements) and identify exosome, CPSF/Integrator/Pol II machinery as the dominant TE-fraction-increasing perturbation classes. Their metric is class-level only (Alu, L1, MIR aggregated). REclaim quantifies repeats at locus, gene_id, family_id, and class_id granularity, so the same dataset supports a higher-resolution view: subfamily- or locus-level fingerprints per perturbation, and validation of REclaim's pseudobulk counts against the published `TE_ratio`.

## Data sources

Two figshare deposits (Replogle and Weissman, CC-BY 4.0):

| Article | DOI | Contents |
|---|---|---|
| 20022944 | 10.25452/figshare.plus.20022944 | SRA/GEO file manifests (`KD6/KD8/KD8_ultima/RD7_raw_files.csv`) |
| 20029387 | 10.25452/figshare.plus.20029387 | Processed h5ad datasets (raw and normalized, bulk and singlecell, K562 essential, K562 gwps, RPE1) |

Three files are needed for the K562 essential pilot. URLs and md5 checksums are pinned in `paper/configs/repogle_2022.yaml`:

| File | Size | Purpose |
|---|---|---|
| `KD6_raw_files.csv` | 76 KB | per-library fastq and BAM filenames, gemgroup ids |
| `K562_essential_raw_bulk_01.h5ad` | 76 MB | per-perturbation pseudobulk + their `obs.TE_ratio` (used for validation only) |
| `K562_essential_raw_singlecell_01.h5ad` | 10 GB | 310,385 cells with per-cell guide assignment in `obs.gene_transcript`, `obs.sgID_AB`, `obs.gene` |

The `_raw_singlecell_` h5ad has true raw UMI counts in `X` (despite float32 storage) but is not used for repeat quantification. Only its `obs` table is consumed: it provides the cell-barcode-to-perturbation mapping that joins to STARsolo output by `(gem_group, cell_barcode)`. Repeat quantification is run from scratch on fastqs downloaded from SRA (BioProject PRJNA831566).

Fastqs and BAMs are not on figshare. SRA holds fastqs only. The figshare manifest maps each `KD6_<gemgroup>_essential` library to its source fastqs by filename, and the BioProject lists the matching SRA accessions.

## Cell selection rules

Defaults in `paper/configs/repogle_2022.yaml`:

- `top_n_perturbations: 30`. top perturbations by Replogle's `TE_ratio`, excluding controls and perturbations with fewer than `min_cells_per_perturbation` cells. Top 30 covers exosome (EXOSC2-9, DIS3), CPSF + 3' end processing (SYMPK, CPSF1-4, CSTF3, WDR33, WDR82), Integrator (INTS2/4/5/8/9), Pol II + Mediator + PAF (POLR2 subunits, MED9/30, PAF1, CTR9, RUVBL1, SUPT5/6H), spliceosome (SNRPD1/2, SNRPF). Each major complex is represented by 5-10 subunits, giving within-complex consistency checks. The TE_ratio dynamic range is narrow (controls 0.013, max 0.035), so pushing past top 50 enters the noise floor.
- `n_controls: 1000`. random sample of non-targeting cells (~10% of available controls). Provides a stable baseline for pseudobulk DE while keeping cell-count balance comparable to Replogle's per-perturbation cell budget (~170 cells per top hit at this default).
- `min_cells_per_perturbation: 50`. filters perturbations with too few cells in their bulk file from the top-N ranking.
- `seed: 42`. reproducibility for the control sample.
- `gemgroup_subset: []`. empty means all 48 KD6 libraries. Set to a small list (e.g. `[1, 2, 3, 4]`) for a development run that touches only a few libraries' worth of fastqs.

These produce ~6,200 cells (5,200 perturbed + 1,000 control) spanning all 48 gemgroups when no subset is set.

## Pipeline phases

Phase 1 (implemented): figshare metadata retrieval and cell selection.

  - `paper/Snakefile_repogle_2022` rules `download_repogle_2022_figshare` and `select_repogle_2022_cells` download the three files, verify size and md5, and run `paper/scripts/select_repogle_2022_cells.py`.
  - Outputs three slim TSVs under `results/paper/repogle_2022/data/`:
    - `selected_perturbations.tsv` (gene_transcript, replogle_te_ratio, n_cells_replogle)
    - `cells_to_perturbation.tsv` (cell_barcode, gem_group, gene, sgID_AB, gene_transcript, umi_count)
    - `selected_libraries.tsv` (gem_group, library, BAM/matrix filenames, mRNA and sgRNA fastq lists)
  - The 10 GB singlecell h5ad is downloaded once per machine; the slim TSVs are committed and reused.

Phase 2 prep (implemented): SRA accession resolution per library.

  - Rule `resolve_repogle_2022_sra` runs `paper/scripts/resolve_repogle_2022_sra.py`, which queries SRA via `pysradb` for BioProject PRJNA831566 and matches each KD6 library's mRNA fastqs by `library_name` regex `KD6_seq\d+_essential_mRNA_lane_<gem_group>_`.
  - Outputs:
    - `selected_libraries_with_srrs.tsv` (the original libraries TSV plus `srrs` and `n_srrs` columns)
    - `repogle_kd6_sc_merged.yaml` (a single complete workflow config with `dataset.samples` injected, used as the only `--configfile` for Phase 2 alignment)
  - The sgRNA SRRs are deliberately excluded by the `mRNA` source pattern; per-cell guide assignments come from Replogle's published `obs.gene_transcript` already in `cells_to_perturbation.tsv`.

Phase 2 (alignment, run separately on the compute machine): SRA fastq download and STARsolo alignment.

  - Driven by `workflow/configs/repogle_kd6_sc.yaml` with `pipeline_type: sc`. Reference is `hg38` Ensembl 112 (shared with `gse230647_sc.yaml`).
  - The yaml encodes `sc_count_mode: tagcount` so the architecture is selected by config alone, no special CLI target. The default `all` rule honours `sc_count_mode`: in `tagcount` mode it requests only the canonical `unique_genes` BAM per library plus the `starsolo_tagcount/` outputs, dropping the per-feature_set fanout. To run the legacy per-feature_set path on Repogle, flip the yaml to `sc_count_mode: per_featureset` (no other change required; outputs land at `starsolo/...` and coexist with any `starsolo_tagcount/` outputs already on disk).
  - In `tagcount` mode: one STARsolo run per library aligned with `genes.gtf` for splice junctions and gene counting; the same BAM is then re-counted by `workflow/scripts/sc_count_features.py` against `genes.gtf`, `genic_repeats.gtf`, and `intergenic_repeats.gtf`. Compresses 132 STAR runs (44 libraries × 3 feature_sets) to 44 STAR runs plus ~132 cheap recount steps. Outputs land at `results/repogle_kd6_sc/starsolo_tagcount/<KD6_N_essential>/<mode>_<feature_set>/`.
  - In `per_featureset` mode: one STAR run per (library, feature_set) with `--soloMultiMappers Unique EM`. The unique and EM count matrices are written from the same alignment, and the EM matrix is mirrored into the `multi_<feature_set>/` sister directory so downstream rules see both modes without re-aligning. Outputs land at `results/repogle_kd6_sc/starsolo/<KD6_N_essential>/<mode>_<feature_set>/`.
  - Run: `cd workflow; snakemake --use-conda --cores N --configfile ../results/paper/repogle_2022/data/repogle_kd6_sc_merged.yaml`. The merged yaml inherits `sc_count_mode` from `repogle_kd6_sc.yaml` via `resolve_repogle_2022_sra`.
  - Per-library STARsolo matrices land at gene_id, family_id, and class_id granularity in either output tree.

Phase 3 (per-perturbation RUVg DE plus per-library alignment QC; needs Phase 2 outputs):

  - `aggregate_star_alignment_qc` parses each library's `unique_genes/Log.final.out` via `workflow/scripts/parse_star_log.py` and writes `report/alignment_qc_per_library.tsv` (mismatch rate, mapping rate, unmapped fractions, average mapped length).
  - `render_repogle_2022_perturbseq_report` runs `paper/scripts/ruv_repogle_2022_perturbseq_report.Rmd`.
  - The Rmd loads STARsolo matrices (multi mode reads `UniqueAndMult-EM.mtx`), joins cells to perturbations via `cells_to_perturbation.tsv`, and builds per-(perturbation, gem_group) pseudobulks. Each perturbation gets ~30-48 columns (one per gemgroup with cells of it), matched by ~48 control columns from the same gemgroups. Gemgroup is the blocking factor. Pseudobulks stay sparse until they enter DGEList.
  - RUVg fit on the gene pseudobulk via `workflow/scripts/ruv_common.R`; the same W is reused across repeat featuresets per mode.
  - Per-perturbation DE under `~ gem_group + W + group`, run at gene_id, family_id, and class_id. Family and class matrices come from `workflow/scripts/granularity_rollup.R` plus the locus_map TSVs at `<indices_base>/indices/all/<feature_set>_locus_map.tsv`.
  - Alignment-QC scatter: one strip per perturbation, one point per contributing library (size = cell count), median highlighted. Library mismatch rates come from STAR `Log.final.out`. Used as a batch-confound check: a perturbation whose median sits high is enriched in poorly-aligning libraries.
  - Family/class summary: per-perturbation logFC heatmap at class_id and logFC distribution at family_id, faceted by polyA-competence (`POLYA_CLASSES` from `workflow/scripts/biology_checks.R`). Class for family_id rows comes from the locus_map; class for gene_id rows comes from `repeat_biology_annotation.tsv`.
  - Validation: REclaim's per-perturbation `te_total_fraction` against Replogle's `te_ratio`, with Pearson and Spearman per mode.

## Why pseudobulk

Replogle's `TE_ratio` is itself per-perturbation pseudobulk, so a one-to-one validation requires the same level of aggregation. Per-cell TE quantification at 10x 3' depth (~10k UMIs per cell) puts most subfamily-level signal below per-cell noise. Pseudobulk lifts SNR and mirrors the existing TDP-43 paper-side workflow (which pseudobulks per cluster), so the same edgeR plus RUVg machinery applies with only the grouping key changed. Single-cell-resolution analysis is a possible follow-up to ask whether perturbation effects are uniform across cells or driven by subsets with stronger knockdown.

## Why both unique and multi

STARsolo is run with `--soloMultiMappers Unique` and `--soloMultiMappers EM` as separate count streams.

- Unique counts skew toward older, uniquely-mappable repeats (MIR, divergent Alu, ancient LINEs) and discard most reads from young elements (L1HS, AluY, SVA) where multi-mapping dominates.
- Multi/EM redistributes multi-mappers across mapping locations and recovers young-element signal at the cost of slightly noisier per-locus assignments.
- Cross-mode agreement on per-perturbation TE response is a robustness check. Mode-specific divergence (e.g. EXOSC knockdown shifts multi-counts more than unique-counts) localizes the response to particular subfamily age strata, which is a story Replogle's class-level metric cannot tell.

The simulation benchmark (`workflow/Snakefile` with `pipeline_type: simulation`) provides an empirical ranking of unique vs multi accuracy on synthetic ground truth. If those numbers strongly favor one mode, the operating mode for downstream figures can be picked accordingly; the other mode stays in the supplement.

## Reproducibility

- Figshare URLs and md5s are pinned in the config.
- Cell selection is seeded.
- The slim TSVs are intended for commit; the h5ad is not.
- Versioning: the `_01` suffix on Replogle's filenames is a chunk index. As of figshare deposit version 2 (2022), only `_01` exists for K562 essential. If they release a re-versioned deposit, the URLs and md5s change and need updating in the config.

## Run

Phase 1 plus Phase 2 prep (figshare downloads, cell selection, SRA accession resolution):

```
source ~/miniconda3/bin/activate
conda activate snakemake
cd paper
snakemake -s Snakefile_repogle_2022 --use-conda --cores 2 \
    --configfile configs/repogle_2022.yaml
```

Phase 2 alignment (run on the compute machine after Phase 1+prep):

```
cd workflow
snakemake --use-conda --cores N \
    --configfile ../results/paper/repogle_2022/data/repogle_kd6_sc_merged.yaml
```

Phase 3 (per-perturbation pseudobulk RUVg DE, after Phase 2 alignment outputs are in place):

```
cd paper
snakemake -s Snakefile_repogle_2022 --use-conda --cores 4 \
    --configfile configs/repogle_2022.yaml \
    -- ../results/paper/repogle_2022/report/ruv_repogle_2022_perturbseq_report.html
```

Outputs land under `../results/paper/repogle_2022/`: `data/` for the slim TSVs and `repogle_kd6_sc_merged.yaml`, `figshare/` for the (uncommitted) h5ad and manifest, `report/` for the Phase 3 HTML and `ruv_rds/` snapshots, and `logs/` for per-rule logs.
