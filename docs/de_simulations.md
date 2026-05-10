# de_simulations

Count-level benchmark for repeat differential expression. Two scenarios share the same code path.

## Two scenarios, two questions

The default scenario (`configs/de_simulations_polymenidou_bulk.yaml`) asks a power question: if a small fraction of repeats are DE, how often does each normalization call them at FDR <= 0.05, and how does that change with depth.

The coordinated-derepression scenario (`configs/de_simulations_coordinated_derepression.yaml`) asks a bias question: if a large fraction of repeats all move in the same direction, does TMM-on-repeats absorb the coordinated signal and compress the recovered logFC. This mirrors the TDP-43 overexpression contrast where many young retrotransposon subfamilies move up together. Gene-library-size transfer (TMM size factors fit on the gene matrix, then transferred to the repeat DGEList) is the comparator: the gene matrix still satisfies the "most features not DE" assumption that TMM relies on, so size factors derived there are unbiased.

## How

1. Load gene and repeat count TSVs from a real run.
2. Fit per-feature NB mean and dispersion with edgeR.
3. Choose a planted-DE set on the repeat matrix:
   - `n_de_repeats` features picked,
   - `signed_fraction` of them positive logFC (rest negative),
   - optional class-weighted selection via `class_map_path` and `class_weights`.
4. For each (gene_lib_scale, repeat_lib_scale) on the grid, simulate `n_iter` count matrices: scale the mean by the cell's scaler, add a shared latent W loading per feature, plant the chosen logFC.
5. Run no-norm, TMM, gene_lib (gene-derived size factors), and RUVg(k=1..k_max) on each replicate. RUVg picks empirical controls from the gene matrix.
6. Score each call set at FDR <= `fdr` against the planted truth and record per-feature recovered logFC for bias-slope estimation.

Planted features that fail the count filter still count as false negatives.

## Methods covered

- `none`: edgeR with `norm.factors = 1`. Lower bound.
- `TMM`: edgeR's TMM size factors fit on the repeat matrix. Standard practice.
- `gene_lib`: TMM size factors fit on the gene matrix and transferred to the repeat DGEList. Sidesteps the "most features not DE" assumption when many repeats co-move.
- `RUVg_k1` .. `RUVg_kK`: RUVSeq with empirical controls picked from the gene matrix as the least-significant genes under a naive F-test on TMM-normalized gene counts.

## Output

- `de_simulations_results.tsv`: one row per (grid cell, replicate, method) with power, FPR, TP/FN/FP/TN.
- `de_simulations_summary.tsv`: per-cell mean and SD across replicates.
- `de_simulations_recovered_logfc.tsv.gz`: one row per (grid cell, replicate, method, feature) with `planted` flag, `planted_logfc`, `recovered_logfc`, `class_id`. Drives the recovered-vs-input bias panel in the report.
- `de_simulations_results.rds`: full results plus fitted NB params, recovered-logfc table, and planted indices.
- `de_simulations_heatmap.pdf`: power and FPR heatmaps faceted by method.
- `de_simulations_report.html`: rendered report with heatmaps, precision, F1, per-replicate spread, NB-params sanity check, and the recovered-logFC bias scatter when the coordinated-derepression flag is active. The Rmd source is `workflow/scripts/de_simulations_report.Rmd`.

## How to read the heatmap

Rows are repeat library scalers, columns are gene library scalers, color is the metric (0 to 1). Bright cells mean the method recovers the planted FC at that depth combination. Dark bottom rows mean repeat depth is the bottleneck. Dark left columns mean gene depth starves RUVg's W estimation.

## How to read the recovered-vs-input scatter

X axis: planted logFC. Y axis: recovered logFC. One panel per method. Under the default scenario, all methods sit on the y = x line (no bias). Under coordinated derepression, TMM on the repeat matrix sits on a line with slope < 1 (compression), gene_lib stays close to y = x (unbiased), and RUVg sits in between depending on whether W absorbed the coordinated signal. The slope is the headline number.

## Inputs

- `gene_counts`: TSV, feature_id rownames, sample columns.
- `repeat_counts`: TSV, same sample columns. Granularity (locus, family, class) is whatever the input rows are.
- `metadata` (optional): TSV with `sample` and `condition`. Empty splits samples in halves (A then B).
- `class_map_path` (optional): either a 2-column `feature_id, class_id` TSV, or the existing 4-column locus_map (`transcript_id, gene_id, family_id, class_id`) produced by `reference.snmk::build_repeat_locus_map`. When passing a locus_map, set `class_map_granularity` to the column whose values match the repeat-matrix rownames.

The two count TSVs come from a real bulk run. This pipeline does not align reads.

## Config keys

Under `config['de_simulation']`:

- `gene_counts`, `repeat_counts`: required input paths.
- `metadata`, `sample_column`, `condition_column`: optional sample table.
- `fc`: planted fold change (default 3).
- `n_iter`: replicates per cell (default 20).
- `n_de_repeats`, `n_de_genes`: planted DE counts.
- `gene_lib_grid`, `repeat_lib_grid`: scaler vectors.
- `sigma_w`: SD of per-feature W loading.
- `k_max`: max RUVg factors.
- `n_controls`: empirical controls per replicate.
- `fdr`: FDR threshold (default 0.05).
- `seed`: master seed.
- `signed_fraction` (default 1.0): fraction of `n_de_repeats` planted with positive logFC; 1.0 = all up (coordinated derepression), 0.5 = balanced, 0.0 = all down.
- `class_map_path` (default ""): TSV mapping feature_id (or locus_map column) to class_id.
- `class_map_granularity` (default `family_id`): which column of the locus_map to use as feature key. Ignored when `class_map_path` is a 2-column file.
- `class_weights` (default {}): named map from class_id to selection weight; requires `class_map_path`.

## Run

Via the Makefile from the project root:

```
make de_polymenidou_bulk CORES=N             # default scenario
make de_coordinated_derepression CORES=N     # coordinated derepression scenario
```

Each invocation runs the main `workflow/Snakefile` with `pipeline_type: de_simulation` (declared inside the configfile), runs the `de_simulations` rule, and renders the HTML report in the same DAG.

Or directly:

```
snakemake --use-conda --cores N --configfile configs/de_simulations_polymenidou_bulk.yaml
snakemake --use-conda --cores N --configfile configs/de_simulations_coordinated_derepression.yaml
```

Outputs go to `{base}/de_simulations/`. The report is rendered by `rule render_de_simulations_report` inside `workflow/modules/de_simulations.snmk` using the `rmarkdown` conda env, and depends on the artifacts above so it is rebuilt whenever the simulation outputs change.

## Scope

This is count-level only. It does not test alignment or multi-mapper assignment. For those, see `workflow/modules/simulations.snmk` and the noise sweep configs.
