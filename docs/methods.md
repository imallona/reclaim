# Methods

## Overview

A single Snakefile drives five pipeline modes selected by `pipeline_type` in the config.

The simulation pipeline generates reads from repeat loci (SmartSeq2 or 10x Chromium), aligns them with four tools, and evaluates quantification accuracy against ground truth.

The bulk pipeline downloads FASTQs from SRA via fasterq-dump. Reads are aligned with STAR, kallisto, and salmon. edgeR produces differential expression reports at gene_id and family_id granularity. Read layout is controlled by `real_data.library_layout`: `paired` (R1 + R2) or `single` (R1 only). Paired-end runs use `bulk_paired.snmk`. Single-end runs use `bulk_single.snmk`, which passes `--single --fragment-length --sd` to kallisto and `--unmatedReads` to salmon.

The sc pipeline downloads 10x Chromium FASTQs from SRA with `fasterq-dump --include-technical`. File numbering follows fasterq-dump conventions. File `_1` is the index read and is discarded. The CB+UMI and cDNA reads are selected via `sc_cbumi_index` and `sc_cdna_index` (defaults 2 and 3, matching 10x Chromium v2/v3). STARsolo (CB_UMI_Simple) and kallisto|bustools quantify genic and intergenic repeat annotations.

The de_simulation pipeline is a count-level power benchmark. It loads pre-computed gene and repeat count TSVs from a real run, fits NB parameters, and resamples count matrices over a 2D grid of gene and repeat library size scalers, planting a known fold change at a fixed feature set. It scores no-norm, TMM, and RUVg(k=1..k_max) for power and FPR. See [docs/de_simulations.md](de_simulations.md) for the full description.

## Reference preparation

Genome sequence (GRCh38), gene annotation (Ensembl GTF), and repeat annotation (RepeatMasker, UCSC RMSK format converted to GTF via makeTEgtf) are downloaded or provided by the user. All shared reference files and aligner indices sit under a single `indices_base` directory. Multiple runs with different simulation parameters (e.g. noise sweeps) reuse the same indices without rebuilding.

Internal files use bare chromosome names without the `chr` prefix (Ensembl convention).  UCSC-sourced files with `chr`-prefixed chromosome names are stripped at decompression time.

When `filter_genic` is enabled, three annotation subsets are produced using bedtools intersect:

- `genic_repeats`: repeat loci overlapping at least one gene body
- `intergenic_repeats`: repeat loci with no overlap to any gene body
- `repeat_harboring_genes`: gene records overlapping at least one repeat locus

These subsets are saved as GTF files for downstream filtering and stratified evaluation.

## Simulation

Reads are generated directly from repeat element loci defined in a RepeatMasker GTF annotation. The genome FASTA is streamed one chromosome at a time. Only chromosomes with at least one sampled locus are loaded into memory, which avoids holding the full genome in RAM.

For each simulated cell a random subset of repeat loci is drawn without replacement. The number of expressed loci per cell varies around a user-specified mean (`n_expressed_per_cell`, default 1000). The per-cell count is sampled uniformly from [0.7x, 1.3x] of that mean. The count assigned to each expressed locus comes from a geometric distribution with mean 5, truncated at 50. The geometric distribution gives the heavy-tailed, low-count profile typical of scRNA-seq without external dependencies.

Each read is a subsequence sampled uniformly at random from the repeat locus sequence. Substitution errors are introduced independently at each base with probability `mutation_rate` (default 0.001). A noise sweep runs four mutation rates (0%, 1%, 5%, 10%) to test how aligners hold up under sequencing errors.

In SmartSeq2 mode each cell is written to a separate gzipped FASTQ file alongside a STARsolo-compatible manifest.

In Chromium mode cell barcodes are drawn at random (16 bp, matching 10x v3). A single paired-end FASTQ pair is generated: read 1 carries the cell barcode and UMI (12 bp), read 2 carries the cDNA.

Ground truth is a long-format TSV with one row per expressed (cell, repeat) pair. `family_id` and `class_id` columns come from the GTF. Only loci whose genome sequence is usable (< 10% N bases, >= 50 bp) appear in the ground truth.

## Quantification granularity

All repeat-aware aligners produce count matrices at configurable aggregation levels:

| Level     | RepeatMasker field | Example value |
|-----------|--------------------|---------------|
| locus     | transcript_id      | AluSz6_dup1   |
| gene_id   | gene_id            | AluSz6        |
| family_id | family_id          | Alu           |
| class_id  | class_id           | SINE          |

The `granularities` config key specifies which levels to compute.  Ground truth is aggregated to the same level before evaluation so metrics are always computed at matching resolution.

## Aligners and quantifiers

### STARsolo

STAR 2.7+ runs in SmartSeq (`soloType SmartSeq`) or Chromium (`soloType CB_UMI_Simple`) mode. A custom `sjdbGTFfile` selects the repeat or gene annotation. Multimapper handling is controlled by `--soloMultiMappers`: `Unique` (only uniquely mapping reads contribute) or `EM` (multimapping reads distributed by the EM algorithm). Both modes are evaluated for all technologies.

For sc Chromium specifically, the alignment+counting architecture is configurable via the top-level config flag `sc_count_mode`. The two modes differ in how the per-feature_set count matrices are produced from the same fastqs. The yaml is the single source of truth: the default `all` rule reads `sc_count_mode` and routes to the right outputs, so no special target on the snakemake CLI is required.

`sc_count_mode: per_featureset` (default; backward-compatible). One STARsolo run per `(sample, feature_set)`; `--sjdbGTFfile` is the feature_set's own GTF; `--soloFeatures Gene` emits the count matrix natively. Multi mode uses STARsolo's iterative `EM` redistribution. Outputs land at `<base>/starsolo/<sample>/<mode>_<feature_set>/`.

`sc_count_mode: tagcount`. One STARsolo run per library, aligned with `genes.gtf` for splice junctions and gene counting. The resulting BAM (with `CB` and `UB` tags for Chromium, `CB` only for SmartSeq2) is then re-counted by `workflow/scripts/sc_count_features.py` against each repeat feature_set GTF. This avoids three independent STAR alignments per library, and splice junctions are inferred from the gene model only (not from repeat GTFs, which would otherwise inject spurious junctions). The recount is constrained to repeat feature_sets (`repeats`, `genic_repeats`, `intergenic_repeats`); at `(genes, gene_id)` STARsolo's own `Solo.out/Gene/{filtered,raw}/matrix.mtx` is canonical and re-counting there would duplicate work. Outputs land at `<base>/starsolo_tagcount/<sample>/<mode>_<feature_set>/` so they coexist on disk with `per_featureset` results from sibling yamls and don't trigger downstream rebuilds.

The same recount path serves both Chromium and SmartSeq2 simulations through two technology-specific modules. `starsolo_tagcount.snmk` handles Chromium: the BAM carries real UMIs in `UB`, the script runs with `--umi-dedup 1mm_all` (default) or `exact`, and the whitelist is the EmptyDrops_CR-passed `Solo.out/Gene/filtered/barcodes.tsv` for real-data sc, or `raw/barcodes.tsv` for simulations (which run `--soloCellFilter None`). `smartseq2_tagcount.snmk` handles SmartSeq2: STARsolo SmartSeq2 emits one BAM per library where each read carries `CB` from the manifest cell-id (no `UB`, no real UMIs); the script runs with `--umi-dedup none`, the whitelist is `raw/barcodes.tsv` (= manifest cell list). The shared BAI rule lives in `bam_index.snmk` so both technologies depend on the same indexer.

The recount script `sc_count_features.py` is the load-bearing piece in `tagcount` mode. Its design choices, exposed as CLI flags and configurable from the workflow config:

- Indexing level. The GTF is parsed at the broadest available feature_type (gene > transcript > exon). Gene-level indexing means a read sitting in an intron still counts toward its parent gene, mirroring STARsolo's `GeneFull` semantics. Repeat GTFs use exon-only entries (one exon per repeat locus) so the entire repeat span is covered by definition.
- Cell barcode whitelist. Restricting to the CBs that STARsolo's gene-counting pass already retained avoids the multi-million-row unfiltered all-CB matrix and aligns the repeat counts with the same set of "real cells" used for genes.
- UMI deduplication (`--umi-dedup`). `1mm_all` is the default and matches STARsolo's `--soloUMIdedup 1MM_All`: per `(CB, gene_id)`, UMIs within Hamming-1 of each other collapse via union-find before counting. `exact` counts distinct UB strings without collapsing and is for debugging only and systematically over-counts unique molecules by ~2-5% because sequencing errors create spurious UMIs. `none` skips UMI bucketing entirely: each BAM record is keyed by a per-record synthetic UB (`(qname, ref_start)`) so every record contributes its weight (1.0 in unique mode, 1/NH in multi mode). This is the SmartSeq2 path: with no real UMIs, no honest collapsing is possible, so every BAM record contributes one count. For paired-end SmartSeq2 BAMs, only R1 is counted (the script skips records with `is_paired` and `is_read2` set) so a fragment is not double-counted across its two ends. Without UMIs the pipeline does not attempt any duplicate removal: in RNA-seq it is impossible to tell technical PCR duplicates from coincident, biologically real reads at the same coordinates without a molecular barcode.
- Multimapper handling (`--multimapper`). `unique` skips `NH > 1` reads and is identical to STARsolo's `Unique`. `multi` includes `NH > 1` reads with weight `1/NH` per alignment record. This is non-iterative and linear-cost, unlike STARsolo's `--soloMultiMappers EM` which iterates a redistribution step. The two definitions differ systematically: EM concentrates weight on genes with strong unique-mapper support, `1/NH` spreads it uniformly across the candidate loci. For repeat-family quantification the uniform spread is the more conservative choice and is widely used in tools like RSEM (no-EM) and scTE.
- Granularity outputs. When the GTF carries `family_id` and/or `class_id` attributes (the repeat GTFs do; Ensembl genes do not), the script emits per-granularity matrices alongside the gene_id matrix, at `Solo.out/Gene_family_id/raw/` and `Solo.out/Gene_class_id/raw/`. The rollup is post-UMI-dedup at the gene_id level.
- Parallelism. One worker per BAM contig; each worker holds its own pysam handle (handles are not fork-safe) and processes one chromosome. Memory is dominated by the per-key UMI dedup state and stays under 1 GB for typical 10x libraries.

Unit tests in `test/unit/test_sc_count_features.py` cover, against synthetic GTFs and BAMs: exonic vs intronic counting, multi-feature ambiguity, missing CB/UB tags, off-whitelist barcodes, multimapper handling, UMI dedup at exact, Hamming-1, and `none` mode (SmartSeq2 path: counts every record, distributes multi-mappers across loci, still applies the whitelist), and family/class rollups. In unique mode with `1mm_all` dedup, the recount is expected to agree exactly with STARsolo's gene matrix produced from the same BAM; the comparison Rmd asserts this as a regression check before the project relies on the recount in production.

### Kallisto

For SmartSeq2 (simulation) and bulk single-end, `kallisto quant` runs per sample against a repeat-sequence pseudo-transcriptome. Single-end bulk passes `--single --fragment-length --sd`. Fragment length and standard deviation are set via `real_data.fragment_length` and `real_data.fragment_sd` (defaults 200 and 20). Per-sample `abundance.tsv` files are merged and aggregated at the requested granularity.

For bulk paired-end, `kallisto quant` is run without `--single`, passing both R1 and R2.

For 10x Chromium, `kallisto bus` runs on the combined R1/R2 files. The resulting BUS file is sorted and corrected with bustools. Per-cell count matrices are produced with `bustools count`.

### Alevin (Salmon)

For SmartSeq2 and bulk paired-end, `salmon quant` runs per sample in quasi-mapping mode. `-1` and `-2` pass the two mates. `--minAssignedFrags 1` stops Salmon from exiting when few reads map.

For bulk single-end, `salmon quant` is run with `--unmatedReads` instead of `-1`/`-2`.

For 10x Chromium, `salmon alevin` is run in Chromium v3 mode using the repeat pseudo-transcriptome index.  Output is converted by `normalize_alevin_chromium.py`.

### Bowtie2: pseudo-genome approach

Bowtie2 is indexed against a pseudo-genome FASTA where each repeat locus is its own reference sequence, named `transcript_id::chrom:start-end(strand)`. Reads align with `bowtie2 -a` (all alignments). Per-locus counts come from `samtools idxstats`, which counts alignments per reference without coordinate lookup.

This pseudo-genome layout is the reason bowtie2 cannot be migrated onto the `sc_count_features.py` tagcount path that Chromium and SmartSeq2 STAR alignments use. Tagcount maps an alignment's `(chrom, ref_start)` to a feature interval in a real-genome GTF; a bowtie2 BAM has its references in pseudo-genome space (one contig per repeat locus), so a real-genome GTF lookup is meaningless. Bowtie2 keeps its own counting path: `samtools idxstats` for the repeat granularities and featureCounts (cell-chunked, GTF-chunked) for the genes feature_set in SmartSeq2.

For SmartSeq2, one BAM is produced per cell; featureCounts is then run on configurable cell chunks to limit peak memory.

For 10x Chromium, read 2 (cDNA) aligns to the repeat index. Cell barcode and UMI tags from read 1 are attached via `umi_tools extract`. UMI deduplication runs via `umi_tools dedup` before counting.

## Normalization

Each aligner's native output is converted to a common feature x cell TSV (rows = repeat features, columns = cell identifiers) by a per-aligner `normalize_*.py` script.  Only non-zero rows are written.  This format is the sole input to the evaluation step.

For 10x Chromium, kallisto and alevin normalization uses sparse accumulators. Only non-zero `{cell_index: count}` pairs are stored per feature group. Memory stays proportional to expressed pairs rather than O(features x cells).

The bowtie2 Chromium counting script streams the CB-tagged deduplicated BAM once via `samtools view`, accumulating `(barcode, feature)` counts in a sparse dict without per-cell BAM splitting.

### Repeat differential expression: TMM, gene_lib, and RUVg

Differential expression on the repeat count matrix is run with edgeR's quasi-likelihood (`glmQLFit` + `glmQLFTest`). Three normalization paths are available and reported side by side in the count-level benchmark and in the paper-side reports:

- **TMM on repeats**: `calcNormFactors(method = "TMM")` fit on the repeat matrix. The standard practice. Assumes most repeat subfamilies are not differentially expressed; when this assumption holds, TMM is unbiased.
- **gene_lib transfer**: `calcNormFactors(method = "TMM")` fit on the *gene* matrix, then transferred to the repeat DGEList by setting `lib.size` and `norm.factors` from the gene fit. The repeat matrix never sees its own size factors. Used when many repeats co-move and the "most features not DE" assumption is suspect (the gene matrix still satisfies it). See `workflow/scripts/gene_lib_size_common.R::run_gene_lib_repeat_de`.
- **RUVg(k)**: empirical controls picked from the gene matrix as the least-significant genes under a naive F-test on TMM-normalized gene counts; `RUVg(set, controls, k)` fits k unwanted-variation factors W; the repeat DE design becomes `~ W + group` (or `~ sample_block + W + group`). See `workflow/scripts/ruv_common.R`.

The count-level benchmark module exercises all three (plus a `none` path with `norm.factors = 1` as a lower bound) and exposes recovered logFC per planted feature so a bias slope can be computed per method.

### Coordinated-derepression scenario

The default count-level scenario (`workflow/configs/de_simulations_polymenidou_bulk.yaml`) plants a small, mixed-direction set of DE repeats: TMM is in its happy regime and the metric of interest is power across the depth grid. The coordinated-derepression scenario (`workflow/configs/de_simulations_coordinated_derepression.yaml`) plants a large, same-sign set of DE repeats with class-weighted selection biased toward young polyA-competent retrotransposons (LINE, SINE, LTR, Retroposon). It mirrors the TDP-43 overexpression biology and is the controlled test of the TMM-compression claim. Recovered logFC is exported per (method, feature, replicate, grid cell) to `de_simulations_recovered_logfc.tsv.gz` and the report Rmd computes per-method bias slopes (recovered ~ planted) at the canonical depth cell. See `docs/de_simulations.md` for the full parameter set.

## Evaluation

`evaluate.py` compares each aligner's count matrix against the simulation ground truth at three levels:

Global (one row per aligner x granularity): Pearson r, Spearman r, log1p RMSE, precision, recall, F1, Jaccard index, specificity.

Per cell: Pearson r and Spearman r computed separately for each cell across all repeat loci.

Per repeat class: All global metrics stratified by `class_id` (LINE, SINE, LTR, DNA transposon, satellite, etc.), revealing which repeat types are accurately quantified.

Compute resource metrics (wall time, CPU time, peak RSS, I/O) are read from Snakemake benchmark files and included in the summary table.

The `aggregate_global_metrics` rule concatenates all per-aligner global metric TSVs into a single `summary_global_metrics.tsv`.

An HTML report (`evaluation_report.html`) is rendered by `evaluation_report.Rmd` using `rmarkdown::render` with ggplot2 + patchwork. The report also writes paper-ready CSV tables to `{eval_dir}/paper_csv/`: `global_metrics_long.csv` (tidy long format), `global_metrics_wide.csv`, `aligner_ranking_per_metric.csv` (rank per feature_set x granularity x metric, direction-aware), `best_aligner_per_metric.csv` (best, second, delta), `per_cell_metrics_summary.csv` (n, mean, median, sd, q25, q75 per aligner x granularity x feature_set), `per_class_metrics.csv`, and `resources_wide.csv`.

### External tool benchmark

When `external_benchmark: true` is set in a simulation config, the module `workflow/modules/external_tools_benchmark.snmk` runs published repeat-quantification tools on the same simulation BAMs and ground truth as REclaim, scored through the same `evaluate.py`. Tool coverage:

| Tool | Modality | Native granularity | Conda env |
|---|---|---|---|
| TEtranscripts (TEcount) | bulk | gene_id (subfamily) | `workflow/envs/tetranscripts.yaml` (TEtranscripts 2.2.3) |
| scTE | single cell | family_id | `workflow/envs/scte.yaml` (scTE 1.0.0) |
| SQuIRE | bulk | locus | `workflow/envs/squire.yaml` (SQuIRE 0.9.9.92, Python 3.6, STAR 2.5.3a; isolated env) |

REclaim is the only quantifier here that natively emits at every granularity (locus, gene_id, family_id, class_id) from a single BAM pass through the GTF feature hierarchy. External tools emit at one native level. Rows scored at non-native granularities are rolled up via the same locus_map REclaim uses, so a comparison at those levels measures the native counting step rather than the rollup logic.

Per tool, the rule chain is:

1. Run rules invoke each upstream tool against the existing STARsolo BAMs. TEcount runs once per simulated SmartSeq2 cell. scTE runs once on the multi-cell Chromium BAM. SQuIRE runs as a chain (`Fetch` -> `Clean` -> per-sample `Count` -> combine); the existing STARsolo BAMs are symlinked into a SQuIRE-shaped `map_folder` instead of running a second alignment. Outputs land under `{base}/external_benchmark/{tool}/`.
2. `harmonize_external` calls `workflow/scripts/harmonize_external_counts.py` to convert the tool's native output to the REclaim feature x cell TSV format. The harmoniser asserts that at least `external_min_overlap_fraction` of the tool's feature IDs map to a known REclaim feature in the locus_map. For scTE on Chromium, the harmoniser also takes a `barcode_to_cell_id.tsv` map produced by the simulator so the h5ad `obs_names` (cell barcodes) are translated into REclaim cell_ids before sample filtering.
3. `score_external_benchmark` runs the same `evaluate.py` used for REclaim's own quantifiers, with the same locus_map and granularity contracts. Wall time and peak RSS are aggregated across per-sample runs for TEcount and SQuIRE so the cost panel in the report carries one value per tool rather than one per sample.

The render rule reads `summary_global_metrics.tsv` from both `external_benchmark/` and `evaluation/` and produces side-by-side accuracy bars per quantifier, granularity, and class. REdiscoverTE is not included; rationale in `docs/external_tool_benchmark.md`.

Sample scope. The benchmark scores every simulation cell (`simulation.n_cells`) by default. The snmk module derives the cell list from the simulation config when `external_samples` is omitted, so the external benchmark and the REclaim evaluation cover the same cells and the per-quantifier accuracy bars use the same sample set.

A separate noise sweep report (`noise_sweep_report.Rmd`) loads `summary_global_metrics.tsv` from each noise-level run and plots metric degradation as a function of mutation rate. It writes paper-ready CSV tables to `{csv_outdir}/paper_csv/` (by default the directory of the HTML output): `noise_metrics_long.csv`, `noise_metrics_wide.csv` (columns per mutation rate), `noise_degradation_slopes.csv` (linear fit `value ~ mutation_rate` with slope, intercept, R^2, n_points per aligner x feature_set x granularity x metric), `noise_robustness_ranking.csv` (ranks aligners by shallowest slope for accuracy metrics), `per_cell_noise_summary.csv`, and `per_class_noise_metrics.csv`.

## Conda environments

| Environment | File                 | Key packages                              |
|-------------|----------------------|-------------------------------------------|
| star        | envs/star.yaml       | STAR, pigz, bedtools, gffread, samtools   |
| kallisto    | envs/kallisto.yaml   | kallisto, bustools                        |
| alevin      | envs/alevin.yaml     | salmon                                    |
| bowtie2     | envs/bowtie2.yaml    | bowtie2, samtools                         |
| umi_tools   | envs/umi_tools.yaml  | umi_tools, samtools (>= 1.12), pysam       |
| evaluation  | envs/evaluation.yaml | python, scipy                             |
| rmarkdown   | envs/rmarkdown.yaml  | R, rmarkdown, ggplot2, patchwork          |
| edger       | envs/edger.yaml      | R, edgeR, rmarkdown, ggplot2              |
| ruvseq      | paper/envs/ruvseq.yaml | R, edgeR, RUVSeq, EDASeq, rmarkdown     |
| tetranscripts | envs/tetranscripts.yaml | TEtranscripts 2.2.3, samtools         |
| scte        | envs/scte.yaml       | scTE 1.0.0, anndata, h5py                 |
| squire      | envs/squire.yaml     | SQuIRE 0.9.9.92 (Python 3.6, STAR 2.5.3a) |

## Paper analyses

The `paper/` directory holds a separate Snakefile for dataset-specific reanalyses used in the manuscript. It is kept out of `workflow/Snakefile` for two reasons. First, paper analyses are cross-dataset, with a single config referencing outputs of multiple method runs. Second, paper rules should not retrigger method jobs via Snakemake provenance tracking.

Entry point: `paper/Snakefile`, driven by `paper/configs/tdp43.yaml`. Counts are consumed from the method-pipeline output directory declared in each dataset entry (`counts_dir` for bulk, `starsolo_base` for sc), so the relevant method pipeline must be run first for each dataset. Dataset keys use last author plus study label plus assay kind; each entry also records the source GEO accession (`gse`) and the method-pipeline config (`source_config`) that produced its inputs.

Paper datasets currently wired in `paper/configs/tdp43.yaml`:

| Paper dataset key | GSE | Assay | Source method config |
|---|---|---|---|
| `lee_nuclear_tdp_retention_bulk` | GSE126543 | bulk | `workflow/configs/gse126543_bulk.yaml` |
| `polymenidou_tdp_kd_oe_cultures_bulk` | GSE230647 | bulk | `workflow/configs/gse230647_bulk.yaml` |
| `polymenidou_tdp_oe_cultures_scrnaseq` | GSE230647 | scRNA-seq | `workflow/configs/gse230647_sc.yaml` |

### TDP-43 bulk differential expression

The `render_ruv_bulk` rule runs a per-dataset RUVg-normalized edgeR analysis on STAR unique-count matrices. A dataset wildcard lets `polymenidou_tdp_kd_oe_cultures_bulk` and `lee_nuclear_tdp_retention_bulk` share a single rule definition. Inputs are all repeats, genic repeats, and intergenic repeats at gene_id and family_id granularity, plus the gene matrix. STAR unique counts outperform kallisto and salmon on repeats in the simulation benchmark. The kallisto and salmon matrices are not consumed by the paper rule. The baseline bulk edgeR reports are kept; the paper rule adds a second, RUVg-normalized analysis rather than replacing them.

Per-dataset procedure, encoded in `paper/scripts/ruv_{dataset}_report.Rmd`:

1. Load the gene matrix and the six repeat matrices. Filter the gene matrix to rows with count >= 10 in at least half the samples; repeat matrices are not filtered at this stage.
2. Fit a naive edgeR quasi-likelihood model on the filtered gene matrix with TMM normalization. The design is two-group for GSE126543 and four-group for GSE230647, with F-test across the three non-reference coefficients. Rank genes by p-value. Take the top `n_controls` with the largest p-values (default 5000) as empirical controls. These are the genes least associated with the contrasts of interest.
3. Run RUVg with k = 1..4 on the gene matrix. Score each k by the mean-absolute-deviation of the relative log expression distribution. The full k = 0..4 table is written to `ruv_rle_diagnostic.tsv` for audit. The operating k is set via `k_selected` (default 2), not picked automatically from the RLE. RUVg produces per-sample W factors that do not depend on feature granularity or featureset, so the same W is reused across all repeat contrasts.
4. For each repeat featureset and contrast: build a DGEList on the raw repeat counts, compute TMM factors, and fit `~ W + group` with `glmQLFit` and `glmQLFTest`. The ranked TSV contains feature_id, logFC, logCPM, F, PValue, and FDR. GSE126543 drops the two unsorted samples and runs a single tdp43pos vs tdp43neg contrast. GSE230647 runs two contrasts: TDP-43 knockdown vs non-targeting and TDP-43 overexpression vs DOX-off control.

Outputs per dataset are written under `results/paper/tdp43/{dataset}/`. The report lives alongside one ranked DE TSV per featureset x contrast named `ruv_k{k}_{featureset}_{label}.tsv`. RDS snapshots (`ruv_state.rds`, `ruv_naive_gene_de.rds`, `ruv_repeat_de.rds`) are saved next to `ruv_empirical_controls.txt`, `ruv_rle_diagnostic.tsv`, and `ruv_summary.tsv`.

Shared helper functions used by both bulk Rmds live in `workflow/scripts/ruv_common.R`: `load_counts_tsv`, `rle_score`, `rank_empirical_controls`, `fit_ruv_range`, `pick_W_matrix`, and `run_ruv_repeat_de`. The same helpers are reused by the single-cell pseudo-bulk RUVg analysis so that control selection, RUVg fitting, and edgeR contrast code are identical across bulk and sc.

A parallel gene-library-size repeat DE is run alongside RUV in the same Rmd. Library-size factors come from TMM on the filtered gene matrix and are transferred to each repeat DGEList; no RUV covariates enter the design. TMM computed directly on the repeat matrix is avoided because TDP-43 perturbation moves many repeat subfamilies in concert, which violates TMM's "most features not DE" assumption and biases repeat logFC toward zero. The gene-library-size path serves as a transparent baseline against which the RUV-normalized analysis is compared in the integrative Rmd. Outputs are `gene_lib_repeat_de.rds` (list of topTags by featureset x contrast) and per-(featureset, contrast) TSVs named `gene_lib_{featureset}_{label}.tsv`. The helper is `run_gene_lib_repeat_de` in `workflow/scripts/gene_lib_size_common.R`, sourced alongside `ruv_common.R`.

### Single-cell pseudo-bulk RUVg (polymenidou_tdp_oe_cultures_scrnaseq)

The `render_ruv_sc` rule dispatches on the `kind: sc` field in the dataset config and renders `paper/scripts/ruv_polymenidou_tdp_oe_cultures_scrnaseq_report.Rmd`. Inputs are the STARsolo raw matrices produced by the sc method pipeline (`starsolo_base`), the GEO cell metadata file, and a sample metadata TSV. No counts are re-derived from raw reads. Pseudo-bulking is per (sample, cluster) pair: each GEO-annotated cluster in each of the three samples contributes its own pseudo-bulk column, dropping any pair with fewer than `min_cells_per_col` cells (default 10). With three samples and 15+ annotated clusters the gene matrix ends up with tens of columns rather than 6, which gives the design comfortable residual degrees of freedom. The binary test contrast stays cluster 12 versus any non-12 cluster, encoded as `group = factor(is_cluster12)` with `TRUE` for cluster-12 columns. Sample of origin enters the design through `sample_block = factor(sample_id)`. The pseudo-bulk gene matrix is filtered to genes with count >= 10 in at least half the columns. A naive `~ sample_block + group` model ranks genes by p-value; the top `n_controls` (default 3000) with the largest p-values become empirical controls. RUVg is fit for k = 1..4. The operating k is fixed via `k_selected` (default 1). Repeat DE uses `~ sample_block + W + group` with the `groupTRUE` coefficient as the test. The same `workflow/scripts/ruv_common.R` helpers used by the bulk reports are reused here. Outputs mirror the bulk layout: ranked TSVs `ruv_k{k}_{mode}_{featureset}_cluster12_vs_other.tsv`, RDS snapshots (`ruv_state.rds`, `ruv_naive_gene_de.rds`, `ruv_repeat_de.rds`, `pseudobulk_counts.rds`, `geo_meta_annotated.rds`), `ruv_empirical_controls_{mode}.txt`, `ruv_rle_diagnostic_{mode}.tsv`, `ruv_summary.tsv`, and CSV exports under `<report_outdir>/csv/` (`cluster_sizes_per_sample.csv`, `group_sizes_per_sample.csv`, `pseudobulk_col_meta.csv`, `pseudobulk_gene_lib_sizes.csv`, `ruv_W_factors.csv`, `ruv_rle_diagnostic.csv`, `ruv_de_long.csv`, `ruv_summary.csv`). `ruv_W_factors.csv` and `pseudobulk_gene_lib_sizes.csv` carry `sample_id`, `cluster_num`, `cluster_name`, `n_cells`, and `is_cluster12` columns so downstream code can map pseudo-bulk columns back to their (sample, cluster) origin.

RUVg is paper-side only, matching the bulk convention: `workflow/scripts/gse230647_sc_report.Rmd` runs gene-library-normalized edgeR plus QC only, and all RUVg-normalized DE for this dataset lives in the paper-side Rmd. Control selection, RUVg fitting, and edgeR contrast code are executed under the paper Snakefile and consume only the method pipeline's STARsolo outputs.

The sc Rmd adds two analyses on top of the RUV cluster-12-vs-other contrast. First, a parallel gene-library-size repeat DE is run for the same cluster-12-vs-other contrast, with library-size factors transferred from the pseudo-bulk gene matrix via `run_gene_lib_repeat_de`. Output is `gene_lib_repeat_de.rds` and per-(mode, featureset) TSVs named `gene_lib_{mode}_{featureset}_cluster12_vs_other.tsv`. Second, a per-cluster-vs-rest sweep runs the same gene-library-size analysis for every cluster with pseudo-bulk support in at least two samples. The sweep exists to empirically test the cluster-12-as-OE-counterpart hypothesis; the integrative Rmd correlates each cluster's repeat logFC profile against Polymenidou bulk OE and ranks clusters by that correlation. Output is `per_cluster_gene_lib_repeat_de.rds`. RUV is not run per cluster because RUV-g empirical controls are selected as genes least associated with the contrast under test, which means a per-cluster RUV analysis would need fresh control ranking per contrast and would mix noise-suppression with the question under study; the gene-library-size path is preferred here.

Run:

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd paper
snakemake --use-conda --cores 4 --configfile configs/tdp43.yaml
```

### Cross-dataset TDP-43 synthesis

The paper aims at a single TDP-43 view combining three inputs, all rendered by `paper/Snakefile`:

- `polymenidou_tdp_kd_oe_cultures_bulk` (GSE230647 bulk): knockdown (`tdp43_kd vs nt_kd_ctrl`) and overexpression (`tdp43_oe vs oe_ctrl`).
- `lee_nuclear_tdp_retention_bulk` (GSE126543 bulk): sorted nuclei, `tdp43pos vs tdp43neg`.
- `polymenidou_tdp_oe_cultures_scrnaseq` (GSE230647 scRNA-seq): cluster-12 vs other cells.

Cluster 12 in the sc dataset matches the TDP-43 overexpression arm of the bulk datasets. Cells in cluster 12 express the HA-TDP43 transgene under DOX-on conditions, while the rest of the DOX-on sc population shows little to no transgene. So the cluster-12 contrast is a neuron-resolved analogue of the bulk `tdp43_oe vs oe_ctrl` comparison in `polymenidou_tdp_kd_oe_cultures_bulk`. The other DOX-on and DOX-off sc samples, and the KD/NT bulk arms, are treated separately.

The synthesis rule `render_integrative_tdp43` in `paper/Snakefile` consumes the three per-dataset `ruv_rds/` directories and produces `results/paper/tdp43/integrative_tdp43/integrative_tdp43_report.html` plus ranked CSVs and RDS snapshots under `integrative_rds/`. Two concordance matrices are reported: one on RUV-normalized repeat DE tables (`integrative_repeat_concordance_summary.csv`), one on gene-library-size repeat DE tables (`integrative_gene_lib_repeat_concordance_summary.csv`). The two are also written side by side in `integrative_ruv_vs_gene_lib_concordance_side_by_side.csv`. Comparing the two separates cross-study signal robust to normalization choice from signal that depends on RUV factor estimation.

Pairs reported across both matrices: `lee_tdppos vs poly_kd`, `lee_tdppos vs poly_oe`, `lee_tdppos vs sc_cluster12`, `poly_kd vs poly_oe`, `poly_oe vs sc_cluster12`, `poly_kd vs sc_cluster12`. Each pair is reported separately for genic and intergenic repeats. Columns: n (features inner-joined on feature_id), pearson, spearman, n_sig_both (significant in both contrasts at `fdr_cutoff`), n_concordant_sig (same sign), n_discordant_sig (opposite sign). Biology-stratified scatter panels and correlation tables follow: polyA competence, downstream A-run (polyA priming proxy), family age. The biology annotation is loaded from any per-dataset `biology_annotation.rds` and joined on `feature_id = gene_id`.

`polya_class` and `a_run_class` are orthogonal axes. `polya_class` is class-level RNA biology (LINE, SINE, LTR, Retroposon generate polyA-tailed RNA during their life cycle; other classes do not). `a_run_class` is locus-level genomic DNA context (fraction of loci with a downstream A-run >= 12 bp) and proxies the risk of internal oligo-dT priming. Signal concentrated in polyA-competent + low-a_run subfamilies is the cleanest; polyA-competent + medium/high-a_run subfamilies (L1HS, SVA_B-F, young LTRs) are biologically interesting but cannot be disambiguated from priming artifact using RNA-seq alone. `a_run_label` in `workflow/scripts/biology_checks.R` excludes `Simple_repeat`, `Low_complexity`, and `Satellite` classes (returns NA) because their repeat units are themselves A-containing so the flag fires tautologically and carries no information about internal-priming of real polyadenylated transcripts; the exclusion keeps the a_run facets focused on classes where the priming question is meaningful.

A per-cluster ranking panel tests the cluster-12 choice empirically. Per-cluster gene-library-size repeat DE tables from the sc Rmd are inner-joined against the Polymenidou bulk OE gene-library-size repeat DE table, and Pearson and Spearman correlations of logFC are computed per cluster. Output is `integrative_per_cluster_vs_bulk_oe_ranking.csv` plus a bar chart of Pearson r per cluster, sorted; cluster 12 is highlighted. If cluster 12 does not top the ranking, the choice of sc counterpart of bulk OE deserves re-examination before further inference.
- Bulk per-contrast ranked TSVs at `results/paper/tdp43/{lee_nuclear_tdp_retention_bulk,polymenidou_tdp_kd_oe_cultures_bulk}/ruv_k{k}_{featureset}_{contrast}.tsv`.
- sc per-(mode, featureset) ranked TSVs at `results/paper/tdp43/polymenidou_tdp_oe_cultures_scrnaseq/ruv_k{k}_{mode}_{featureset}_cluster12_vs_other.tsv` and the long-format CSV `csv/ruv_de_long.csv`.


### Replogle 2022 K562 essential Perturb-seq (`repogle_2022`)

A second paper-side use case reanalyses the K562 essential Perturb-seq screen from Replogle et al. 2022 (Cell, doi 10.1016/j.cell.2022.05.013) at higher repeat granularity than the published class-level `TE_ratio`. Detailed procedure, data sources, and selection rules are in [repogle_2022.md](repogle_2022.md).

The Replogle dataset includes a per-perturbation `TE_ratio` (fraction of UMIs in repetitive elements, aggregated across Alu, L1, MIR). Their analysis identifies exosome subunits, CPSF, Integrator, and Pol II machinery as the dominant TE-fraction-increasing perturbation classes. REclaim re-quantifies the same fastqs at locus, gene_id, family_id, and class_id granularity, supporting subfamily- or locus-level fingerprints per perturbation alongside a one-to-one validation against the published `TE_ratio`.

Entry point: `paper/Snakefile_repogle_2022`, driven by `paper/configs/repogle_2022.yaml`. The use case is split into three phases. Phase 1 (implemented) downloads three figshare files (KD6 manifest, K562 essential raw bulk h5ad, K562 essential raw singlecell h5ad), verifies size and md5, and produces three slim TSVs that downstream phases consume: `selected_perturbations.tsv` (top 30 perturbations by `TE_ratio` plus their published value, for validation), `cells_to_perturbation.tsv` (~6,200 cells joined to gemgroup, gene_transcript, sgID_AB), and `selected_libraries.tsv` (per-gemgroup fastq and BAM filenames). Only the slim TSVs are tracked; the 10 GB singlecell h5ad lives outside the repo. Selection is seeded for reproducibility.

Cell selection defaults are `top_n_perturbations: 30`, `n_controls: 1000`, `min_cells_per_perturbation: 50`. Top 30 covers exosome (EXOSC2-9, DIS3), CPSF and 3' end processing (SYMPK, CPSF1-4, CSTF3, WDR33, WDR82), Integrator (INTS2/4/5/8/9), Pol II + Mediator + PAF (POLR2 subunits, MED9/30, PAF1, CTR9, RUVBL1, SUPT5/6H), and spliceosome subunits (SNRPD1/2, SNRPF). Each major complex is represented by 5-10 subunits, giving within-complex consistency checks. The `TE_ratio` dynamic range is narrow (~2.5x between controls and the top hit), so pushing past top 50 enters the noise floor. The 1000 control sample is comparable to the per-top-hit cell count (~170 cells/perturbation on average) and provides a stable pseudobulk baseline. A `gemgroup_subset` config knob restricts to a subset of the 47 KD6 libraries for development; setting it empty uses all libraries.

Phase 2 (alignment, run separately on the compute machine) reuses `pipeline_type: sc` with a new `workflow/configs/repogle_kd6_sc.yaml`. Aligners are restricted to STARsolo only (kallisto and salmon are skipped, halving alignment compute) since the simulation benchmark already favors STAR for repeats. Both `unique` and `multi` (EM) multimapper modes are kept: unique skews toward older, uniquely-mappable repeats and gives the cleanest baseline comparable to the TDP-43 RUVg pipeline; multi/EM redistributes multimappers and recovers young-element signal (L1HS, AluY, SVA). Cross-mode agreement on per-perturbation TE response is a robustness check; mode-specific divergence localizes the response to particular subfamily age strata. Once the simulation benchmark is rerun on barbara, the operating mode for downstream figures can be picked accordingly and the other mode kept in supplement.

Phase 3 (implemented in `paper/Snakefile_repogle_2022` and `paper/scripts/ruv_repogle_2022_perturbseq_report.Rmd`) builds per-(perturbation, gem_group) pseudobulks by summing UMIs across cells of each perturbation within each gemgroup. This gives ~30-48 columns per perturbation plus matched control_pool columns from the same gemgroups, which lets the per-contrast design fit `~ gem_group + W + group` with proper residual degrees of freedom (a single per-perturbation pseudobulk would have been 1 vs 1 against the control pool and degenerate). DE under `~ W + group` with RUVg fitted on the gene matrix using the existing `workflow/scripts/ruv_common.R` helpers, mirroring the TDP-43 bulk RUVg pattern. The validation step correlates REclaim's per-perturbation `te_total_fraction` against `replogle_te_ratio` from the bulk h5ad. Pseudobulk is chosen over single-cell because Replogle's `TE_ratio` is itself per-perturbation pseudobulk (one-to-one comparison required), per-cell TE quantification at 10x 3' depth is below per-cell noise for the relevant dynamic range, and the existing paper-side machinery is pseudobulk-based.

Reproducibility note: figshare URLs and md5s are pinned in the config. The `_01` suffix in Replogle's filenames is a chunk index; only `_01` exists for K562 essential as of figshare deposit version 2 (2022). A re-versioned deposit would change the URLs and md5s and require a config update.
