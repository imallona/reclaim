<img src="docs/reclaim-logo.svg" alt="REclaim" height="20"/>

**Repetitive element quantification in bulk and single-cell RNA-seq.**

REclaim (REpetitive Elements + claim) is a Snakemake pipeline for systematic analysis of repeat RNA expression. The name foregrounds what the tool does: bring repetitive elements back into analyses that routinely discard them.

## Status

Active development (2024-2026). Previous exploratory version: [v0.1](https://github.com/imallona/repeats/releases/tag/v0.1).

## Acknowledgements

We are extremely grateful to SNF for their funding in ca. 2020.

## Overview

A single `Snakefile` drives five pipeline modes selected by `pipeline_type` in the config:

- `simulation` - simulate reads from repeat loci (SmartSeq2 or 10x Chromium), align, and evaluate against ground truth
- `bulk` - download FASTQs from SRA (paired or single, set via `real_data.library_layout`), align with STAR/kallisto/salmon, quantify repeats, render a differential expression report
- `sc` - download 10x Chromium FASTQs from SRA, run STARsolo and kallisto|bustools, render a single-cell report
- `noise_report` - render a noise sweep HTML report across multiple simulation runs
- `de_simulation` - count-level power benchmark for repeat differential expression under varying gene and repeat library sizes (see [docs/de_simulations.md](docs/de_simulations.md))

For method details see [docs/methods.md](docs/methods.md). Dataset-specific paper reanalyses live under [paper/](paper/) with their own Snakefile. Workflow diagrams (mermaid, renders on GitHub) are in [docs/diagrams.md](docs/diagrams.md).

## Requirements

- Snakemake 8+
- conda / mamba (for --use-conda)

All other dependencies are installed automatically via per-rule conda environments.

Exact package pins (platform-locked explicit exports) are in `workflow/envs/explicit/`:

| File | Environment | Key packages |
|---|---|---|
| `repeats_star.txt` | repeats_star | STAR, samtools, bedtools, pigz, gffread |
| `repeats_kallisto.txt` | repeats_kallisto | kallisto, bustools |
| `repeats_alevin.txt` | repeats_alevin | salmon |
| `repeats_bowtie2.txt` | repeats_bowtie2 | bowtie2, samtools, subread (featureCounts) |
| `repeats_umi_tools.txt` | repeats_umi_tools | umi_tools, samtools >= 1.12, pysam |
| `repeats_evaluation.txt` | repeats_evaluation | python, scipy, pysam |
| `rmarkdown.txt` | rmarkdown | R, rmarkdown, ggplot2, patchwork |

Two additional envs without explicit pins: `workflow/envs/edger.yaml` (used by the baseline bulk edgeR report) and `paper/envs/ruvseq.yaml` (used by the paper TDP-43 RUVg report, see Paper analyses below).

## Running

### Via Makefile (recommended)

A `Makefile` at the project root orchestrates all configs and renders reports.

```
# Run everything (all simulations + noise sweep + reports)
make

# Individual targets
make simulation_smartseq2      # SmartSeq2 base simulation
make simulation_chromium       # Chromium base simulation
make noise_smartseq2           # SmartSeq2 noise sweep (0%, 1%, 5%, 10% mutation rate)
make noise_chromium            # Chromium noise sweep
make report_noise_smartseq2    # Render SmartSeq2 noise sweep HTML report
make report_noise_chromium     # Render Chromium noise sweep HTML report
make de_polymenidou_bulk       # Count-level DE library size and normalization check (TSV/RDS, heatmap PDF, HTML report)
make help                      # List all targets
```

Tune parallelism with `make CORES=20`. The Makefile activates the `snakemake` conda environment automatically.

### Memory-aware scheduling

Alignment rules and the report rule declare `resources: mem_mb` from `resources.max_mem_mb` (override the report budget with `resources.report_mem_mb`). Pass `--resources mem_mb=N` to cap the concurrent sum, for example `--resources mem_mb=16000` on a 16 GB host.

### Manually

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd workflow

# simulation
snakemake --configfile configs/simulation_smartseq2.yaml --use-conda --cores 10
snakemake --configfile configs/simulation_chromium.yaml  --use-conda --cores 10

# bulk reanalysis
snakemake --configfile configs/gse230647_bulk.yaml --use-conda --cores 10
snakemake --configfile configs/gse126543_bulk.yaml --use-conda --cores 10

# single-cell reanalysis
snakemake --configfile configs/gse230647_sc.yaml --use-conda --cores 10
```

The simulation evaluation report is written to `{base}/evaluation/evaluation_report.html`. Bulk and SC reports are written to `{base}/{report.output}`.

## Configuration

Simulation configs are under `workflow/configs/`:

| File | Technology | Cells | Expressed loci/cell | Mutation rate |
|---|---|---|---|---|
| `simulation_smartseq2.yaml` | SmartSeq2 | 100 | 500 | 0.1% |
| `simulation_chromium.yaml`  | 10x Chromium | 100 | 500 | 0.1% |
| `simulation_smartseq2_noise_0pct.yaml` | SmartSeq2 | 100 | 500 | 0% |
| `simulation_smartseq2_noise_1pct.yaml` | SmartSeq2 | 100 | 500 | 1% |
| `simulation_smartseq2_noise_5pct.yaml` | SmartSeq2 | 100 | 500 | 5% |
| `simulation_smartseq2_noise_10pct.yaml` | SmartSeq2 | 100 | 500 | 10% |
| `simulation_chromium_noise_0pct.yaml` | 10x Chromium | 100 | 500 | 0% |
| `simulation_chromium_noise_1pct.yaml` | 10x Chromium | 100 | 500 | 1% |
| `simulation_chromium_noise_5pct.yaml` | 10x Chromium | 100 | 500 | 5% |
| `simulation_chromium_noise_10pct.yaml` | 10x Chromium | 100 | 500 | 10% |

Unused real-data configs have been moved to `workflow/configs/old/`.

### Bulk RNA-seq reanalysis

Two publicly available bulk RNA-seq datasets are reanalysed for repeat element quantification:

| GSE | Description | Design | Config |
|---|---|---|---|
| GSE126543 | NeuN+ neuronal nuclei, human frontal cortex (FTD/ALS), TDP-43 nuclear positive vs negative, 7 donors | paired by donor | `configs/gse126543_bulk.yaml` |
| GSE230647 | iPSC-derived neural cultures, TDP-43 OE (HA-TDP43 DOX on vs off, n=4) and KD (TDP-43 shRNA vs NT shRNA, n=4) | two independent experiments | `configs/gse230647_bulk.yaml` |

Both pipelines download paired-end FASTQs from SRA (`library_layout: paired`). Reads align with STAR, kallisto, and salmon. Repeats are quantified at gene_id and family_id granularity. An edgeR differential expression report is rendered at the end.

### Single-cell reanalysis

A single-cell analysis covers the TDP-43-HA overexpression scRNA-seq experiment from GSE230647. Three DOX-on 10x Chromium 3' v3 samples are included (two sequencing lanes each): GSM7230453 (2-week culture), GSM7230454 and GSM7230455 (4-week culture, replicates A and B). The DOX-off control and unlabeled samples are excluded.

| Item | Value |
|---|---|
| Config | `configs/gse230647_sc.yaml` |
| Sample metadata | `configs/gse230647_sc_sample_metadata.tsv` |
| Report script | `scripts/gse230647_sc_report.Rmd` |
| Output | `../results/gse230647_sc/gse230647_sc_report.html` |

The pipeline downloads each SRR with `fasterq-dump --include-technical`. For 10x Chromium, file `_2` carries the CB+UMI read and file `_3` carries the cDNA read (set via `sc_cbumi_index` and `sc_cdna_index`). STARsolo (CB_UMI_Simple, 10x v3 geometry) and kallisto|bustools run against genic and intergenic repeat annotations. The GEO supplementary file `GSE230647_single_cell_metadata_tdp43ha.txt.gz` maps each cell barcode to its cluster (17 clusters). Cluster 12 marks the TDP-43-HA-expressing cells.

The report builds pseudo-bulk count matrices by aggregating STARsolo raw counts for cluster-12 cells against all other annotated cells. edgeR (quasi-likelihood, sample as blocking factor) identifies differentially expressed repeat elements. A second, RUVg-normalized pass on the same pseudo-bulk matrices provides a sensitivity check.

Indices are shared with the bulk GSE230647 pipeline via `indices_base: ../results/shared`, so building them once is sufficient for both analyses.

Key config parameters (shared across pipeline types):

- `pipeline_type`: `simulation`, `bulk`, `sc`, `noise_report`, or `de_simulation`.
- `base`: run-specific output directory.
- `indices_base`: shared directory for aligner indices. Multiple runs can share one index.
- `feature_sets`: which repeat subsets to quantify (`repeats`, `genic_repeats`, `intergenic_repeats`).
- `granularities`: aggregation levels (`gene_id`, `family_id`, `class_id`).

Simulation-specific:

- `simulation.mutation_rate`: per-base substitution rate applied to simulated reads.
- `aligner_params.{aligner}.multimapper_modes`: `unique` or `multi` (EM).

Bulk-specific (`real_data`):

- `library_layout`: `paired` (default) or `single`.

SC-specific (`real_data`):

- `sc_cbumi_index`: fasterq-dump file number containing CB+UMI (default `2`).
- `sc_cdna_index`: fasterq-dump file number containing cDNA (default `3`).
- `cb_length`, `umi_length`: barcode and UMI lengths. Single source of truth for both STARsolo (`--soloCBlen`/`--soloUMIlen`) and kallisto bus (`-x 0,0,cb:0,cb,cb+umi:1,0,0` custom geometry), so both aligners stay in sync if you change the chemistry.

## Paper analyses

A separate Snakefile under `paper/` drives manuscript-specific reanalyses across datasets produced by the method pipelines. It is kept separate from `workflow/Snakefile` because paper analyses consume outputs of multiple method runs and should not retrigger method jobs via Snakemake provenance tracking. The paper rule adds a RUVg-normalized edgeR report per dataset alongside the baseline method reports. It does not replace the baseline reports and does not re-run alignment.

Dataset keys in `paper/configs/tdp43.yaml` are named by last author plus study label plus assay kind. Each entry also records the source GEO accession and the method-pipeline config that produced its inputs.

| Paper dataset key | GSE | Assay | Source method config | Paper report script |
|---|---|---|---|---|
| `lee_nuclear_tdp_retention_bulk` | GSE126543 | bulk | `workflow/configs/gse126543_bulk.yaml` | `paper/scripts/ruv_lee_nuclear_tdp_retention_bulk_report.Rmd` |
| `polymenidou_tdp_kd_oe_cultures_bulk` | GSE230647 | bulk | `workflow/configs/gse230647_bulk.yaml` | `paper/scripts/ruv_polymenidou_tdp_kd_oe_cultures_bulk_report.Rmd` |
| `polymenidou_tdp_oe_cultures_scrnaseq` | GSE230647 | scRNA-seq | `workflow/configs/gse230647_sc.yaml` | `paper/scripts/ruv_polymenidou_tdp_oe_cultures_scrnaseq_report.Rmd` |

| Item | Value |
|---|---|
| Entry point | `paper/Snakefile` |
| Config | `paper/configs/tdp43.yaml` |
| Conda env | `paper/envs/ruvseq.yaml` (R, edgeR, RUVSeq, EDASeq) |
| Output base | `results/paper/tdp43/{dataset}/` |

The bulk reports run a per-dataset RUVg-normalized edgeR analysis on STAR unique counts. Empirical controls come from the gene matrix by edgeR F-test ranking. RUVg is fit on those controls. The per-sample W factors go into `~ W + group` (GSE126543) or `~ sample_block + W + group` (GSE230647 four-group design) models for every repeat featureset and granularity. See [docs/methods.md](docs/methods.md) for the full procedure.

The single-cell RUVg report (`polymenidou_tdp_oe_cultures_scrnaseq`) is the paper-side counterpart of the sc method pipeline's RUVg block. Each (sample, GEO-annotated cluster) pair becomes its own pseudo-bulk column built from the STARsolo raw matrices at `results/gse230647_sc/starsolo/`, giving tens of columns across the three samples. The binary contrast is cluster 12 versus any non-12 cluster, encoded as `group = factor(is_cluster12)` in a `~ sample_block + W + group` design with `sample_block = factor(sample_id)`. Cluster 12 captures the DOX-on TDP-43-HA overexpressing neurons and acts as the neuron-resolved counterpart of the bulk `tdp43_oe vs oe_ctrl` contrast in `polymenidou_tdp_kd_oe_cultures_bulk`.

A planned cross-dataset synthesis rule will pull CSVs from `results/paper/tdp43/{lee_nuclear_tdp_retention_bulk,polymenidou_tdp_kd_oe_cultures_bulk,polymenidou_tdp_oe_cultures_scrnaseq}/` and merge them into a long-format table for cross-contrast intersections (FDR < 0.05 concordance, logFC direction agreement). Not yet wired.

Run after the relevant method pipelines complete:

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd paper
snakemake --use-conda --cores 4 --configfile configs/tdp43.yaml
```

### Replogle 2022 K562 essential Perturb-seq (`repogle_2022`)

A second paper-side use case reanalyses the Replogle et al. 2022 K562 essential Perturb-seq screen at higher repeat granularity than the published class-level `TE_ratio`. See [docs/repogle_2022.md](docs/repogle_2022.md) for the full procedure, data sources (figshare DOIs 20022944 and 20029387), cell selection rules, and rationale.

Phase 1 retrieves figshare metadata and selects cells:

```
cd paper
snakemake -s Snakefile_repogle_2022 --use-conda --cores 2 \
    --configfile configs/repogle_2022.yaml
```

Outputs three slim TSVs under `results/paper/repogle_2022/data/`: `selected_perturbations.tsv` (top 30 perturbations by `TE_ratio` plus their published value, for validation), `cells_to_perturbation.tsv` (~6,200 cells joined to perturbation and gemgroup), `selected_libraries.tsv` (per-gemgroup fastq and BAM filenames parsed from the KD6 manifest).

Phase 2 (planned): SRA fastq download and STARsolo alignment via a new `workflow/configs/repogle_kd6_sc.yaml` reusing `pipeline_type: sc`. Aligners restricted to STARsolo with both `unique` and `multi` (EM) multimapper modes. Other aligners are skipped to halve compute since the simulation benchmark already shows STAR is best for repeats.

Phase 3 (planned): per-perturbation pseudobulk DE under `~ W + group` with RUVg, mirroring the TDP-43 paper-side workflow. Validation step correlates REclaim's per-perturbation `te_total_fraction` against Replogle's `replogle_te_ratio` from `selected_perturbations.tsv`.

## External tool benchmark

A separate config flag `external_benchmark: true` activates a comparison of REclaim against published repeat-quantification tools on the same simulation BAMs and ground truth, scored by the same `evaluate.py`. SmartSeq2 runs TEtranscripts and SQuIRE (both bulk-style, per-cell BAM treated as a sample). Chromium runs scTE (single-cell). See [docs/external_tool_benchmark.md](docs/external_tool_benchmark.md) for tool choices, version pins, and harmonisation contract.

```
make external_benchmark_smartseq2 CORES=N    # TEtranscripts and SQuIRE on the SmartSeq2 simulation
make external_benchmark_chromium  CORES=N    # scTE on the Chromium simulation
```

Outputs land at `results/{simulation_*}/external_benchmark/external_benchmark_report.html` alongside the per-tool harmonised counts and accuracy tables. The report compares REclaim quantifiers and the external tools side by side at every granularity, with a coverage table that marks each cell as native or rolled up.

## Implementation notes

Chromium normalization (kallisto, alevin) uses sparse accumulators to avoid allocating a dense cells x features matrix. Only non-zero (cell_index, count) pairs are stored per feature group, keeping memory proportional to expressed pairs rather than O(features x cells).

Bowtie2 Chromium counting streams the CB-tagged BAM once and accumulates per-(barcode, locus) counts directly, without splitting into per-cell BAM files.

See [docs/methods.md](docs/methods.md) for full details.

## Testing

Unit tests, integration tests, and Snakemake dry-run tests are in the `test/` directory. See [test/README.md](test/README.md) for details on design, coverage, and how to run the tests.

## Contact

izaskun dot mallona dot work at gmail.com

## License

GNU General Public License (GPL)
