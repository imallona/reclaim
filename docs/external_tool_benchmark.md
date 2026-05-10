# external_tool_benchmark

Head-to-head benchmark of REclaim against published repeat-quantification tools on the same simulated data, scored by the same `evaluate.py` against the same ground truth at the same granularities. The result feeds the manuscript's external-comparison panel.

## Why a benchmark on the existing simulation

REclaim's evaluation report at `results/{sim}/evaluation/evaluation_report.html` already scores STARsolo (unique and multi/EM), kallisto, alevin, and bowtie2 against the simulated ground truth at locus, gene_id, family_id, and class_id. The external benchmark reuses the same fastqs and BAMs and feeds external-tool counts through the same `evaluate.py` so every number is directly comparable.

## Tools included

| Tool | Modality | Native granularity | Why included | Conda env |
|---|---|---|---|---|
| TEtranscripts (TEcount) | bulk | subfamily (gene_id) | Most cited bulk tool. Per-cell SmartSeq2 BAMs are treated as per-sample bulk replicates of one group. | `workflow/envs/tetranscripts.yaml` |
| scTE | single cell | family | The natural sc counterpart. Consumes a STARsolo BAM with CB and UB tags and emits a cell x family matrix. | `workflow/envs/scte.yaml` |
| SQuIRE | bulk | locus | Optional: pinned to Python 3.6 and STAR 2.5.3a; the upstream is unmaintained. Drop it from the config when the conda env fails to build. | `workflow/envs/squire.yaml` |

## Tools deliberately skipped

REdiscoverTE bundles a kallisto-based custom transcriptome plus a hand-curated subfamily mapping. Integrating it cleanly would mean a second pseudoalignment branch parallel to REclaim's existing kallisto path, which inflates scope without adding a new methodological axis. The REclaim `kallisto` quantifier is already in the comparison.

scTEcount, scNanoGPS, and other tool variants released after 2023 were considered. None met both criteria: published with a versioned conda recipe and accepting BAM input, so the benchmark stays on stable ground.

## How the comparison is wired

For each tool, three rules in `workflow/modules/external_tools_benchmark.snmk` produce the final per-(tool, granularity) accuracy table:

1. `run_<tool>` invokes the upstream tool on the existing simulation BAMs. TEcount runs once per simulated SmartSeq2 cell (one BAM each); scTE runs once on the multi-cell Chromium BAM. Outputs land under `{base}/external_benchmark/{tool}/raw/`.
2. `harmonize_external` calls `workflow/scripts/harmonize_external_counts.py` to convert the tool's native output to the REclaim feature x cell TSV format. The harmonisation step asserts that at least `external_min_overlap_fraction` of the tool's feature IDs map to a known REclaim feature in the locus_map; below that, it raises and the tool is dropped from the report.
3. `score_external_benchmark` runs the standard `workflow/scripts/evaluate.py` on the harmonised counts. The output TSVs follow the same naming convention as the REclaim per-aligner outputs, so the comparison report can read both with one ingestion path.

The render rule then reads `summary_global_metrics.tsv` from both `external_benchmark/` and `evaluation/` and produces side-by-side accuracy bars per quantifier, granularity, and class.

## Ground-truth contract

`evaluate.py` consumes a per-locus ground truth TSV with columns `cell_id, locus_id, repeat_id, family_id, class_id, true_count` and aggregates inside the script to the requested granularity. Every external tool feeds in at the same granularity contract so there is no double-aggregation hidden in the comparison.

## Tool versions

Pinned in the conda env yamls:

| Tool | Version | Notes |
|---|---|---|
| TEtranscripts | 2.2.3 | pip-installed; depends on pysam and HTSeq |
| scTE | 1.0.0 | pip-installed; uses anndata for output |
| SQuIRE | 0.9.9.92 | pip-installed; pins Python 3.6 and STAR 2.5.3a |

Cite these versions in the methods section of the manuscript.

## Running

```
## After the base simulations have produced BAMs and ground truth:
make simulation_smartseq2 CORES=N
make simulation_chromium CORES=N

## Then the external benchmark on each:
make external_benchmark_smartseq2 CORES=N
make external_benchmark_chromium CORES=N
```

Outputs land at `results/simulation_smartseq2/external_benchmark/external_benchmark_report.html` and `results/simulation_chromium/external_benchmark/external_benchmark_report.html`.

## Harmonisation logic

`harmonize_external_counts.py` accepts `--tool tetranscripts|scte|squire`, parses the tool-native format, and writes a feature_id x cell TSV.

For TEtranscripts, the feature ID is `subfamily:family:class`; the helper strips to the subfamily component (matches REclaim gene_id) and rolls up to family or class via the locus_map for the requested granularity.

For scTE, the feature ID is family-level; gene_id rollup is not possible (the tool aggregates upstream of subfamily) so the harmoniser warns and emits family rows. Family and class outputs are exact.

For SQuIRE, locus rows aggregate to subfamily then up to family or class via the locus_map.

Every map step counts dropped features. The drop-count is logged to stderr and visible in the snakemake log.

## What can go wrong, and what the asserts catch

| Failure mode | Symptom | Caught by |
|---|---|---|
| Tool's feature naming drifts | Below-threshold overlap | `check_overlap` assert in `harmonize_external_counts.py` |
| Sample list mismatched | TEtranscripts header missing requested column | `parse_tetranscripts_table` assert |
| Locus map missing | Empty mapping, all features dropped | `aggregate_to_granularity` returns empty rolled dict; `write_count_matrix` asserts |
| SQuIRE conda env fails | Snakemake rule fails before harmoniser | drop SQuIRE from `external_tools` list |

The unit tests in `test/unit/test_harmonize_external_counts.py` and `test/unit/test_external_benchmark_helpers.py` lock these contracts.

## Scope

The benchmark is on simulation data only. It is not designed to score external tools on the TDP-43 application data (no ground truth there). For the application-side comparison the manuscript uses the within-Polymenidou KD/OE concordance as the empirical readout.
