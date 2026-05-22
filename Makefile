SHELL   := /bin/bash
.DEFAULT_GOAL := all
CORES   ?= 8
WF      := workflow
RESULTS := results

CONDA_RUN := source ~/miniconda3/etc/profile.d/conda.sh && conda activate snakemake &&
SM        := cd $(WF) && $(CONDA_RUN) ulimit -Sv 209715200 && snakemake --use-conda --cores $(CORES) --rerun-triggers mtime
SM_REPORT := cd $(WF) && $(CONDA_RUN) snakemake -s Snakefile_noise_report --use-conda --cores 1

# noise sweep eval dirs (relative to workflow/, comma-separated for the Rmd param)
NOISE_EDIRS_SS2 := $(CURDIR)/$(RESULTS)/simulation_smartseq2_noise_0pct/evaluation,$(CURDIR)/$(RESULTS)/simulation_smartseq2_noise_1pct/evaluation,$(CURDIR)/$(RESULTS)/simulation_smartseq2_noise_5pct/evaluation,$(CURDIR)/$(RESULTS)/simulation_smartseq2_noise_10pct/evaluation
NOISE_EDIRS_CHR := $(CURDIR)/$(RESULTS)/simulation_chromium_noise_0pct/evaluation,$(CURDIR)/$(RESULTS)/simulation_chromium_noise_1pct/evaluation,$(CURDIR)/$(RESULTS)/simulation_chromium_noise_5pct/evaluation,$(CURDIR)/$(RESULTS)/simulation_chromium_noise_10pct/evaluation

# output HTML files (relative to project root)
NOISE_REPORT_SS2 := $(RESULTS)/noise_sweep_smartseq2.html
NOISE_REPORT_CHR := $(RESULTS)/noise_sweep_chromium.html

.PHONY: all \
        simulation_smartseq2 simulation_chromium \
        noise_smartseq2 noise_smartseq2_0pct noise_smartseq2_1pct noise_smartseq2_5pct noise_smartseq2_10pct \
        noise_chromium noise_chromium_0pct noise_chromium_1pct noise_chromium_5pct noise_chromium_10pct \
        report_noise_smartseq2 report_noise_chromium reports_noise \
        de_polymenidou_bulk de_coordinated_derepression \
        external_benchmark_smartseq2 external_benchmark_chromium \
        help

# -----------------------------------------------------------------------
# base simulations
# -----------------------------------------------------------------------

simulation_smartseq2:
	$(SM) --configfile configs/simulation_smartseq2.yaml

simulation_chromium:
	$(SM) --configfile configs/simulation_chromium.yaml

# -----------------------------------------------------------------------
# noise sweep - SmartSeq2
# note: genome refs are reused from the base simulation_smartseq2 run;
#       alignment indices are shared via results/shared/ (never rebuilt)
# -----------------------------------------------------------------------

noise_smartseq2_0pct:
	$(SM) --configfile configs/simulation_smartseq2_noise_0pct.yaml

noise_smartseq2_1pct:
	$(SM) --configfile configs/simulation_smartseq2_noise_1pct.yaml

noise_smartseq2_5pct:
	$(SM) --configfile configs/simulation_smartseq2_noise_5pct.yaml

noise_smartseq2_10pct:
	$(SM) --configfile configs/simulation_smartseq2_noise_10pct.yaml

noise_smartseq2: noise_smartseq2_0pct noise_smartseq2_1pct noise_smartseq2_5pct noise_smartseq2_10pct

# -----------------------------------------------------------------------
# noise sweep - Chromium
# -----------------------------------------------------------------------

noise_chromium_0pct:
	$(SM) --configfile configs/simulation_chromium_noise_0pct.yaml

noise_chromium_1pct:
	$(SM) --configfile configs/simulation_chromium_noise_1pct.yaml

noise_chromium_5pct:
	$(SM) --configfile configs/simulation_chromium_noise_5pct.yaml

noise_chromium_10pct:
	$(SM) --configfile configs/simulation_chromium_noise_10pct.yaml

noise_chromium: noise_chromium_0pct noise_chromium_1pct noise_chromium_5pct noise_chromium_10pct

# -----------------------------------------------------------------------
# noise sweep reports
# file targets: only (re)rendered when the HTML does not yet exist;
# use  make -B report_noise_smartseq2  to force a re-render
# -----------------------------------------------------------------------

$(NOISE_REPORT_SS2):
	mkdir -p $(RESULTS)
	$(SM_REPORT) --config eval_dirs='$(NOISE_EDIRS_SS2)' output_file='$(CURDIR)/$(NOISE_REPORT_SS2)'

$(NOISE_REPORT_CHR):
	mkdir -p $(RESULTS)
	$(SM_REPORT) --config eval_dirs='$(NOISE_EDIRS_CHR)' output_file='$(CURDIR)/$(NOISE_REPORT_CHR)'

report_noise_smartseq2: $(NOISE_REPORT_SS2)
report_noise_chromium:  $(NOISE_REPORT_CHR)
reports_noise: report_noise_smartseq2 report_noise_chromium

# -----------------------------------------------------------------------
# DE library size and normalization check (count-level)
# requires gene and repeat count TSVs from a prior bulk run; see
# workflow/configs/de_simulations_polymenidou_bulk.yaml for the input paths
# -----------------------------------------------------------------------

de_polymenidou_bulk:
	$(SM) --configfile configs/de_simulations_polymenidou_bulk.yaml

de_coordinated_derepression:
	$(SM) --configfile configs/de_simulations_coordinated_derepression.yaml

# -----------------------------------------------------------------------
# External tool benchmark (TEtranscripts, scTE, optional SQuIRE)
# requires the corresponding base simulation to have been run first so the
# STAR-aligned BAMs and ground truth exist
# -----------------------------------------------------------------------

external_benchmark_smartseq2:
	$(SM) --configfile configs/external_benchmark_smartseq2.yaml

external_benchmark_chromium:
	$(SM) --configfile configs/external_benchmark_chromium.yaml

# -----------------------------------------------------------------------
# all
# -----------------------------------------------------------------------

all: simulation_smartseq2 simulation_chromium noise_smartseq2 noise_chromium reports_noise

# -----------------------------------------------------------------------
# help
# -----------------------------------------------------------------------

help:
	@echo "Usage: make [target] [CORES=N]   (default CORES=$(CORES))"
	@echo ""
	@echo "Base simulation runs:"
	@echo "  simulation_smartseq2             SmartSeq2 base run (mutation_rate 0.1%)"
	@echo "  simulation_chromium              Chromium base run  (mutation_rate 0.1%)"
	@echo ""
	@echo "Noise sweep runs (base simulation must have been run first):"
	@echo "  noise_smartseq2                  all 4 noise levels for SmartSeq2"
	@echo "  noise_smartseq2_{0,1,5,10}pct   individual SmartSeq2 noise level"
	@echo "  noise_chromium                   all 4 noise levels for Chromium"
	@echo "  noise_chromium_{0,1,5,10}pct    individual Chromium noise level"
	@echo ""
	@echo "Reports (file targets; use -B to force re-render):"
	@echo "  report_noise_smartseq2           $(NOISE_REPORT_SS2)"
	@echo "  report_noise_chromium            $(NOISE_REPORT_CHR)"
	@echo "  reports_noise                    both noise sweep reports"
	@echo ""
	@echo "DE library size and normalization check:"
	@echo "  de_polymenidou_bulk              default scenario: count-level power benchmark for the gse230647 bulk gene/repeat counts"
	@echo "                                   (TSV/RDS, heatmap PDF, and HTML report under results/de_simulations_polymenidou_bulk/de_simulations/)"
	@echo "  de_coordinated_derepression      coordinated-derepression scenario: large fraction of repeats moving in same direction;"
	@echo "                                   tests TMM-on-repeats compression vs gene-library-size transfer (recovered-vs-input logFC)"
	@echo "                                   (outputs under results/de_simulations_coordinated_derepression/de_simulations/)"
	@echo ""
	@echo "External tool benchmark (TEtranscripts, scTE, optional SQuIRE):"
	@echo "  external_benchmark_smartseq2     run TEtranscripts on the SmartSeq2 simulation BAMs and score against ground truth"
	@echo "                                   (outputs under results/simulation_smartseq2/external_benchmark/)"
	@echo "  external_benchmark_chromium      run scTE on the Chromium simulation BAM and score against ground truth"
	@echo "                                   (outputs under results/simulation_chromium/external_benchmark/)"
	@echo ""
	@echo "  all                              run everything in sequence"
