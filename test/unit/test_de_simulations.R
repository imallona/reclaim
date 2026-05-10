#!/usr/bin/env Rscript
## Unit tests for workflow/scripts/de_simulations.R
##
## Sources the helpers from the production script and uses them on small
## synthetic matrices. Exits non-zero on any failure. Designed to require only
## edgeR, RUVSeq, EDASeq, Biobase, and matrixStats from the edger conda env.

suppressPackageStartupMessages({
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(Biobase)
  library(matrixStats)
})

script_path <- Sys.getenv("DE_SIM_SCRIPT", unset = "")
if (!nzchar(script_path)) {
  script_path <- file.path(getwd(), "workflow", "scripts", "de_simulations.R")
}
stopifnot(file.exists(script_path))

env <- new.env(parent = globalenv())
source(script_path, local = env)


report <- function(name, ok) {
  cat(sprintf("%s %s\n", if (ok) "PASS" else "FAIL", name))
  if (!ok) quit(status = 1)
}


## simulate_count_matrix: shape and non-negativity
m <- env$simulate_count_matrix(mu = c(10, 50, 200, 1000),
                               dispersion = rep(0.1, 4),
                               n_samples = 6, lib_scale = 1, seed = 1)
report("simulate_count_matrix shape and >= 0",
       is.matrix(m) && nrow(m) == 4 && ncol(m) == 6 && all(m >= 0))

## lib_scale increases mean counts roughly proportionally
m1 <- env$simulate_count_matrix(rep(100, 200), rep(0.1, 200), 10, 1, seed = 2)
m5 <- env$simulate_count_matrix(rep(100, 200), rep(0.1, 200), 10, 5, seed = 2)
report("lib_scale of 5 raises mean counts above 3x of lib_scale 1",
       mean(m5) > 3 * mean(m1))

## planted DE recovered by score_calls when FDR is small only on planted rows
features <- paste0("f", seq_len(100))
fake_tt <- data.frame(FDR = c(rep(0.01, 10), rep(0.5, 90)), row.names = features)
sc <- env$score_calls(fake_tt, planted_features = features[1:10], fdr = 0.05)
report("score_calls perfect recovery",
       sc$power == 1 && sc$fpr == 0 && sc$tp == 10 && sc$fn == 0)

## complete miss: significant rows are entirely outside the planted set
sc2 <- env$score_calls(fake_tt, planted_features = features[91:100], fdr = 0.05)
report("score_calls complete miss",
       sc2$power == 0 && sc2$fp == 10 && sc2$tn == 80)

## planted features that are not in tt_df count as FN, not silently dropped
truth_with_missing <- c(features[1:10], "missing_feature")
sc3 <- env$score_calls(fake_tt, planted_features = truth_with_missing, fdr = 0.05)
report("score_calls penalizes filtered-out planted features",
       sc3$power == 10 / 11 && sc3$fn == 1)

## pick_empirical_controls returns the requested number of unique features
set.seed(7)
counts <- matrix(rpois(50 * 8, lambda = 100), nrow = 50)
rownames(counts) <- paste0("g", seq_len(50))
colnames(counts) <- paste0("s", seq_len(8))
design <- model.matrix(~ factor(rep(c("A", "B"), each = 4)))
ctrl <- env$pick_empirical_controls(counts, design, n_controls = 10)
report("pick_empirical_controls returns 10 unique gene IDs",
       length(ctrl) == 10 && !anyDuplicated(ctrl) && all(ctrl %in% rownames(counts)))

## end-to-end run on a tiny synthetic problem
gene_params <- list(
  mu = rep(c(50, 200), length.out = 60),
  dispersion = rep(0.1, 60),
  feature_id = paste0("g", seq_len(60)))
repeat_params <- list(
  mu = rep(c(20, 80), length.out = 30),
  dispersion = rep(0.2, 30),
  feature_id = paste0("r", seq_len(30)))
samples <- paste0("s", seq_len(8))
condition <- rep(c("A", "B"), each = 4)
out <- env$run_one_grid_cell(
  gene_params, repeat_params, samples, condition,
  gene_lib_scale = 1, repeat_lib_scale = 1,
  fc = 5, n_de_repeats = 10, n_de_genes = 15,
  sigma_w = 0.3, n_controls = 20, k_max = 1, fdr = 0.05,
  n_iter = 2, seed_base = 42)
res <- out$scores
report("run_one_grid_cell returns scores frame with expected columns",
       is.list(out)
       && is.data.frame(res)
       && all(c("power", "fpr", "method", "iter") %in% colnames(res))
       && all(c("none", "TMM", "gene_lib", "RUVg_k1") %in% unique(res$method))
       && nrow(res) >= 8)

## select_de_features: signed_fraction = 1.0 means all positive logfc
fids <- paste0("f", seq_len(50))
plant_up <- env$select_de_features(fids, n_de = 20, fc = 3,
                                   signed_fraction = 1.0, seed = 1L)
report("select_de_features signed_fraction=1 yields all positive logfc",
       length(plant_up$idx) == 20
       && sum(plant_up$planted) == 20
       && all(plant_up$logfc[plant_up$idx] > 0)
       && all(plant_up$logfc[!plant_up$planted] == 0))

## select_de_features: signed_fraction = 0.0 means all negative
plant_dn <- env$select_de_features(fids, n_de = 20, fc = 3,
                                   signed_fraction = 0.0, seed = 1L)
report("select_de_features signed_fraction=0 yields all negative logfc",
       all(plant_dn$logfc[plant_dn$idx] < 0))

## select_de_features: signed_fraction = 0.5 splits roughly evenly
plant_mix <- env$select_de_features(fids, n_de = 50, fc = 3,
                                    signed_fraction = 0.5, seed = 2L)
n_pos <- sum(plant_mix$logfc > 0)
n_neg <- sum(plant_mix$logfc < 0)
report("select_de_features signed_fraction=0.5 splits 25/25",
       n_pos == 25 && n_neg == 25)

## select_de_features: class_weights bias selection toward weighted classes
class_id <- rep(c("LINE", "SINE", "LTR", "DNA"), length.out = 50)
weights <- c(LINE = 1, SINE = 1, LTR = 1, DNA = 0)
plant_w <- env$select_de_features(fids, n_de = 20, fc = 3,
                                  class_id = class_id, class_weights = weights,
                                  seed = 3L)
picked_classes <- class_id[plant_w$idx]
report("select_de_features class_weights respect zero-weighted class",
       !any(picked_classes == "DNA"))

## select_de_features: zero n_de returns empty plant set
plant_zero <- env$select_de_features(fids, n_de = 0, fc = 3, seed = 1L)
report("select_de_features n_de=0 returns no plant",
       length(plant_zero$idx) == 0
       && all(plant_zero$logfc == 0)
       && !any(plant_zero$planted))

## select_de_features: invalid signed_fraction triggers stopifnot
err1 <- tryCatch(env$select_de_features(fids, n_de = 5, fc = 3,
                                        signed_fraction = 1.5, seed = 1L),
                 error = function(e) e)
report("select_de_features rejects signed_fraction > 1",
       inherits(err1, "error"))

## select_de_features: invalid (all zero) class_weights triggers stopifnot
err2 <- tryCatch(env$select_de_features(
                   fids, n_de = 5, fc = 3,
                   class_id = rep("DNA", 50),
                   class_weights = c(LINE = 1, SINE = 0, LTR = 0, DNA = 0),
                   seed = 1L),
                 error = function(e) e)
report("select_de_features rejects all-zero effective weights",
       inherits(err2, "error"))

## run_de_gene_lib: gene-derived size factors give a topTags-shaped table
set.seed(11)
n_samp <- 8
samples_v <- paste0("s", seq_len(n_samp))
g_counts <- matrix(rpois(80 * n_samp, lambda = 200), nrow = 80,
                   dimnames = list(paste0("g", seq_len(80)), samples_v))
r_counts <- matrix(rpois(40 * n_samp, lambda = 60), nrow = 40,
                   dimnames = list(paste0("r", seq_len(40)), samples_v))
des <- model.matrix(~ factor(rep(c("A", "B"), each = n_samp / 2)))
gl_tt <- env$run_de_gene_lib(r_counts, g_counts, des)
report("run_de_gene_lib returns a topTags table aligned to repeat features",
       is.data.frame(gl_tt)
       && all(c("logFC", "PValue", "FDR") %in% colnames(gl_tt))
       && nrow(gl_tt) == nrow(r_counts))

## run_de_gene_lib: column-mismatched matrices trigger stopifnot
g_swap <- g_counts[, rev(seq_len(n_samp)), drop = FALSE]
err3 <- tryCatch(env$run_de_gene_lib(r_counts, g_swap, des),
                 error = function(e) e)
report("run_de_gene_lib rejects column-name mismatch",
       inherits(err3, "error"))

## recovered_logfc_table: aligns to feature_id space and NA for missing rows
fid <- paste0("r", seq_len(5))
planted_lfc <- c(log(3), 0, log(3), 0, 0)
planted_flag <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
fake_tt <- data.frame(logFC = c(1.1, -0.05, 0.9),
                      PValue = c(1e-3, 0.5, 1e-3),
                      FDR = c(1e-3, 0.5, 1e-3),
                      row.names = c("r1", "r2", "r3"))
rec <- env$recovered_logfc_table(fake_tt, fid, planted_lfc, planted_flag)
report("recovered_logfc_table aligns and NAs missing",
       nrow(rec) == 5
       && rec$feature_id[5] == "r5"
       && is.na(rec$recovered_logfc[4])
       && is.na(rec$recovered_logfc[5])
       && abs(rec$recovered_logfc[1] - 1.1) < 1e-9)

## TMM compresses recovered logFC under high-fraction same-sign DE
## while gene_lib stays close to the planted value. This is the central
## scientific claim of the coordinated-derepression panel; this test
## locks it in at the unit level.
set.seed(42)
n_g <- 600
n_r <- 200
gene_p <- list(mu = rep(c(40, 200, 800), length.out = n_g),
               dispersion = rep(0.15, n_g),
               feature_id = paste0("g", seq_len(n_g)))
repeat_p <- list(mu = rep(c(20, 80, 200), length.out = n_r),
                 dispersion = rep(0.25, n_r),
                 feature_id = paste0("r", seq_len(n_r)))
samples_e <- paste0("s", seq_len(8))
condition_e <- rep(c("A", "B"), each = 4)
co_out <- env$run_one_grid_cell(
  gene_p, repeat_p, samples_e, condition_e,
  gene_lib_scale = 1, repeat_lib_scale = 1,
  fc = 3, n_de_repeats = 100, n_de_genes = 100,
  sigma_w = 0.2, n_controls = 100, k_max = 1, fdr = 0.05,
  n_iter = 3, seed_base = 99,
  signed_fraction = 1.0,
  record_recovered_logfc = TRUE)

report("run_one_grid_cell returns scores + recovered tables",
       is.list(co_out)
       && all(c("scores", "recovered", "planted_features") %in% names(co_out))
       && is.data.frame(co_out$scores)
       && is.data.frame(co_out$recovered)
       && "gene_lib" %in% co_out$scores$method)

planted_logfc_value <- log(3)
mean_recovered_lfc <- function(rec_df, m) {
  d <- rec_df[rec_df$method == m & rec_df$planted, ]
  d <- d[!is.na(d$recovered_logfc), ]
  if (nrow(d) == 0) NA_real_ else mean(d$recovered_logfc)
}
mean_tmm <- mean_recovered_lfc(co_out$recovered, "TMM")
mean_gl <- mean_recovered_lfc(co_out$recovered, "gene_lib")
report("under coordinated derepression, TMM compresses logFC vs gene_lib",
       is.finite(mean_tmm) && is.finite(mean_gl)
       && mean_tmm < mean_gl
       && mean_gl > 0.5 * planted_logfc_value)

## load_repeat_class_map: empty path returns NULL
report("load_repeat_class_map empty path returns NULL",
       is.null(env$load_repeat_class_map("", paste0("r", seq_len(5)))))

## load_repeat_class_map: roundtrip
tmp_map <- tempfile(fileext = ".tsv")
write.table(data.frame(feature_id = c("r1", "r2", "r3"),
                       class_id = c("LINE", "SINE", "LTR")),
            tmp_map, sep = "\t", quote = FALSE, row.names = FALSE)
mm <- env$load_repeat_class_map(tmp_map, c("r1", "r2", "r3", "r99"))
report("load_repeat_class_map roundtrips and NA-fills missing",
       length(mm) == 4
       && mm[1] == "LINE" && mm[2] == "SINE" && mm[3] == "LTR"
       && is.na(mm[4]))
unlink(tmp_map)

## summarize_grid aggregates iterations correctly across all methods
sc_df <- co_out$scores
sc_df$gene_lib_scale <- 1
sc_df$repeat_lib_scale <- 1
sm <- env$summarize_grid(sc_df)
report("summarize_grid covers all method levels including gene_lib",
       nrow(sm) == length(unique(sc_df$method))
       && "gene_lib" %in% sm$method
       && all(c("power_mean", "power_sd", "fpr_mean", "fpr_sd") %in% colnames(sm)))

cat("OK: all de_simulations unit tests passed\n")
