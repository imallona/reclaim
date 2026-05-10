#!/usr/bin/env Rscript
## de_simulations.R
##
## Count-level power benchmark for repeat differential expression as a function
## of gene library size and repeat library size. Reads two pre-computed count
## matrices (rows = features, cols = samples), fits negative-binomial mean and
## dispersion from each, then sweeps a 2D grid of multiplicative library-size
## scalers. At each grid cell it simulates new counts with a known fold change
## planted at a fixed set of repeat features, runs several normalization
## strategies, and reports power and false positive rate over n_iter replicates.
##
## Two scenarios are supported via config keys:
##   - default: a small fraction of repeats DE, mixed direction. Probes
##     power and FPR; the "most features not DE" assumption underlying TMM
##     holds and TMM should be unbiased.
##   - coordinated derepression: a large fraction of repeats DE, all in the
##     same direction (signed_fraction = 1). Mirrors the TDP-43 OE biology
##     where many young retrotransposon subfamilies move up together. The
##     "most features not DE" assumption breaks; TMM-on-repeats compresses
##     the recovered logFC. Gene-library-size transfer fits norm factors on
##     the gene matrix (where the assumption still holds) and is unbiased.
##
## Both scenarios share the same code path. Method coverage in both:
##   none, TMM, RUVg_k1..k_max, gene_lib.
##
## Outputs: power/FPR per method per grid cell (de_simulations_results.tsv),
## per-iteration summary (de_simulations_summary.tsv), and per-feature
## recovered-vs-planted logFC for bias-slope estimation
## (de_simulations_recovered_logfc.tsv, gzipped).
##
## See docs/de_simulations.md for the plain-English description.

suppressPackageStartupMessages({
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(Biobase)
  library(matrixStats)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(viridis)
})


load_count_matrix <- function(path, samples = NULL) {
  stopifnot(file.exists(path))
  m <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE)
  m <- as.matrix(m)
  if (!is.null(samples)) {
    stopifnot(all(samples %in% colnames(m)))
    m <- m[, samples, drop = FALSE]
  }
  storage.mode(m) <- "integer"
  stopifnot(all(is.finite(m)), all(m >= 0))
  m
}


fit_nb_params <- function(counts, design) {
  stopifnot(is.matrix(counts), ncol(counts) == nrow(design))
  min_present <- max(2L, ncol(counts) %/% 4L)
  keep <- rowSums(counts >= 1) >= min_present
  counts <- counts[keep, , drop = FALSE]
  stopifnot(nrow(counts) >= 2)
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  mu <- pmax(rowMeans(counts), 1e-6)
  disp <- dge$tagwise.dispersion
  if (is.null(disp) || length(disp) != length(mu)) {
    disp <- rep_len(dge$common.dispersion, length(mu))
  }
  list(mu = mu, dispersion = disp, feature_id = rownames(counts))
}


simulate_count_matrix <- function(mu, dispersion, n_samples, lib_scale = 1,
                                  w_loadings = NULL, w_factor = NULL,
                                  de_logfc = NULL, condition = NULL,
                                  seed = 1L) {
  set.seed(seed)
  n_feat <- length(mu)
  stopifnot(length(dispersion) == n_feat, n_samples >= 2)
  log_mu <- matrix(log(mu * lib_scale), nrow = n_feat, ncol = n_samples)
  if (!is.null(w_loadings) && !is.null(w_factor)) {
    stopifnot(length(w_loadings) == n_feat, length(w_factor) == n_samples)
    log_mu <- log_mu + outer(w_loadings, w_factor)
  }
  if (!is.null(de_logfc) && !is.null(condition)) {
    stopifnot(is.factor(condition), length(condition) == n_samples,
              length(de_logfc) == n_feat, length(levels(condition)) == 2)
    treated <- as.numeric(condition == levels(condition)[2])
    log_mu <- log_mu + outer(de_logfc, treated)
  }
  size <- 1 / pmax(dispersion, 1e-4)
  mu_mat <- exp(log_mu)
  counts <- matrix(rnbinom(n_feat * n_samples,
                           size = rep(size, n_samples),
                           mu = as.vector(mu_mat)),
                   nrow = n_feat)
  storage.mode(counts) <- "integer"
  counts
}


## Plant DE assignments at the feature level. Returns the picked indices, a
## per-feature logFC vector matching the full feature space, and a logical
## flag indicating which features carry signal. The signed_fraction argument
## controls how many of the picked features get a positive logFC (vs negative).
## A signed_fraction of 1.0 means full coordinated derepression (all up); 0.5
## means balanced; values in between scale linearly. class_weights, if set,
## must be a named numeric vector with one entry per class_id; selection
## probability is proportional to class weight.
##
## Asserts on bounds, length matching, and weight positivity guard against
## silent miswiring from upstream YAML.
select_de_features <- function(feature_id, n_de, fc, signed_fraction = 1.0,
                               class_id = NULL, class_weights = NULL,
                               seed = 1L) {
  stopifnot(is.character(feature_id), length(feature_id) >= 2L,
            !anyDuplicated(feature_id),
            is.numeric(n_de), length(n_de) == 1L, n_de >= 0,
            is.numeric(fc), length(fc) == 1L, fc > 0,
            is.numeric(signed_fraction), length(signed_fraction) == 1L,
            signed_fraction >= 0, signed_fraction <= 1)
  n_feat <- length(feature_id)
  n_de <- as.integer(min(n_de, n_feat))

  if (!is.null(class_id) || !is.null(class_weights)) {
    stopifnot(!is.null(class_id), !is.null(class_weights),
              length(class_id) == n_feat,
              !is.null(names(class_weights)),
              all(class_weights >= 0), any(class_weights > 0))
    weights <- class_weights[class_id]
    weights[is.na(weights)] <- 0
    stopifnot(any(weights > 0))
  } else {
    weights <- rep(1, n_feat)
  }

  set.seed(seed)
  if (n_de == 0L) {
    return(list(idx = integer(0),
                logfc = numeric(n_feat),
                planted = rep(FALSE, n_feat)))
  }
  idx <- sample.int(n_feat, n_de, prob = weights / sum(weights))

  n_up <- as.integer(round(signed_fraction * n_de))
  n_dn <- n_de - n_up
  signs <- sample(c(rep(1, n_up), rep(-1, n_dn)))
  stopifnot(length(signs) == n_de)

  logfc <- numeric(n_feat)
  logfc[idx] <- signs * log(fc)

  planted <- rep(FALSE, n_feat)
  planted[idx] <- TRUE

  list(idx = idx, logfc = logfc, planted = planted)
}


## Resolve the per-feature class_id vector for the repeat matrix from an
## optional class_map TSV. Returns NULL when class_map_path is empty.
## Two TSV formats are accepted:
##   - 2-column: feature_id, class_id (used as-is).
##   - locus_map (4-column with transcript_id, gene_id, family_id, class_id):
##     deduplicated on (granularity_column, class_id) and rolled up to a
##     feature_id->class_id map. The granularity_column argument selects
##     gene_id, family_id, or transcript_id, matching whatever level the
##     repeat count matrix uses.
## Missing features become NA class. The granularity_column argument is
## ignored when the input is already a 2-column class_map. Kept in the
## pure-R section so the helper is testable without Bioconductor.
load_repeat_class_map <- function(class_map_path, feature_id,
                                  granularity_column = "family_id") {
  if (is.null(class_map_path) || !nzchar(class_map_path)) {
    return(NULL)
  }
  stopifnot(file.exists(class_map_path),
            granularity_column %in% c("gene_id", "family_id", "transcript_id"))
  m <- read.table(class_map_path, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("class_id" %in% colnames(m))
  if ("feature_id" %in% colnames(m)) {
    key_col <- "feature_id"
  } else {
    stopifnot(granularity_column %in% colnames(m))
    key_col <- granularity_column
  }
  ## Deduplicate on (key_col, class_id) and warn if the same key has more
  ## than one class. The pipeline's locus_map should not produce that, so
  ## the warning here is a tripwire.
  m <- unique(m[, c(key_col, "class_id"), drop = FALSE])
  dup_keys <- m[[key_col]][duplicated(m[[key_col]])]
  if (length(dup_keys) > 0) {
    warning(sprintf(
      "load_repeat_class_map: %d feature(s) map to >1 class; first wins (e.g. %s)",
      length(unique(dup_keys)), paste(head(dup_keys, 3), collapse = ", ")))
    m <- m[!duplicated(m[[key_col]]), , drop = FALSE]
  }
  out <- m$class_id[match(feature_id, m[[key_col]])]
  stopifnot(length(out) == length(feature_id))
  out
}


## edgeR DE on a repeat count matrix using TMM size factors fit on the gene
## matrix. Sample ordering must match between the two matrices; an assert
## guards against silent reordering from caller code. No RUV covariates;
## the design is whatever the caller passes.
run_de_gene_lib <- function(repeat_counts, gene_counts, design) {
  stopifnot(is.matrix(repeat_counts), is.matrix(gene_counts),
            ncol(repeat_counts) == ncol(gene_counts),
            identical(colnames(repeat_counts), colnames(gene_counts)),
            qr(design)$rank == ncol(design),
            nrow(design) == ncol(repeat_counts))
  gene_dge <- DGEList(counts = gene_counts)
  gene_dge <- calcNormFactors(gene_dge, method = "TMM")
  dge <- DGEList(counts = repeat_counts,
                 lib.size = gene_dge$samples$lib.size,
                 norm.factors = gene_dge$samples$norm.factors)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  topTags(qlf, n = Inf, sort.by = "none")$table
}


## Build a recovered-logFC table aligned to the full feature_id space.
## Features filtered out before testing get NA recovered_logfc; the planted
## flag and planted_logfc are still populated so the caller can compute bias
## slope conditional on the feature being testable. Used for the
## TMM-compression bias panel.
recovered_logfc_table <- function(tt_df, feature_id, planted_logfc, planted) {
  stopifnot(is.character(feature_id),
            length(feature_id) == length(planted_logfc),
            length(feature_id) == length(planted),
            is.numeric(planted_logfc), is.logical(planted))
  recovered <- rep(NA_real_, length(feature_id))
  if (!is.null(tt_df) && nrow(tt_df) > 0) {
    m <- match(feature_id, rownames(tt_df))
    recovered[!is.na(m)] <- tt_df$logFC[m[!is.na(m)]]
  }
  data.frame(feature_id = feature_id,
             planted = planted,
             planted_logfc = planted_logfc,
             recovered_logfc = recovered,
             stringsAsFactors = FALSE)
}


pick_empirical_controls <- function(counts, design, n_controls) {
  stopifnot(qr(design)$rank == ncol(design),
            n_controls >= 1, n_controls <= nrow(counts))
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  rownames(tt)[order(tt$PValue, decreasing = TRUE)][seq_len(n_controls)]
}


run_de <- function(counts, design) {
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  topTags(qlf, n = Inf, sort.by = "none")$table
}


run_de_no_norm <- function(counts, design) {
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- rep(1, ncol(counts))
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  topTags(qlf, n = Inf, sort.by = "none")$table
}


## Power and FPR scoring. planted_features is the FULL planted-truth set,
## including features that might have been filtered out before testing.
## Filtered planted features count as FN, which is the correct behavior when
## comparing across grid cells with very different repeat library sizes.
score_calls <- function(tt_df, planted_features, fdr = 0.05) {
  called <- rownames(tt_df)[!is.na(tt_df$FDR) & tt_df$FDR <= fdr]
  truth <- planted_features
  null <- setdiff(rownames(tt_df), truth)
  tp <- sum(called %in% truth)
  fn <- length(truth) - tp
  fp <- sum(called %in% null)
  tn <- length(null) - fp
  list(power = tp / max(1L, tp + fn),
       fpr = fp / max(1L, fp + tn),
       tp = tp, fn = fn, fp = fp, tn = tn)
}


## Run all methods on a (filtered) repeat matrix and return a list of
## per-method topTags tables. RUVg keys are "RUVg_k1".."RUVg_kN" by k.
## Returns NULL entries for methods that could not run (e.g. rank-deficient
## RUV designs at small sample sizes); callers must handle NULL.
run_methods_on_replicate <- function(repeat_counts_f, gene_counts_f, design,
                                     k_max, n_controls) {
  stopifnot(is.matrix(repeat_counts_f), is.matrix(gene_counts_f),
            ncol(repeat_counts_f) == ncol(gene_counts_f),
            identical(colnames(repeat_counts_f), colnames(gene_counts_f)))
  out <- list()
  out[["none"]] <- if (nrow(repeat_counts_f) >= 2)
    run_de_no_norm(repeat_counts_f, design) else NULL
  out[["TMM"]] <- if (nrow(repeat_counts_f) >= 2)
    run_de(repeat_counts_f, design) else NULL
  out[["gene_lib"]] <- if (nrow(repeat_counts_f) >= 2 && nrow(gene_counts_f) >= 2)
    tryCatch(run_de_gene_lib(repeat_counts_f, gene_counts_f, design),
             error = function(e) NULL) else NULL

  n_ctrl <- min(n_controls, max(0L, nrow(gene_counts_f) - 1L))
  if (n_ctrl >= 1 && nrow(gene_counts_f) >= 2) {
    controls <- pick_empirical_controls(gene_counts_f, design, n_ctrl)
    pheno <- AnnotatedDataFrame(data.frame(
      row.names = colnames(gene_counts_f)))
    set1 <- newSeqExpressionSet(counts = gene_counts_f, phenoData = pheno)
    for (k in seq_len(k_max)) {
      key <- paste0("RUVg_k", k)
      out[[key]] <- NULL
      ruv_fit <- tryCatch(RUVg(set1, controls, k = k), error = function(e) NULL)
      if (!is.null(ruv_fit)) {
        W <- as.matrix(pData(ruv_fit)[, grep("^W_", colnames(pData(ruv_fit))),
                                      drop = FALSE])
        design_ruv <- cbind(design[, -ncol(design), drop = FALSE], W,
                            design[, ncol(design), drop = FALSE])
        if (qr(design_ruv)$rank == ncol(design_ruv) && nrow(repeat_counts_f) >= 2) {
          out[[key]] <- tryCatch(run_de(repeat_counts_f, design_ruv),
                                 error = function(e) NULL)
        }
      }
    }
  } else {
    for (k in seq_len(k_max)) {
      out[[paste0("RUVg_k", k)]] <- NULL
    }
  }
  out
}


run_one_grid_cell <- function(gene_params, repeat_params, samples, condition,
                              gene_lib_scale, repeat_lib_scale,
                              fc, n_de_repeats, n_de_genes,
                              sigma_w, n_controls, k_max, fdr,
                              n_iter, seed_base,
                              signed_fraction = 1.0,
                              repeat_class_id = NULL,
                              repeat_class_weights = NULL,
                              record_recovered_logfc = TRUE) {
  n_samples <- length(samples)
  stopifnot(length(condition) == n_samples)
  cond_factor <- factor(condition, levels = sort(unique(condition)))
  stopifnot(length(levels(cond_factor)) == 2)
  design <- model.matrix(~ cond_factor)

  set.seed(seed_base)
  w_factor <- rnorm(n_samples)
  gene_w_loadings <- rnorm(length(gene_params$mu), 0, sigma_w)
  repeat_w_loadings <- rnorm(length(repeat_params$mu), 0, sigma_w)

  repeat_plant <- select_de_features(
    feature_id = repeat_params$feature_id,
    n_de = n_de_repeats, fc = fc,
    signed_fraction = signed_fraction,
    class_id = repeat_class_id,
    class_weights = repeat_class_weights,
    seed = seed_base + 1L)
  de_repeat_features <- repeat_params$feature_id[repeat_plant$idx]
  repeat_logfc <- repeat_plant$logfc

  gene_plant <- select_de_features(
    feature_id = gene_params$feature_id,
    n_de = n_de_genes, fc = 2,
    signed_fraction = 0.5,
    seed = seed_base + 2L)
  gene_logfc <- gene_plant$logfc

  method_keys <- c("none", "TMM", "gene_lib", paste0("RUVg_k", seq_len(k_max)))

  iter_results <- list()
  recovered_results <- list()
  for (i in seq_len(n_iter)) {
    seed_i <- seed_base + 100L + i
    gene_counts <- simulate_count_matrix(
      mu = gene_params$mu, dispersion = gene_params$dispersion,
      n_samples = n_samples, lib_scale = gene_lib_scale,
      w_loadings = gene_w_loadings, w_factor = w_factor,
      de_logfc = gene_logfc, condition = cond_factor, seed = seed_i)
    rownames(gene_counts) <- gene_params$feature_id
    colnames(gene_counts) <- samples

    repeat_counts <- simulate_count_matrix(
      mu = repeat_params$mu, dispersion = repeat_params$dispersion,
      n_samples = n_samples, lib_scale = repeat_lib_scale,
      w_loadings = repeat_w_loadings, w_factor = w_factor,
      de_logfc = repeat_logfc, condition = cond_factor, seed = seed_i + 1L)
    rownames(repeat_counts) <- repeat_params$feature_id
    colnames(repeat_counts) <- samples

    keep_g <- rowSums(gene_counts >= 1) >= 2L
    gene_counts_f <- gene_counts[keep_g, , drop = FALSE]
    keep_r <- rowSums(repeat_counts >= 1) >= 2L
    repeat_counts_f <- repeat_counts[keep_r, , drop = FALSE]

    tt_per_method <- run_methods_on_replicate(
      repeat_counts_f, gene_counts_f, design, k_max, n_controls)

    method_scores <- list()
    method_recovered <- list()
    for (m in method_keys) {
      tt <- tt_per_method[[m]]
      if (is.null(tt)) {
        method_scores[[m]] <- list(power = NA_real_, fpr = NA_real_,
                                   tp = NA_integer_, fn = NA_integer_,
                                   fp = NA_integer_, tn = NA_integer_)
        if (record_recovered_logfc) {
          method_recovered[[m]] <- recovered_logfc_table(
            NULL, repeat_params$feature_id, repeat_logfc, repeat_plant$planted)
        }
      } else {
        method_scores[[m]] <- score_calls(tt, de_repeat_features, fdr)
        if (record_recovered_logfc) {
          method_recovered[[m]] <- recovered_logfc_table(
            tt, repeat_params$feature_id, repeat_logfc, repeat_plant$planted)
        }
      }
    }

    iter_df <- do.call(rbind, lapply(method_keys, function(m) {
      r <- method_scores[[m]]
      data.frame(method = m, iter = i,
                 power = r$power, fpr = r$fpr,
                 tp = r$tp, fn = r$fn, fp = r$fp, tn = r$tn,
                 stringsAsFactors = FALSE)
    }))
    iter_results[[i]] <- iter_df

    if (record_recovered_logfc) {
      rec_df <- do.call(rbind, lapply(method_keys, function(m) {
        d <- method_recovered[[m]]
        d$method <- m
        d$iter <- i
        d
      }))
      recovered_results[[i]] <- rec_df
    }
  }
  list(scores = do.call(rbind, iter_results),
       recovered = if (record_recovered_logfc) do.call(rbind, recovered_results)
                   else NULL,
       planted_features = de_repeat_features,
       planted_logfc = repeat_logfc[repeat_plant$idx])
}


run_grid <- function(gene_counts_path, repeat_counts_path, metadata,
                     gene_lib_grid, repeat_lib_grid, fc,
                     n_de_repeats, n_de_genes, sigma_w,
                     n_controls, k_max, fdr, n_iter, seed,
                     signed_fraction = 1.0,
                     class_map_path = NULL,
                     class_map_granularity = "family_id",
                     class_weights = NULL,
                     record_recovered_logfc = TRUE) {
  samples <- metadata$sample
  condition <- metadata$condition
  gene_mat <- load_count_matrix(gene_counts_path, samples)
  repeat_mat <- load_count_matrix(repeat_counts_path, samples)
  cond_factor <- factor(condition, levels = sort(unique(condition)))
  design <- model.matrix(~ cond_factor)
  gene_params <- fit_nb_params(gene_mat, design)
  repeat_params <- fit_nb_params(repeat_mat, design)

  repeat_class_id <- load_repeat_class_map(class_map_path,
                                           repeat_params$feature_id,
                                           granularity_column = class_map_granularity)
  if (!is.null(class_weights) && is.null(repeat_class_id)) {
    stop("class_weights set but class_map_path missing or empty")
  }

  cells <- expand.grid(gene_lib_scale = gene_lib_grid,
                       repeat_lib_scale = repeat_lib_grid,
                       KEEP.OUT.ATTRS = FALSE)
  results <- vector("list", nrow(cells))
  recovered <- vector("list", nrow(cells))
  for (ci in seq_len(nrow(cells))) {
    g <- cells$gene_lib_scale[ci]
    r <- cells$repeat_lib_scale[ci]
    cell_seed <- as.integer(seed) * 1000L + ci * 10L
    cat(sprintf("[grid %d/%d] gene_lib_scale=%g repeat_lib_scale=%g\n",
                ci, nrow(cells), g, r))
    cell_out <- run_one_grid_cell(
      gene_params, repeat_params, samples, condition,
      gene_lib_scale = g, repeat_lib_scale = r,
      fc = fc, n_de_repeats = n_de_repeats, n_de_genes = n_de_genes,
      sigma_w = sigma_w, n_controls = n_controls, k_max = k_max,
      fdr = fdr, n_iter = n_iter, seed_base = cell_seed,
      signed_fraction = signed_fraction,
      repeat_class_id = repeat_class_id,
      repeat_class_weights = class_weights,
      record_recovered_logfc = record_recovered_logfc)
    df <- cell_out$scores
    df$gene_lib_scale <- g
    df$repeat_lib_scale <- r
    results[[ci]] <- df
    if (record_recovered_logfc && !is.null(cell_out$recovered)) {
      rec <- cell_out$recovered
      rec$gene_lib_scale <- g
      rec$repeat_lib_scale <- r
      if (!is.null(repeat_class_id)) {
        rec$class_id <- repeat_class_id[match(rec$feature_id,
                                              repeat_params$feature_id)]
      }
      recovered[[ci]] <- rec
    }
  }
  list(results = do.call(rbind, results),
       recovered = if (record_recovered_logfc) do.call(rbind, recovered)
                   else NULL,
       gene_params = gene_params,
       repeat_params = repeat_params,
       repeat_class_id = repeat_class_id)
}


summarize_grid <- function(results_df) {
  results_df %>%
    group_by(method, gene_lib_scale, repeat_lib_scale) %>%
    summarise(power_mean = mean(power, na.rm = TRUE),
              power_sd = sd(power, na.rm = TRUE),
              fpr_mean = mean(fpr, na.rm = TRUE),
              fpr_sd = sd(fpr, na.rm = TRUE),
              n_iter = n(),
              .groups = "drop")
}


plot_heatmap <- function(summary_df, metric = c("power_mean", "fpr_mean"),
                         title = NULL) {
  metric <- match.arg(metric)
  ggplot(summary_df, aes(x = factor(gene_lib_scale),
                         y = factor(repeat_lib_scale),
                         fill = .data[[metric]])) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", .data[[metric]])), size = 2.5) +
    facet_wrap(~ method) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = "gene library size scaler",
         y = "repeat library size scaler",
         fill = metric, title = title) +
    theme_minimal()
}


# Snakemake glue. Only runs when invoked via the snakemake script: directive.
if (exists("snakemake")) {
  log_path <- snakemake@log[[1]]
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(log_path, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  params <- snakemake@params
  outdir <- params$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  gene_counts_path <- snakemake@input$gene_counts
  repeat_counts_path <- snakemake@input$repeat_counts

  if (nzchar(params$metadata)) {
    md_raw <- read.table(params$metadata, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
    metadata <- data.frame(
      sample = md_raw[[params$sample_column]],
      condition = md_raw[[params$condition_column]],
      stringsAsFactors = FALSE)
  } else {
    header <- read.table(gene_counts_path, header = TRUE, sep = "\t",
                         nrows = 1, check.names = FALSE)
    samples <- colnames(header)[-1]
    half <- length(samples) %/% 2
    metadata <- data.frame(
      sample = samples,
      condition = c(rep("A", half), rep("B", length(samples) - half)),
      stringsAsFactors = FALSE)
  }

  signed_fraction <- if (!is.null(params$signed_fraction))
    as.numeric(params$signed_fraction) else 1.0
  stopifnot(signed_fraction >= 0, signed_fraction <= 1)

  class_map_path <- if (!is.null(params$class_map_path))
    as.character(params$class_map_path) else ""
  class_weights <- NULL
  if (!is.null(params$class_weights) && length(params$class_weights) > 0) {
    cw <- params$class_weights
    if (is.list(cw)) cw <- unlist(cw)
    stopifnot(is.numeric(cw), !is.null(names(cw)),
              all(cw >= 0), any(cw > 0))
    class_weights <- cw
  }

  grid_out <- run_grid(
    gene_counts_path = gene_counts_path,
    repeat_counts_path = repeat_counts_path,
    metadata = metadata,
    gene_lib_grid = as.numeric(unlist(params$gene_lib_grid)),
    repeat_lib_grid = as.numeric(unlist(params$repeat_lib_grid)),
    fc = as.numeric(params$fc),
    n_de_repeats = as.integer(params$n_de_repeats),
    n_de_genes = as.integer(params$n_de_genes),
    sigma_w = as.numeric(params$sigma_w),
    n_controls = as.integer(params$n_controls),
    k_max = as.integer(params$k_max),
    fdr = as.numeric(params$fdr),
    n_iter = as.integer(params$n_iter),
    seed = as.integer(params$seed),
    signed_fraction = signed_fraction,
    class_map_path = class_map_path,
    class_map_granularity = if (!is.null(params$class_map_granularity))
      as.character(params$class_map_granularity) else "family_id",
    class_weights = class_weights,
    record_recovered_logfc = TRUE)

  results <- grid_out$results
  summary_df <- summarize_grid(results)

  write.table(results, snakemake@output$results_tsv,
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_df, snakemake@output$summary_tsv,
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (!is.null(snakemake@output$recovered_logfc_tsv)) {
    rec <- grid_out$recovered
    stopifnot(!is.null(rec))
    rec_path <- snakemake@output$recovered_logfc_tsv
    if (grepl("\\.gz$", rec_path)) {
      gz <- gzfile(rec_path, "w")
      on.exit(close(gz), add = TRUE)
      write.table(rec, gz, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      write.table(rec, rec_path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  saveRDS(list(results = results, summary = summary_df,
               metadata = metadata,
               gene_params = grid_out$gene_params,
               repeat_params = grid_out$repeat_params,
               recovered = grid_out$recovered,
               repeat_class_id = grid_out$repeat_class_id,
               params = as.list(params)),
          snakemake@output$results_rds)

  pdf(snakemake@output$heatmap_pdf, width = 11, height = 7)
  print(plot_heatmap(summary_df, "power_mean",
                     sprintf("Power at FC=%g, FDR=%g",
                             as.numeric(params$fc), as.numeric(params$fdr))))
  print(plot_heatmap(summary_df, "fpr_mean",
                     sprintf("False positive rate at FDR=%g",
                             as.numeric(params$fdr))))
  dev.off()

  sink(type = "message")
  sink(type = "output")
  close(log_con)
}
