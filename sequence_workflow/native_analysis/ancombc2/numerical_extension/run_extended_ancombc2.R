args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, args[3] %in% c("16S", "ITS"))
indir <- normalizePath(args[1]); outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(ANCOMBC))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(packageVersion("ANCOMBC")) == "2.14.0", getRversion() >= "4.5.0")
cfg <- fromJSON(file.path(indir, "FROZEN_CONFIG.json"), simplifyVector = FALSE)
stopifnot(cfg$formula == "locality + zone")
ids <- c("B1", "B2", "B3", "M1", "M2", "M3", "A1", "A2", "A3")
md <- read.delim(file.path(indir, "metadata_model.tsv"), check.names = FALSE,
                 stringsAsFactors = FALSE, row.names = 1)
stopifnot(identical(rownames(md), ids))
md$zone <- factor(md$zone, levels = c("Low degradation", "Moderate degradation", "High degradation"))
md$locality <- factor(md$locality, levels = c("J", "SI"))
x <- model.matrix(~ locality + zone, md)
stopifnot(qr(x)$rank == 4L, nrow(x) - qr(x)$rank == 5L)
write.table(cbind(sample_id = rownames(x), as.data.frame(x)), file.path(outdir, "MODEL_MATRIX.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "SESSION_INFO.txt"))
writeLines(capture.output(packageDescription("ANCOMBC")), file.path(outdir, "ANCOMBC_DESCRIPTION.txt"))
ns <- asNamespace("ANCOMBC")
for (fn in c("ancombc2", ".ancombc2_core", ".iter_mle", ".bias_em", ".ancombc_global_F")) {
  writeLines(capture.output(print(get(fn, ns))), file.path(outdir, paste0("SOURCE_NATIVE_", gsub("^\\.", "", fn), ".txt")))
}


# In-memory diagnostic instrumentation, numerical expressions unchanged.
.add_progress <- function(node) {
 if (!is.call(node)) return(node)
 if (identical(node[[1]], as.name("while"))) {
   stopifnot(identical(node[[3]][[1]], as.name("{")))
   node[[3]] <- as.call(c(as.list(node[[3]]), list(quote({
     if (iterNum %% 500L == 0L) message("ANCOM_AUDIT_PROGRESS marker=", .GlobalEnv$.audit_marker,
       " pseudo=", .GlobalEnv$.audit_pseudo, " iter=", iterNum,
       " epsilon=", format(epsilon, digits=16), " tol=", tol)
   }))))
   return(node)
 }
 for (j in seq_along(node)[-1]) if (is.call(node[[j]])) node[[j]] <- .add_progress(node[[j]])
 node
}
.bias_original <- get(".bias_em", ns)
.bias_instrumented <- .bias_original
body(.bias_instrumented) <- .add_progress(body(.bias_original))
unlockBinding(".bias_em", ns); assign(".bias_em", .bias_instrumented, envir=ns); lockBinding(".bias_em", ns)
writeLines(capture.output(print(.bias_instrumented)), file.path(outdir, "BIAS_EM_PROGRESS_ONLY_INSTRUMENTED.txt"))

# Diagnostics only: capture scalar stopping state in memory; mathematical
# function source and package files are not changed on disk.
.audit_marker <- NA_character_; .audit_pseudo <- NA_real_
.audit_path <- file.path(outdir, "CONVERGENCE_DIAGNOSTICS.jsonl")
.audit_counter <- 0L
.audit_exit <- function(name, frame) {
  val <- function(n) if (exists(n, envir = frame, inherits = FALSE)) get(n, envir = frame, inherits = FALSE) else NULL
  .GlobalEnv$.audit_counter <- .GlobalEnv$.audit_counter + 1L
  eps <- val("epsilon"); lim <- val("max_iter"); it <- val("iterNum"); tol <- val("tol")
  finite <- !is.null(eps) && length(eps) == 1L && is.finite(eps)
  rec <- list(event = .audit_counter, marker = .audit_marker, pseudo = .audit_pseudo,
              function_name = name, stage = val(".audit_stage"),
              iterations = it, max_iter = lim, epsilon = eps, tolerance = tol,
              converged = if (is.null(it)) NA else isTRUE(finite && eps <= tol),
              reached_limit = if (is.null(it)) NA else isTRUE(it >= lim),
              n_taxa = val("n_tax"), utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  cat(toJSON(rec, auto_unbox = TRUE, null = "null", na = "null", digits = NA), "\n",
      file = .audit_path, append = TRUE, sep = "")
}
trace(".ancombc2_core", where = ns,
      tracer = quote(.GlobalEnv$.audit_pseudo <- pseudo), print = FALSE)
trace(".iter_mle", where = ns,
      tracer = quote(.audit_stage <- if (is.null(theta)) "estimate_sampling_bias" else "given_sampling_bias"),
      exit = quote(.GlobalEnv$.audit_exit("iter_mle", environment())), print = FALSE)
trace(".bias_em", where = ns,
      exit = quote(.GlobalEnv$.audit_exit("bias_em", environment())), print = FALSE)
writeLines(capture.output(print(get(".iter_mle", ns))), file.path(outdir, "DIAGNOSTIC_TRACE_ITER_MLE.txt"))
writeLines(capture.output(print(get(".bias_em", ns))), file.path(outdir, "DIAGNOSTIC_TRACE_BIAS_EM.txt"))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
save_tsv <- function(obj, path) write.table(obj, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
all_status <- list()
for (marker in args[3]) {
  .audit_marker <- marker
  started <- Sys.time(); mp <- file.path(outdir, marker)
  dir.create(mp, showWarnings = FALSE)
  counts <- read.delim(file.path(indir, paste0(marker, "_counts_eligible.tsv")),
                        row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts); storage.mode(counts) <- "double"
  expected <- if (marker == "16S") 356L else 179L
  stopifnot(nrow(counts) == expected, identical(colnames(counts), ids),
            !anyDuplicated(rownames(counts)), all(is.finite(counts)),
            all(counts >= 0), all(counts == floor(counts)))
  for (i in seq_len(nrow(counts))) {
    pos <- counts[i, ] > 0
    stopifnot(sum(pos) >= 6L, sum(counts[i, ]) >= 100L,
              qr(x[pos, , drop = FALSE])$rank == 4L,
              sum(pos) - qr(x[pos, , drop = FALSE])$rank >= 2L)
  }
  warning_records <- character()
  set.seed(as.integer(cfg$seed))
  result <- tryCatch(withCallingHandlers(
    ANCOMBC::ancombc2(
      data = counts, taxa_are_rows = TRUE, meta_data = md, tax_level = NULL,
      fix_formula = "locality + zone", rand_formula = NULL,
      p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
      group = "zone", struc_zero = FALSE, neg_lb = FALSE,
      alpha = 0.05, n_cl = 1, verbose = TRUE,
      global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
      iter_control = list(tol = 0.01, max_iter = 100, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 10000)),
    warning = function(w) {
      warning_records <<- c(warning_records, conditionMessage(w))
      # Keep native warning emission; do not muffle or change the fit.
    }), error = function(e) e)
  writeLines(warning_records, file.path(mp, "WARNINGS.txt"))
  if (inherits(result, "error")) {
    status <- list(marker = marker, status = "FIT_FAILED", error = conditionMessage(result),
                    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    write_json(status, file.path(mp, "STATUS.json"), pretty = TRUE, auto_unbox = TRUE)
    all_status[[marker]] <- status
    next
  }
  saveRDS(result, file.path(mp, "ANCOMBC2_NATIVE_RESULT.rds"), compress = "xz")
  for (nm in c("res", "res_global", "ss_tab")) save_tsv(result[[nm]], file.path(mp, paste0("NATIVE_", nm, ".tsv")))
  save_tsv(data.frame(sample_id = names(result$samp_frac), sampling_fraction = as.numeric(result$samp_frac)),
           file.path(mp, "SAMPLING_FRACTIONS.tsv"))
  save_tsv(data.frame(parameter = colnames(x), delta_em = as.numeric(result$delta_em),
                      delta_wls = as.numeric(result$delta_wls)), file.path(mp, "BIAS_ESTIMATES.tsv"))
  stopifnot(nrow(result$res) == expected, nrow(result$res_global) == expected,
            setequal(result$res$taxon, rownames(counts)),
            isTRUE(all.equal(unname(as.matrix(result$feature_table)), unname(counts))))
  diag_connection <- file(.audit_path, "r")
  diag <- stream_in(diag_connection, verbose = FALSE)
  close(diag_connection)
  dm <- diag[diag$marker == marker, , drop = FALSE]
  iterative <- !is.na(dm$iterations)
  convergence_ok <- nrow(dm) >= 24L && isTRUE(all(dm$converged[iterative]))
  sampling_ok <- all(is.finite(result$samp_frac)) && all(is.finite(result$delta_em)) && all(is.finite(result$delta_wls))
  joined <- merge(result$res, result$res_global, by = "taxon", sort = FALSE,
                    suffixes = c("_coefficient", "_global"))
  joined <- joined[match(rownames(counts), joined$taxon), , drop = FALSE]
  id <- read_tsv(file.path(indir, paste0(marker, "_feature_identity.tsv")))
  id$nominal_genus <- tolower(as.character(id$nominal_genus)) == "true"
  names(id)[1] <- "taxon"
  joined <- merge(id, joined, by = "taxon", sort = FALSE)
  joined <- joined[match(rownames(counts), joined$taxon), , drop = FALSE]
  reg <- read_tsv(file.path(indir, paste0(marker, "_support_registry.tsv")))
  rr <- reg[match(joined$taxon, reg$feature_id), ]
  target_cols <- grep("^(lfc|se)_zone", names(joined), value = TRUE)
  stopifnot(length(target_cols) == 4L)
  finite_coef <- apply(joined[, target_cols, drop = FALSE], 1, function(z) all(is.finite(z)))
  positive_se <- apply(joined[, grep("^se_zone", names(joined)), drop = FALSE], 1, function(z) all(z > 0))
  joined$valid_numeric_result <- finite_coef & positive_se & is.finite(joined$W) &
    is.finite(joined$p_val) & is.finite(joined$q_val) & sampling_ok & convergence_ok
  joined$failure_or_uncertainty <- ifelse(joined$valid_numeric_result, "", "NONFINITE_OR_NONPOSITIVE_SE_OR_BIAS_OR_UNVERIFIED_CONVERGENCE")
  joined$possible_native_failure_encoded_as_one <- !finite_coef | !is.finite(joined$W)
  joined$robust_exploratory_candidate <- joined$valid_numeric_result & !is.na(joined$passed_ss) &
    joined$passed_ss & !is.na(joined$diff_robust_abn) & joined$diff_robust_abn & joined$q_val <= 0.05
  joined$marker <- marker; joined$family_size <- expected
  joined$total_count <- rr$total_count; joined$positive_samples <- rr$positive_samples
  joined$positive_design_rank <- rr$positive_design_rank; joined$positive_residual_df <- rr$positive_residual_df
  joined$mean_full_library_fraction <- rr$mean_full_library_fraction
  save_tsv(joined, file.path(mp, "ALL_ANNOTATION_GROUP_RESULTS.tsv"))
  save_tsv(joined[joined$nominal_genus, , drop = FALSE], file.path(mp, "NOMINAL_GENUS_RESULTS_FULL_FAMILY_Q.tsv"))
  save_tsv(joined[joined$robust_exploratory_candidate, , drop = FALSE], file.path(mp, "ROBUST_EXPLORATORY_CANDIDATES.tsv"))
  status <- list(marker = marker, status = if (all(joined$valid_numeric_result)) "EXECUTED_NUMERICALLY_VALID_EXPLORATORY" else "EXECUTED_WITH_NUMERIC_FLAGS",
                 hypotheses = expected, native_global_q_le_005 = sum(joined$q_val <= 0.05, na.rm = TRUE),
                 native_global_robust = sum(joined$diff_robust_abn, na.rm = TRUE),
                 valid_robust_candidates = sum(joined$robust_exploratory_candidate),
                 valid_nominal_robust_candidates = sum(joined$robust_exploratory_candidate & joined$nominal_genus),
                 invalid_rows = sum(!joined$valid_numeric_result),
                 diagnostic_events = nrow(dm), convergence_verified = convergence_ok,
                 finite_bias_and_sampling_fractions = sampling_ok,
                 seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                 limitation = "n=3/category; small-sample calibration, field independence and ecological confounding remain unresolved; no causal or measured absolute-abundance conclusion.")
  write_json(status, file.path(mp, "STATUS.json"), pretty = TRUE, auto_unbox = TRUE)
  all_status[[marker]] <- status
}
for (fn in c(".ancombc2_core", ".iter_mle", ".bias_em")) untrace(fn, where = ns)
write_json(all_status, file.path(outdir, "ANCOMBC2_EXECUTION_SUMMARY.json"), pretty = TRUE, auto_unbox = TRUE)
if (any(vapply(all_status, function(s) identical(s$status, "FIT_FAILED"), logical(1)))) quit(status = 2L)
