# ===========================================================================
# PAPER 2 - STEP 6b v2 summary: reads results/step6b_v2/*.rds written by R/29
# (copy of R/26 with: job names taken from the files, a STUCK check -
# fits whose Kss moved < 1% from its start - and an optional folder argument).
# Safe to run at any time, also while R/25 is still running (uses finished
# fits only). Pairs A_M3 / D_M3 / D_M1 by replicate.
#
#   Per job and parameter:
#     n         - successful fits (error-free)
#     bias%     - median relative error, 100*(est/true - 1)
#     rRMSE%    - 100*sqrt(mean((est/true - 1)^2))
#     empSD%    - SD of log(est) x 100 (empirical precision)
#     meanRSE%  - mean reported SE(log est) x 100 (model-based precision)
#   Also: covariance-step success, % BLQ, fit time.
#   KSS head-to-head on replicates where all three jobs succeeded.
#   PopED (FO, CRAN 0.7.0) expected Kss RSE: A 78.9%, D 33.6% (R/21).
#   Writes no files. Run from the project root.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args)) args[1] else "results/step6b_v2"
f <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
if (!length(f)) stop("no results yet")
R <- lapply(f, readRDS)
ok <- !vapply(R, function(r) !is.null(r$error), logical(1))
cat(sprintf("%d fits found, %d failed\n", length(R), sum(!ok)))
if (any(!ok)) for (r in R[!ok]) cat("  FAILED", r$job, "rep", r$rep, ":", r$error, "\n")
R <- R[ok]
truth <- R[[1]]$truth; pn <- names(R[[1]]$est)
cat("versions:", paste(names(R[[1]]$versions), R[[1]]$versions, collapse = ", "), "\n\n")

job_names <- sort(unique(vapply(R, function(r) r$job, "")))
for (j in job_names) {
  Rj <- Filter(function(r) r$job == j, R)
  if (!length(Rj)) next
  E <- t(vapply(Rj, function(r) unname(r$est / r$truth), numeric(length(pn))))
  S <- t(vapply(Rj, function(r) unname(r$se_log), numeric(length(pn))))
  covok <- vapply(Rj, function(r) !is.null(r$cov) && nzchar(r$cov) && !grepl("fail", r$cov, ignore.case = TRUE), logical(1))
  stuck <- vapply(Rj, function(r) !is.null(r$start) && abs(log(r$est[["lkss"]] / r$start[["KSS"]])) < 0.01, logical(1))
  cat(sprintf("=== %s: n = %d | covariance ok %d | Kss STUCK at start %d | BLQ %.1f%% | time median %.0f min ===\n",
              j, length(Rj), sum(covok), sum(stuck),
              100 * mean(vapply(Rj, function(r) r$n_blq / r$n_obs, numeric(1))),
              median(vapply(Rj, function(r) r$time_s, numeric(1))) / 60))
  print(data.frame(`bias%` = round(100 * apply(E - 1, 2, median), 1),
                   `rRMSE%` = round(100 * sqrt(colMeans((E - 1)^2)), 1),
                   `empSD%` = if (nrow(E) > 1) round(100 * apply(log(E), 2, sd), 1) else NA,
                   `meanRSE%` = round(100 * colMeans(S, na.rm = TRUE), 1),
                   row.names = toupper(sub("^l", "", pn)), check.names = FALSE))
  cat("\n")
}

# paired Kss
key <- function(r) paste(r$job, r$rep)
kss <- setNames(vapply(R, function(r) r$est[["lkss"]], numeric(1)), vapply(R, key, ""))
reps <- sort(unique(vapply(R, function(r) r$rep, numeric(1))))
P <- t(vapply(reps, function(k) kss[paste(job_names, k)], numeric(length(job_names))))
P <- P[stats::complete.cases(P), , drop = FALSE]
if (nrow(P)) {
  colnames(P) <- job_names
  cat(sprintf("=== Kss, paired replicates (n = %d; truth %.2f) ===\n", nrow(P), truth[["KSS"]]))
  print(round(apply(P, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)), 3))
  cat("Share within +/-30% of truth:", paste(sprintf("%s %.0f%%", job_names,
      100 * colMeans(abs(P / truth[["KSS"]] - 1) <= 0.3)), collapse = " | "), "\n")
}
cat("\nPopED (FO) expected Kss RSE: A 78.9%, D 33.6%. Compare with empSD% above.\n")
