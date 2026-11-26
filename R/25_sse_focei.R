# ===========================================================================
# PAPER 2 - STEP 6b: simulation-estimation (SSE), design A (v2) vs design D.
# Read this block first.
#
#   DECISIONS FROM THE DIAGNOSTICS (R/22-R/24, one dataset, not results):
#   - nlmixr2 SAEM is not reproducible here even with a fixed seed (same data
#     gave Kss 3.06-4.24), and SAEM with censoring (M3) collapsed.
#   - FOCEi with M3 gave sensible estimates and a successful covariance step.
#   -> ESTIMATOR for every fit: FOCEi (deterministic), start = truth.
#
#   JOBS (each replicate = a new simulated trial of 80 subjects; the same
#   subjects are used for A and D within a replicate = paired):
#     A_M3 : design A, BLQ kept as censored (M3)       - primary
#     D_M3 : design D, BLQ kept as censored (M3)       - primary
#     D_M1 : design D, BLQ dropped (M1, as Choi 2025)  - size of the BLQ-drop bias
#   (Design A has ~0.2% BLQ, so A_M1 would equal A_M3 and is not run.)
#   Jobs are ordered replicate by replicate, so stopping early still leaves a
#   balanced set. Each finished fit is saved at once to
#   results/step6b/<job>_rep<NNN>.rds; re-running SKIPS finished fits
#   (safe to stop and restart at any time).
#
#   SIMULATION, MODEL, CHECKS: R/23's setup lines (R/20 model block, R/22
#   nlmixr2 model + CHECK 0, R/23 simulator), re-run inside each worker.
#   PARALLEL: PSOCK cluster, N_WORKERS workers (4 = the M1's performance
#   cores; R/24 peak RAM ~0.9 GB for 3 fits). ~70 min per FOCEi fit.
#   N_REP 50 -> 150 fits -> ~45 h with 4 workers. Summaries: R/26.
#
#   USAGE (project root):  Rscript R/25_sse_focei.R [N_REP] [N_WORKERS]
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
N_REP     <- if (length(args) >= 1) as.integer(args[1]) else 50
N_WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 4

suppressMessages(library(parallel))
out_dir <- "results/step6b"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

jobs <- do.call(rbind, lapply(seq_len(N_REP), function(r)
  data.frame(rep = r, job = c("A_M3", "D_M3", "D_M1"))))
jobs$file <- file.path(out_dir, sprintf("%s_rep%03d.rds", jobs$job, jobs$rep))
todo <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Step 6b: %d replicates x 3 jobs = %d fits; %d already done; %d to run on %d workers\n",
            N_REP, nrow(jobs), nrow(jobs) - nrow(todo), nrow(todo), N_WORKERS))
cat("Started:", format(Sys.time()), "\n\n")
if (nrow(todo) == 0) quit(save = "no")

setup_code <- function() {
  s <- readLines("R/23_blq_diagnostic.R")
  invisible(capture.output(eval(parse(text = s[grep("^suppressMessages", s)[1]:(grep("^fit_one <-", s) - 1)]),
                                envir = globalenv())))
}

run_job <- function(i) {
  rep <- todo$rep[i]; job <- todo$job[i]; dn <- substr(job, 1, 1); meth <- substr(job, 3, 4)
  obs <- sim_obs(designs[[dn]], rep)
  dose <- expand.grid(ID = 1:80, TIME = dose_times)
  dr <- data.frame(ID = dose$ID, TIME = dose$TIME, EVID = 1, AMT = dose_nmol, CMT = "depot", DV = NA)
  if (meth == "M1") {
    o <- obs[!obs$BLQ, ]
    orows <- data.frame(ID = o$ID, TIME = o$TIME, EVID = 0, AMT = 0, CMT = "ctot", DV = o$DV)
  } else {
    dr$CENS <- 0
    orows <- data.frame(ID = obs$ID, TIME = obs$TIME, EVID = 0, AMT = 0, CMT = "ctot",
                        DV = ifelse(obs$BLQ, LLOQ, obs$DV), CENS = as.integer(obs$BLQ))
  }
  d <- rbind(dr, orows); d <- d[order(d$ID, d$TIME, -d$EVID), ]
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(nlmixr2(f_model, d, est = "focei", control = foceiControl(print = 0)))),
    error = function(e) e))[["elapsed"]]
  r <- list(rep = rep, job = job, design = dn, method = meth, time_s = tm,
            n_obs = nrow(obs), n_blq = sum(obs$BLQ), truth = truth, om_true = om_true,
            versions = c(nlmixr2 = as.character(packageVersion("nlmixr2")),
                         nlmixr2est = as.character(packageVersion("nlmixr2est")),
                         rxode2 = as.character(packageVersion("rxode2"))))
  if (inherits(fit, "error")) {
    r$error <- conditionMessage(fit)
  } else {
    r$est <- exp(fixef(fit)[pn]); r$se_log <- fit$parFixedDf[pn, "SE"]
    r$om <- diag(fit$omega)[paste0("eta.", tolower(names(om_true)))]
    r$cov <- fit$covMethod; r$objf <- fit$objf
    r$sigma <- fit$theta[c("add.sd", "prop.sd")]
  }
  saveRDS(r, todo$file[i])
  msg <- sprintf("%s rep %d: %s (%.0f s)%s", job, rep, if (is.null(r$error)) "done" else "FAILED",
          tm, if (is.null(r$error)) sprintf(" Kss %.3f cov %s", r$est[["lkss"]], r$cov) else paste(":", r$error))
  cat(format(Sys.time(), "%d %b %H:%M"), msg, "\n", file = file.path(out_dir, "progress.log"), append = TRUE)
  msg
}

cl <- makeCluster(N_WORKERS)
clusterCall(cl, setwd, getwd())
clusterExport(cl, c("setup_code", "todo", "out_dir"))
invisible(clusterEvalQ(cl, setup_code()))   # each worker: R/23 setup incl. CHECK 0 and data checks
msgs <- parLapplyLB(cl, seq_len(nrow(todo)), run_job)
stopCluster(cl)

cat(unlist(msgs), sep = "\n")
cat("\nFinished:", format(Sys.time()), "- summarise with: Rscript R/26_sse_summary.R\n")
