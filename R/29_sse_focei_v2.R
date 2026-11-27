# ===========================================================================
# PAPER 2 - STEP 6b (v2): SSE with methods that WORK in nlmixr2 FOCEi, and
# starting values AWAY from the truth. Supersedes R/25 (truth starts; M3).
# Read this block first.
#
#   WHY v2 (diagnostics R/22-R/28, one trial each, not results):
#   - M3 (CENS column) in nlmixr2 7.0.x FOCEi: Kss and R0 stay at their
#     starting value with every outer optimiser tried (bobyqa default,
#     nlminb "false convergence", lbfgsb3c) -> M3 is not used here; M3 will be
#     done separately in Monolix (censored likelihood under SAEM).
#   - R/25 started every fit at the truth -> M3 fits returned the truth.
#   - M1 (BLQ dropped) and M5 (BLQ -> LLOQ/2) fits move and converge normally,
#     but the M1 likelihood is flat in Kss (rep 14: Kss 4.08 vs 4.74, objf
#     difference 0.83) -> single-start estimates carry optimiser noise.
#
#   JOBS per replicate (paired subjects, as R/25):
#     A_M1 : design A, BLQ dropped (A has ~0.2% BLQ, so ~= all data)
#     D_M1 : design D, BLQ dropped (Choi 2025 practice)
#     D_M5 : design D, BLQ -> LLOQ/2 (Beal 2001 M5)
#   STARTS: drawn per replicate (seed 900000 + rep), shared by the 3 jobs:
#     each structural parameter x exp(U(-0.4, 0.4)); Kss x exp(U(-0.7, 0.7))
#     (0.5x-2x). Q stays fixed. IIV/residual starts as R/22. Starts are saved.
#   ESTIMATOR: FOCEi (nlmixr2 default outer optimiser).
#   OUTPUT: results/step6b_v2/<job>_rep<NNN>.rds + progress.log; reruns skip
#   finished fits; jobs in replicate order (chunk.size = 1).
#   DEFAULT: 30 replicates, 2 workers (heat). Summary: R/30.
#   USAGE (project root):  Rscript R/29_sse_focei_v2.R [N_REP] [N_WORKERS]
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
N_REP     <- if (length(args) >= 1) as.integer(args[1]) else 30
N_WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 2

suppressMessages(library(parallel))
out_dir <- "results/step6b_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

jobs <- do.call(rbind, lapply(seq_len(N_REP), function(r)
  data.frame(rep = r, job = c("A_M1", "D_M1", "D_M5"))))
jobs$file <- file.path(out_dir, sprintf("%s_rep%03d.rds", jobs$job, jobs$rep))
todo <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Step 6b v2: %d replicates x 3 jobs = %d fits; %d already done; %d to run on %d workers\n",
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
  if (meth == "M5") {
    orows <- data.frame(ID = obs$ID, TIME = obs$TIME, EVID = 0, AMT = 0, CMT = "ctot",
                        DV = ifelse(obs$BLQ, LLOQ / 2, obs$DV))
  } else if (meth == "M1") {
    o <- obs[!obs$BLQ, ]
    orows <- data.frame(ID = o$ID, TIME = o$TIME, EVID = 0, AMT = 0, CMT = "ctot", DV = o$DV)
  } else {
    dr$CENS <- 0
    orows <- data.frame(ID = obs$ID, TIME = obs$TIME, EVID = 0, AMT = 0, CMT = "ctot",
                        DV = ifelse(obs$BLQ, LLOQ, obs$DV), CENS = as.integer(obs$BLQ))
  }
  d <- rbind(dr, orows); d <- d[order(d$ID, d$TIME, -d$EVID), ]
  set.seed(900000 + rep)
  u <- runif(7, -0.4, 0.4); uk <- runif(1, -0.7, 0.7)
  st <- c(KA = bpop[["KA"]] * exp(u[1]), VC = bpop[["VC"]] * exp(u[2]), VP = bpop[["VP"]] * exp(u[3]),
          CL = bpop[["CL"]] * exp(u[4]), KINT = bpop[["KINT"]] * exp(u[5]), KSS = bpop[["KSS"]] * exp(uk),
          KSYN = bpop[["KSYN"]] * exp(u[6]), R0 = bpop[["R0"]] * exp(u[7]))
  m0 <- ui |> ini(lka = log(st[["KA"]]), lvc = log(st[["VC"]]), lvp = log(st[["VP"]]), lcl = log(st[["CL"]]),
                  lkint = log(st[["KINT"]]), lkss = log(st[["KSS"]]), lksyn = log(st[["KSYN"]]), lr0 = log(st[["R0"]]))
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(nlmixr2(m0, d, est = "focei", control = foceiControl(print = 0)))),
    error = function(e) e))[["elapsed"]]
  r <- list(rep = rep, job = job, design = dn, method = meth, time_s = tm,
            n_obs = nrow(obs), n_blq = sum(obs$BLQ), truth = truth, start = st, om_true = om_true,
            versions = c(nlmixr2 = as.character(packageVersion("nlmixr2")),
                         nlmixr2est = as.character(packageVersion("nlmixr2est")),
                         rxode2 = as.character(packageVersion("rxode2"))))
  if (inherits(fit, "error")) {
    r$error <- conditionMessage(fit)
  } else {
    r$est <- exp(fixef(fit)[pn]); r$se_log <- fit$parFixedDf[pn, "SE"]
    r$om <- diag(fit$omega)[paste0("eta.", tolower(names(om_true)))]
    r$cov <- fit$covMethod; r$objf <- fit$objf
    r$msg <- tryCatch(paste(fit$message, collapse = " "), error = function(e) NA)
    r$sigma <- fit$theta[c("add.sd", "prop.sd")]
  }
  saveRDS(r, todo$file[i])
  msg <- sprintf("%s rep %d: %s (%.0f s)%s", job, rep, if (is.null(r$error)) "done" else "FAILED",
          tm, if (is.null(r$error)) sprintf(" Kss %.3f (start %.3f) cov %s", r$est[["lkss"]], st[["KSS"]], r$cov) else paste(":", r$error))
  cat(format(Sys.time(), "%d %b %H:%M"), msg, "\n", file = file.path(out_dir, "progress.log"), append = TRUE)
  msg
}

cl <- makeCluster(N_WORKERS)
clusterCall(cl, setwd, getwd())
clusterExport(cl, c("setup_code", "todo", "out_dir"))
invisible(clusterEvalQ(cl, setup_code()))   # each worker: R/23 setup incl. CHECK 0 and data checks
msgs <- parLapplyLB(cl, seq_len(nrow(todo)), run_job, chunk.size = 1)
stopCluster(cl)

cat(unlist(msgs), sep = "\n")
cat("\nFinished:", format(Sys.time()), "- summarise with: Rscript R/30_sse_summary_v2.R\n")
