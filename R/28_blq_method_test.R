# ===========================================================================
# PAPER 2 - STEP 6b-check 2: can ANY BLQ method in nlmixr2 FOCEi move Kss?
# Read this block first.
#
#   R/27 (rep 14, one trial): with M3 (CENS column) Kss and R0 stayed at
#   whatever value they started from (0.78 -> 0.798; 3.12 -> 3.214; truth ->
#   1.56000) with tiny SEs, while other parameters moved. Without the CENS
#   column (A_M1) Kss moved (1.56 -> 2.036, SE 0.070). All fits used bobyqa
#   (the nlmixr2 default), so R/27's "bobyqa" fit tested nothing new.
#
#   TEST on the SAME trial (rep 14, design D), ALL fits start AWAY from the
#   truth (structural x1.3, Kss x2 = 3.12; as R/27 D_M3_high):
#     1 D_M3_nlminb    M3, gradient-based outer optimiser nlminb
#     2 D_M5           Beal M5: BLQ -> LLOQ/2, no CENS column, default optimiser
#     3 D_M3_lbfgsb3c  M3, gradient-based outer optimiser lbfgsb3c
#     4 D_M1           BLQ dropped, default optimiser
#                      (R/25 truth start gave Kss 4.07653 -> same value = M1 moves
#                      to a real optimum, independent of start)
#   READING:
#     1 or 3 leave 3.12 (SE ~0.1-0.4) -> use that optimiser for M3 in Step 6b.
#     1 and 3 stay near 3.12 -> M3 is not usable in nlmixr2 7.0.x for this
#        model; Step 6b uses M1 vs M5 (+ ORACLE-type reference) and the paper
#        states the limitation.
#     4 returns to ~4.08 -> the M1 bias is a property of the data, not the start.
#
#   2 PSOCK workers, jobs in order. ~3-4 h. Writes no files.
#   Run from the project root.
# ===========================================================================

suppressMessages(library(parallel))

setup_code <- function() {
  s <- readLines("R/23_blq_diagnostic.R")
  invisible(capture.output(eval(parse(text = s[grep("^suppressMessages", s)[1]:(grep("^fit_one <-", s) - 1)]),
                                envir = globalenv())))
}

make_data <- function(dn, rep, meth) {
  obs <- sim_obs(designs[[dn]], rep)
  dose <- expand.grid(ID = 1:80, TIME = dose_times)
  dr <- data.frame(ID = dose$ID, TIME = dose$TIME, EVID = 1, AMT = dose_nmol, CMT = "depot", DV = NA)
  if (meth == "M5") {   # Beal M5: every BLQ value replaced by LLOQ/2, no CENS column
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
  d <- rbind(dr, orows); d[order(d$ID, d$TIME, -d$EVID), ]
}

run_job <- function(job) {
  setup_code()
  b <- bpop
  m_high <- ui |> ini(lka = log(b[["KA"]] * 1.3), lvc = log(b[["VC"]] * 1.3), lvp = log(b[["VP"]] * 1.3),
                      lcl = log(b[["CL"]] * 1.3), lkint = log(b[["KINT"]] * 1.3),
                      lkss = log(b[["KSS"]] * 2), lksyn = log(b[["KSYN"]] * 1.3), lr0 = log(b[["R0"]] * 1.3))
  spec <- switch(job,
    D_M3_nlminb   = list(d = make_data("D", 14, "M3"), ctl = foceiControl(print = 0, outerOpt = "nlminb")),
    D_M5          = list(d = make_data("D", 14, "M5"), ctl = foceiControl(print = 0)),
    D_M3_lbfgsb3c = list(d = make_data("D", 14, "M3"), ctl = foceiControl(print = 0, outerOpt = "lbfgsb3c")),
    D_M1          = list(d = make_data("D", 14, "M1"), ctl = foceiControl(print = 0)))
  st <- exp(m_high$theta[pn])
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(nlmixr2(m_high, spec$d, est = "focei", control = spec$ctl))),
    error = function(e) e))[["elapsed"]]
  if (inherits(fit, "error")) return(list(job = job, err = conditionMessage(fit), tm = tm, start = st))
  list(job = job, tm = tm, start = st, est = exp(fixef(fit)[pn]),
       se_kss = fit$parFixedDf["lkss", "SE"], objf = fit$objf, cov = fit$covMethod,
       n_obs = sum(spec$d$EVID == 0), cens = if (is.null(spec$d$CENS)) "none" else sum(spec$d$CENS),
       iter = tryCatch(nrow(fit$parHistData), error = function(e) NA),
       msg = tryCatch(paste(fit$message, collapse = " "), error = function(e) NA))
}

jobs <- c("D_M3_nlminb", "D_M5", "D_M3_lbfgsb3c", "D_M1")
cl <- makeCluster(2)
clusterCall(cl, setwd, getwd())
clusterExport(cl, c("setup_code", "make_data"))
res <- parLapplyLB(cl, jobs, run_job, chunk.size = 1)
stopCluster(cl)

for (r in res) {
  cat(sprintf("=== %s ===\n", r$job))
  if (!is.null(r$err)) { cat("  FIT FAILED:", r$err, sprintf("(%.0f s)\n\n", r$tm)); next }
  cat(sprintf("  time %.0f s | obs %d | censored %s | cov %s | objf %.2f | outer iter %s\n",
              r$tm, r$n_obs, r$cens, r$cov, r$objf, r$iter))
  cat("  optimiser message:", r$msg, "\n")
  print(data.frame(start = signif(unname(r$start), 6), est = signif(unname(r$est), 6),
                   row.names = toupper(sub("^l", "", names(r$est)))))
  cat(sprintf("  SE(log Kss) %.4g\n\n", r$se_kss))
}
cat("[rep 14: R/27 D_M3 default optimiser from these starts -> Kss 3.21375 objf 4294.74 | R/25 D_M1 truth start -> Kss 4.07653]\n")
cat("DIAGNOSTIC ONLY - one trial; do not report as a result.\n")
