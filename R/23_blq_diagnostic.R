# ===========================================================================
# PAPER 2 - STEP 6a-2: why did design D give Kss ~4.25 (truth 1.56)?
# Read this block first.
#
#   PILOT (R/22): design D gave Kss 4.237 (rep 1) and 4.262 (rep 2) - nearly
#   the same wrong value from two different datasets = systematic, not noise.
#   HYPOTHESIS: dropping BLQ samples (M1) removes the LOWEST late
#   concentrations; the remaining late data look too high (slow decline at
#   low concentration); the model explains slow decline by weaker binding =
#   a HIGHER Kss. Design A has almost no BLQ, so it is barely affected.
#
#   TEST: the SAME dataset (rep 1, design D; 113 of 1920 obs BLQ) fitted 3
#   ways, in parallel (3 cores):
#     M1     - BLQ dropped (repeat of the pilot; must reproduce Kss 4.237)
#     ORACLE - nothing dropped (every simulated value kept; impossible in a
#              real trial, but isolates the effect of dropping)
#     M3     - BLQ kept as censored (CENS = 1, DV = LLOQ; Beal 2001 M3)
#   READING: ORACLE and M3 near 1.56 while M1 ~4.24 -> M1 causes the bias.
#   All three ~4.24 -> the bias is NOT from BLQ (look at the design/model).
#
#   Same simulation, seeds, model, SAEM settings and starting values as R/22.
#   CHECK: BLQ count for rep 1 design D must be 113 (as in the pilot).
#   Expected run time ~1 h (fits run in parallel). Writes no files.
#   Parallel fits use a PSOCK cluster (forked mclapply workers fail to build the model on macOS).
#   Run from the project root (denosumab-mpbpk/).
# ===========================================================================

suppressMessages({ library(deSolve); library(nlmixr2); library(rxode2); library(parallel) })
cat("nlmixr2", as.character(packageVersion("nlmixr2")),
    "| nlmixr2est", as.character(packageVersion("nlmixr2est")),
    "| rxode2", as.character(packageVersion("rxode2")),
    "| cores", detectCores(), "\n\n")

src <- readLines("R/20_poped_baseline.R")
eval(parse(text = src[grep("^MW <-", src):(grep("^for \\(nm in", src) - 1)]))
s22 <- readLines("R/22_sse_pilot.R")   # model, CHECK 0, simulator (validated in the pilot)
eval(parse(text = s22[grep("^LLOQ <-", s22):(grep("^for \\(rep in", s22) - 1)]))

# full simulated observations (same seeds and code path as sim_dataset in R/22)
sim_obs <- function(offsets, rep, n_id = 80) {
  set.seed(1000 + rep)
  eta <- sapply(om_true, function(w) rnorm(n_id, 0, sqrt(w)))
  set.seed(500000 + rep)
  xt <- sort(as.vector(outer(dose_times, offsets, "+")))
  rows <- vector("list", n_id)
  for (i in seq_len(n_id)) {
    p_i <- c(KA = bpop[["KA"]] * exp(eta[i, "KA"]), VC = bpop[["VC"]] * exp(eta[i, "VC"]),
             VP = bpop[["VP"]] * exp(eta[i, "VP"]), CL = bpop[["CL"]] * exp(eta[i, "CL"]),
             Q = bpop[["Q"]], KINT = bpop[["KINT"]] * exp(eta[i, "KINT"]),
             KSS = bpop[["KSS"]] * exp(eta[i, "KSS"]), KSYN = bpop[["KSYN"]] * exp(eta[i, "KSYN"]),
             R0 = bpop[["R0"]] * exp(eta[i, "R0"]))
    p_i <- setNames(as.numeric(p_i), c("KA", "VC", "VP", "CL", "Q", "KINT", "KSS", "KSYN", "R0"))
    ipred <- ff(1, xt, p_i, NULL)$f[, 1]
    dv <- ipred * (1 + 0.07 * rnorm(length(xt))) + 0.72 * rnorm(length(xt))
    rows[[i]] <- data.frame(ID = i, TIME = xt, DV = dv, BLQ = dv < LLOQ)
  }
  do.call(rbind, rows)
}

obs <- sim_obs(designs$D, 1)
cat("CHECK: rep 1 design D BLQ =", sum(obs$BLQ), "of", nrow(obs), "[expected 113 of 1920]\n")
stopifnot(sum(obs$BLQ) == 113)
ref <- sim_dataset(designs$D, 1)$data
stopifnot(isTRUE(all.equal(ref$DV[ref$EVID == 0], obs$DV[!obs$BLQ])))
cat("CHECK: M1 data identical to the pilot's - passed\n\n")

dose <- expand.grid(ID = 1:80, TIME = dose_times)
dose_rows <- data.frame(ID = dose$ID, TIME = dose$TIME, EVID = 1, AMT = dose_nmol, CMT = "depot", DV = NA, CENS = 0)
mk <- function(o) { d <- rbind(dose_rows, o); d[order(d$ID, d$TIME, -d$EVID), ] }
obs_rows <- function(o, dv, cens) data.frame(ID = o$ID, TIME = o$TIME, EVID = 0, AMT = 0, CMT = "ctot", DV = dv, CENS = cens)
datasets <- list(
  M1     = mk(obs_rows(obs[!obs$BLQ, ], obs$DV[!obs$BLQ], 0)),
  ORACLE = mk(obs_rows(obs, obs$DV, 0)),
  M3     = mk(obs_rows(obs, ifelse(obs$BLQ, LLOQ, obs$DV), as.integer(obs$BLQ))))
cat("Observation rows: M1", sum(datasets$M1$EVID == 0), "| ORACLE", sum(datasets$ORACLE$EVID == 0),
    "| M3", sum(datasets$M3$EVID == 0), "(censored", sum(datasets$M3$CENS), ")\n\n")

fit_one <- function(nm) {
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(
      nlmixr2(f_model, datasets[[nm]], est = "saem",
              control = saemControl(seed = 1, print = 0, nBurn = 300, nEm = 200)))),
    error = function(e) e))[["elapsed"]]
  if (inherits(fit, "error")) return(list(nm = nm, err = conditionMessage(fit), tm = tm))
  list(nm = nm, tm = tm, cov = fit$covMethod, est = exp(fixef(fit)[pn]),
       se = fit$parFixedDf[pn, "SE"], om = diag(fit$omega)[paste0("eta.", tolower(names(om_true)))])
}
# PSOCK workers (fresh R sessions): forked mclapply workers fail to build the
# nlmixr2 model on macOS. Each worker re-runs the setup above (same code, seeds, data).
cl <- makeCluster(3)
clusterCall(cl, setwd, getwd())
clusterEvalQ(cl, { s <- readLines("R/23_blq_diagnostic.R")
  invisible(capture.output(eval(parse(text = s[grep("^suppressMessages", s)[1]:(grep("^cl <- makeCluster", s) - 1)]))))
  NULL })
res <- parLapply(cl, names(datasets), function(nm) fit_one(nm))
stopCluster(cl)

for (r in res) {
  cat(sprintf("=== %s ===\n", r$nm))
  if (!is.null(r$err)) { cat("  FIT FAILED:", r$err, sprintf("(%.0f s)\n\n", r$tm)); next }
  cat(sprintf("  time %.0f s | covariance: %s\n", r$tm, r$cov))
  print(data.frame(true = signif(unname(truth), 4), est = signif(unname(r$est), 4),
                   err_pct = round(100 * (unname(r$est) / unname(truth) - 1), 1),
                   rse_pct = round(100 * r$se, 1), row.names = names(truth)))
  cat("  omega KSS (est / true):", sprintf("%.3f / %.3f", r$om[["eta.kss"]], om_true[["KSS"]]), "\n\n")
}
cat("[pilot M1 Kss for this dataset: 4.237]\n")
cat("DIAGNOSTIC ONLY - one dataset; do not report as a result.\n")
