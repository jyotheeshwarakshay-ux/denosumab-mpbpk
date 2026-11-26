# ===========================================================================
# PAPER 2 - STEP 6a: simulation-estimation PILOT, design A (v2) vs design D.
# Read this block first.
#
#   WHY: PopED (Steps 5a/5b) predicts Kss RSE 78.9% (A) vs 33.6% (D) using
#   the FO approximation. FO can be wrong for a strongly nonlinear model.
#   Simulation-estimation (SSE) is the real test: simulate trials, refit
#   them, and see how close the estimates land to the truth.
#   This PILOT only checks that the pipeline works and measures the time per
#   fit, so the full Step 6b (many replicates) can be sized. Do NOT report
#   pilot numbers as results.
#
#   TRUTH: same model and values as R/20 (Choi 2025 Table 3 typical values
#   and IIV; Q fixed, no IIV; residual add 0.72 nM + prop 7%; 60 mg SC at
#   months 0/6/12; 80 subjects). Simulation uses R/20's deSolve model
#   (already checked against Python), so the truth is the validated model.
#   DESIGNS (h after each dose):
#     A v2: 12, 24, 72, 168, 336, 720, 1440, 2160
#     D   : 12, 24, 168, 720, 2160, 3000, 3600, 4200
#   Same random subjects (etas) in A and D for each replicate (paired).
#   BLQ: observed DV < LLOQ 0.1361 nM (20 ng/mL, Choi 2025 Sec 2.2) are
#   dropped (M1), as in Choi 2025 (BLQ set to MDV = 1). M3 = later sensitivity.
#   ESTIMATION: nlmixr2 SAEM (as Paper 1), same structural model, IIV on the
#   same 8 parameters as simulated, Q fixed, initial values = truth
#   (standard for design SSE: tests the information in the design, not the
#   search from poor starting values).
#
#   CHECK 0: the nlmixr2 model's typical Ctot must equal R/20's deSolve values
#   (38.3063 8.5498 40.1514 40.2174 9.0170 nM) to ~4 decimals. If not, STOP:
#   the fitted model is not the simulated model.
#
#   REQUIRES: nlmixr2 (+ rxode2), deSolve. Writes no files.
#   Run from the project root (denosumab-mpbpk/).
# ===========================================================================

N_REP <- 2   # pilot only

suppressMessages({ library(deSolve); library(nlmixr2); library(rxode2) })
cat("nlmixr2", as.character(packageVersion("nlmixr2")),
    "| nlmixr2est", as.character(packageVersion("nlmixr2est")),
    "| rxode2", as.character(packageVersion("rxode2")),
    "| deSolve", as.character(packageVersion("deSolve")), "\n\n")

# --- model block from R/20 (validated; prints its own typical-profile check)
src <- readLines("R/20_poped_baseline.R")
eval(parse(text = src[grep("^MW <-", src):(grep("^for \\(nm in", src) - 1)]))

LLOQ <- 0.1361
designs <- list(A = c(12, 24, 72, 168, 336, 720, 1440, 2160),
                D = c(12, 24, 168, 720, 2160, 3000, 3600, 4200))
om_true <- c(KA = cv2(56.57), VC = cv2(61.69), VP = cv2(15.55), CL = cv2(26.39),
             KINT = cv2(7.88), KSS = cv2(58.12), KSYN = cv2(22.46), R0 = cv2(158.3))

# --- fitted model (nlmixr2); literal initial values = truth ----------------
f_model <- function() {
  ini({
    lka   <- -4.853632   # log(0.0078)
    lvc   <- 0.457425    # log(1.58)
    lvp   <- 1.801710    # log(6.06)
    lcl   <- -5.115996   # log(0.006)
    lq    <- fixed(-1.609438)  # log(0.20), fixed as in Paper 1 v2
    lkint <- -3.816713   # log(0.022)
    lkss  <- 0.444686    # log(1.56)
    lksyn <- -4.605170   # log(0.01)
    lr0   <- 2.723267    # log(15.23)
    eta.ka   ~ 0.277644
    eta.vc   ~ 0.322493
    eta.vp   ~ 0.0238925
    eta.cl   ~ 0.0673251
    eta.kint ~ 0.00619024
    eta.kss  ~ 0.291022
    eta.ksyn ~ 0.049214
    eta.r0   ~ 1.25444
    add.sd  <- 0.72
    prop.sd <- 0.07
  })
  model({
    ka   <- exp(lka + eta.ka)
    vc   <- exp(lvc + eta.vc)
    vp   <- exp(lvp + eta.vp)
    cl   <- exp(lcl + eta.cl)
    q    <- exp(lq)
    kint <- exp(lkint + eta.kint)
    kss  <- exp(lkss + eta.kss)
    ksyn <- exp(lksyn + eta.ksyn)
    r0   <- exp(lr0 + eta.r0)
    kdeg <- ksyn / r0
    b    <- ctot - rtot - kss
    disc <- b * b + 4 * kss * ctot
    if (disc < 0) disc <- 0
    cfree <- 0.5 * (b + sqrt(disc))
    rc   <- ctot - cfree
    rtot(0) <- r0
    d/dt(depot) <- -ka * depot
    d/dt(ctot)  <- ka * depot / vc - (cl + q) * cfree / vc - kint * rc + q * cp / vp
    d/dt(cp)    <- q / vc * cfree - q / vp * cp
    d/dt(rtot)  <- ksyn - kdeg * rtot - (kint - kdeg) * rc
    ctot ~ add(add.sd) + prop(prop.sd)
  })
}

# initial values must equal R/20's truth
ui <- nlmixr2(f_model)
th <- ui$theta; omd <- diag(ui$omega)
stopifnot(abs(exp(th[c("lka","lvc","lvp","lcl","lq","lkint","lkss","lksyn","lr0")]) -
              bpop[c("KA","VC","VP","CL","Q","KINT","KSS","KSYN","R0")]) < 1e-5 * bpop,
          abs(omd - om_true) < 1e-5)

# --- CHECK 0: nlmixr2 model typical profile vs deSolve ---------------------
ev0 <- et(amt = dose_nmol, time = dose_times, cmt = "depot") |>
  et(c(300, 2160, 4683, 9066, 10926))
s0 <- rxSolve(ui, ev0, omega = NA, sigma = NA, returnType = "data.frame")
s0 <- s0[s0$time %in% c(300, 2160, 4683, 9066, 10926), ]
cat("CHECK 0 nlmixr2 typical Ctot:", sprintf("%.4f", s0$ctot), "\n")
cat("        R/20 deSolve expected: 38.3063 8.5498 40.1514 40.2174 9.0170\n")
if (any(abs(s0$ctot - c(38.3063, 8.5498, 40.1514, 40.2174, 9.0170)) > 1e-3)) stop("CHECK 0 FAILED")
cat("CHECK 0 passed\n\n")

# --- simulate one replicate (paired subjects across designs) --------------
sim_dataset <- function(offsets, rep, n_id = 80) {
  set.seed(1000 + rep)
  eta <- sapply(om_true, function(w) rnorm(n_id, 0, sqrt(w)))
  set.seed(500000 + rep)   # residual stream
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
  obs <- do.call(rbind, rows)
  dose <- expand.grid(ID = seq_len(n_id), TIME = dose_times)
  list(n_obs = nrow(obs), n_blq = sum(obs$BLQ),
       data = rbind(data.frame(ID = dose$ID, TIME = dose$TIME, EVID = 1, AMT = dose_nmol, CMT = "depot", DV = NA),
                    data.frame(ID = obs$ID[!obs$BLQ], TIME = obs$TIME[!obs$BLQ], EVID = 0, AMT = 0,
                               CMT = "ctot", DV = obs$DV[!obs$BLQ])) |>
         (\(d) d[order(d$ID, d$TIME, -d$EVID), ])())
}

pn <- c("lka","lvc","lvp","lcl","lkint","lkss","lksyn","lr0")
truth <- bpop[c("KA","VC","VP","CL","KINT","KSS","KSYN","R0")]

for (rep in seq_len(N_REP)) for (dn in names(designs)) {
  sd <- sim_dataset(designs[[dn]], rep)
  cat(sprintf("=== rep %d, design %s: %d obs, %d BLQ dropped (%.1f%%) ===\n",
              rep, dn, sd$n_obs, sd$n_blq, 100 * sd$n_blq / sd$n_obs))
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(
      nlmixr2(f_model, sd$data, est = "saem",
              control = saemControl(seed = rep, print = 0, nBurn = 300, nEm = 200)))),
    error = function(e) e))[["elapsed"]]
  if (inherits(fit, "error")) { cat("  FIT FAILED:", conditionMessage(fit), sprintf("(%.0f s)\n\n", tm)); next }
  est <- exp(fixef(fit)[pn]); se <- fit$parFixedDf[pn, "SE"]
  cat(sprintf("  time %.0f s | covariance: %s\n", tm, fit$covMethod))
  print(data.frame(true = signif(unname(truth), 4), est = signif(unname(est), 4),
                   err_pct = round(100 * (unname(est) / unname(truth) - 1), 1),
                   rse_pct = round(100 * se, 1), row.names = names(truth)))
  cat("  omega (est / true):", paste(sprintf("%s %.3f/%.3f", names(om_true),
                                             diag(fit$omega)[paste0("eta.", tolower(names(om_true)))], om_true),
                                     collapse = "; "), "\n\n")
}
cat("PILOT ONLY - use the timings to size Step 6b; do not report these estimates.\n")
