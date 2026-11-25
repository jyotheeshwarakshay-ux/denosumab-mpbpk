# ===========================================================================
# PAPER 2 - STEP 5a: how well can Paper 1's sampling designs estimate Kss?
# Read this block first.
#
#   WHAT: PopED (Nyberg et al. 2012, Comput Methods Programs Biomed
#   108:789-805) computes the Fisher information matrix (FIM) of a design
#   and from it the expected relative standard error (RSE%) of every
#   parameter, before any data are simulated or fitted.
#
#   MODEL: the clinical 2-compartment QSS-TMDD model with Choi et al. 2025
#   Table 3 typical values and between-subject variability (as in Paper 1).
#   Q fixed (no IIV), exactly as in the Paper 1 v2 fit. Residual error:
#   additive 0.72 nM + proportional 7% (Choi Table 3).
#   DOSING: 60 mg SC at months 0, 6, 12 (30.44-day months), 80 subjects.
#   DESIGNS (samples after EACH dose, hours):
#     v2 (Paper 1 final): 12, 24, 72, 168, 336, 720, 1440, 2160  (24 samples)
#     v1 (Paper 1 first): 12, 24, 168, 720, 2160, 3600, 4380     (21 samples)
#   FIM: PopED default (FO approximation, reduced FIM). FO can be inaccurate
#   for strongly nonlinear models; Step 5b will check FOCE, and Step 6 will
#   check by simulation and re-estimation.
#
#   REQUIRES (once): install.packages(c("PopED", "deSolve"))
#
#   CHECKS + EXPECTED RESULTS (run in R 4.3.3, PopED 0.7.0.9002 from GitHub,
#   deSolve 1.40; the model was also checked against an independent Python
#   LSODA solver, agreeing to 4 decimals):
#     typical Ctot (nM) at 300, 2160, 4683, 9066, 10926 h:
#       38.3063  8.5498  40.1514  40.2174  9.0170
#     v2 RSE% (fixed effects): KA 6.2, VC 9.2, VP 2.3, CL 3.3, KINT 8.4,
#       KSS 79.8, KSYN 6.8, R0 14.1
#     v2 RSE% (IIV variances): KA 16.2, VC 17.9, VP 18.7, CL 16.4,
#       KINT 841.9, KSS 713.0, KSYN 54.4, R0 17.4; sigma prop 6.3, add 18.9
#     v1 RSE% (fixed effects): KA 6.1, VC 7.0, VP 1.9, CL 3.0, KINT 13.8,
#       KSS 35.4, KSYN 4.3, R0 14.5
#     v1 RSE% (IIV variances): KA 16.5, VC 17.9, VP 18.2, CL 16.5,
#       KINT 2555.3, KSS 352.6, KSYN 43.3, R0 17.3; sigma prop 7.1, add 6.8
#   Small differences with a different PopED version are possible; report
#   them, do not tune.
# CONFIRMED on the author's machine with CRAN PopED 0.7.0 + deSolve 1.42 (26 Nov 2026):
# typical Ctot identical to 4 decimals; v2 KSS 78.9, v1 KSS 35.1;
# d_KSS v2 722.7, v1 340.3; d_KINT v2 845.8, v1 2640.5.
# REFERENCE FOR THE MANUSCRIPT = CRAN PopED 0.7.0 values.
#
# Run from the project root (denosumab-mpbpk/).
# ===========================================================================

suppressMessages({ library(PopED); library(deSolve) })
cat("PopED", as.character(packageVersion("PopED")), "| deSolve", as.character(packageVersion("deSolve")), "\n\n")

MW <- 147000; dose_nmol <- 60e-3 / MW * 1e9; month_h <- 30.44 * 24
dose_times <- c(0, 6, 12) * month_h

ode_qss <- function(t, y, p) {
  with(as.list(c(y, p)), {
    b <- Ctot - Rtot - Kss; C <- 0.5 * (b + sqrt(b * b + 4 * Kss * Ctot)); RC <- Ctot - C
    list(c(-KA * depot,
           KA * depot / VC - (CL + Q) * C / VC - KINT * RC + Q * Cp / VP,
           Q / VC * C - Q / VP * Cp,
           KSYN - KDEG * Rtot - (KINT - KDEG) * RC))
  })
}
ff <- function(model_switch, xt, parameters, poped.db) {
  with(as.list(parameters), {
    p <- c(KA = KA, VC = VC, VP = VP, CL = CL, Q = Q, KINT = KINT, Kss = KSS, KSYN = KSYN, KDEG = KSYN / R0)
    times <- sort(unique(c(0, dose_times, xt)))
    ev <- data.frame(var = "depot", time = dose_times, value = dose_nmol, method = "add")
    out <- suppressWarnings(ode(c(depot = 0, Ctot = 0, Cp = 0, Rtot = R0), times, ode_qss, p,
                                events = list(data = ev), method = "lsoda", rtol = 1e-8, atol = 1e-10))
    list(f = matrix(out[match(xt, out[, "time"]), "Ctot"], ncol = 1), poped.db = poped.db)
  })
}
sfg <- function(x, a, bpop, b, bocc) {
  c(KA = bpop[1] * exp(b[1]), VC = bpop[2] * exp(b[2]), VP = bpop[3] * exp(b[3]),
    CL = bpop[4] * exp(b[4]), Q = bpop[5], KINT = bpop[6] * exp(b[5]),
    KSS = bpop[7] * exp(b[6]), KSYN = bpop[8] * exp(b[7]), R0 = bpop[9] * exp(b[8]))
}
cv2 <- function(cv) log(1 + (cv / 100)^2)
bpop <- c(KA = 0.0078, VC = 1.58, VP = 6.06, CL = 0.006, Q = 0.20, KINT = 0.022, KSS = 1.56, KSYN = 0.01, R0 = 15.23)

# --- Check: typical profile ---------------------------------------------
chk <- ff(1, c(300, 2160, 4683, 9066, 10926), bpop, NULL)$f[, 1]
cat("Typical Ctot (nM) at 300/2160/4683/9066/10926 h:", sprintf("%.4f", chk), "\n")
cat("[expected: 38.3063 8.5498 40.1514 40.2174 9.0170]\n\n")

make_db <- function(offsets) {
  xt <- sort(as.vector(outer(dose_times, offsets, "+")))
  create.poped.database(
    ff_fun = ff, fg_fun = sfg, fError_fun = feps.add.prop,
    bpop = bpop, notfixed_bpop = c(1, 1, 1, 1, 0, 1, 1, 1, 1),
    d = c(KA = cv2(56.57), VC = cv2(61.69), VP = cv2(15.55), CL = cv2(26.39), KINT = cv2(7.88),
          KSS = cv2(58.12), KSYN = cv2(22.46), R0 = cv2(158.3)),
    sigma = c(prop = 0.07^2, add = 0.72^2),
    m = 1, groupsize = 80, xt = xt, minxt = 0, maxxt = max(xt))
}

for (nm in c("v2", "v1")) {
  offs <- if (nm == "v2") c(12, 24, 72, 168, 336, 720, 1440, 2160) else c(12, 24, 168, 720, 2160, 3600, 4380)
  ev <- evaluate_design(make_db(offs))
  cat(sprintf("=== Design %s (%d samples per subject) - expected RSE%% ===\n", nm, 3 * length(offs)))
  print(round(ev$rse, 1))
  cat("\n")
}
cat("[expected KSS: v2 79.8, v1 35.4 | d_KSS: v2 713.0, v1 352.6 - see header for all values]\n")
