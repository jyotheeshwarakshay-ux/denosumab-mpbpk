# ===========================================================================
# PAPER 2 - STEP 3: compare the PBPK-TMDD model (2a, RANKL in plasma) with
# the published clinical 2-compartment QSS-TMDD model across doses.
# Read this block first.
#
#   REFERENCE MODEL: Choi et al. 2025 (Front Pharmacol 16:1631034, Table 3,
#   PMO typical values, CC-BY). Used because the Sutjandra 2011 full text
#   (Vp, Q, kint) is not available. Choi's CL, Vc, Vp, Q are apparent
#   (divided by F) and the dose enters without F, so its concentrations
#   are true concentrations and can be compared directly.
#   Sanity anchor: Choi 60 mg SC AUC0-16wk = 324 ug*day/mL vs Prolia label
#   316 ug*day/mL (label, mean of 46 subjects).
#
#   DOSES: 6, 60, 210 mg SC (the fixed-dose range in Sutjandra 2011).
#   Also a sensitivity variant: SAME TARGET AMOUNT as Choi's central
#   compartment (R0 = 15.23 x 1.58 / 3.126 = 7.70 nM in plasma).
#
#   EXPECTED RESULTS (independent Python, LSODA, rtol 1e-10, atol 1e-14):
#   Base PBPK (R0 15.23 nM) vs Choi:
#     6 mg   Cmax 0.278 vs 0.651 ug/mL; AUC0-300d 2.7 vs 11.3; ratio median 0.02
#     60 mg  Cmax 4.476 vs 5.631;       AUC0-300d 252.8 vs 344.1; ratio median 0.67
#     210 mg Cmax 18.645 vs 20.627;     AUC0-300d 1416.5 vs 1363.9; ratio median 1.21
#     60 mg AUC0-16wk: PBPK 249, Choi 324 (label 316)
#   Same-target-amount variant (R0 7.70 nM):
#     6 mg   Cmax 0.289, ratio median 0.09, AUC0-16wk 4 vs 11
#     60 mg  Cmax 5.065, ratio median 1.24, AUC0-16wk 316 vs 324
#     210 mg Cmax 19.237, ratio median 1.69, AUC0-16wk 1281 vs 1223
#   ("ratio" = PBPK / Choi plasma total concentration, from day 1 while
#   Choi > 0.01 ug/mL.)
#
#   WHAT THIS SHOWS (Python):
#   - Low dose (6 mg): PBPK far below Choi under both target assumptions.
#     The PBPK plasma volume (3.126 L) is about twice Choi's apparent
#     central volume (1.58 L), so early concentrations are lower and the
#     target removes a larger share of a small dose.
#   - High dose (210 mg, target saturated): PBPK above Choi late, because
#     the PBPK distribution volume (Vss 5.87 L) gives a longer linear tail.
#   - The mismatch has two separate structural causes (central volume and
#     distribution), which is exactly what Step 4 (reduction) must explain.
#   - Do NOT tune parameters to close these gaps.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/14_sc_dosing_rankl_plasma.R")   # gives mod2, p_on, MW, Kss, R0, kint, ksyn, kdeg

cat("\n\n=======================================================================\n")
cat("=== R/16_compare_clinical_model.R - Paper 2, Step 3 ===\n")
cat("=======================================================================\n\n")

# --- Choi 2025 two-compartment QSS-TMDD (apparent parameters) -------------
choi <- rxode2({
  bq <- Ctot - Rtot - Kss
  C <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot))
  RC <- Ctot - C
  d/dt(Depot) <- -cka * Depot
  d/dt(Ctot)  <- cka * Depot / cVc - (cCL + cQ) * C / cVc - kint * RC + cQ * Ap / cVp / cVc
  d/dt(Ap)    <- cQ * (C - Ap / cVp)
  d/dt(Rtot)  <- ksyn - kdeg * Rtot - (kint - kdeg) * RC
  Rtot(0) <- R0
})
p_choi <- c(cka = 0.0078, cVc = 1.58, cVp = 6.06, cCL = 0.006, cQ = 0.20,
            Kss = Kss, R0 = R0, kint = kint, ksyn = ksyn, kdeg = kdeg)

R0_amt <- 15.23 * 1.58 / V_central
p_amt  <- p_on; p_amt["R0"] <- R0_amt; p_amt["ksyn"] <- kdeg * R0_amt
cat(sprintf("Same-target-amount variant: R0 = %.2f nM   [Python: 7.70]\n\n", R0_amt))

tt <- seq(0, 24 * 300, by = 1)
run_pbpk <- function(dose_mg, p) {
  ev <- et(amt = dose_mg * 1e-3 / MW * 1e9, cmt = "Depot", time = 0) |> et(tt)
  r <- rxSolve(mod2, p, ev, atol = 1e-14, rtol = 1e-10)
  r$Central / V_central
}
run_choi <- function(dose_mg) {
  ev <- et(amt = dose_mg * 1e-3 / MW * 1e9, cmt = "Depot", time = 0) |> et(tt)
  r <- rxSolve(choi, p_choi, ev, atol = 1e-14, rtol = 1e-10)
  r$Ctot
}
auc <- function(c, upto = Inf) {
  k <- tt <= upto
  sum(diff(tt[k]) * (head(c[k], -1) + tail(c[k], -1)) / 2) / 24 * ug_per_mL
}
summ <- function(a, b) {
  w <- tt >= 24 & b * ug_per_mL > 0.01
  q <- a[w] / b[w]
  sprintf("ratio median %.2f (min %.2f, max %.2f)", median(q), min(q), max(q))
}

for (dm in c(6, 60, 210)) {
  a  <- run_pbpk(dm, p_on)
  a2 <- run_pbpk(dm, p_amt)
  b  <- run_choi(dm)
  cat(sprintf("=== %d mg SC ===\n", dm))
  cat(sprintf("Base PBPK : Cmax %.3f ug/mL day %.1f | AUC0-300d %.1f | AUC0-16wk %.0f | %s\n",
              max(a) * ug_per_mL, tt[which.max(a)] / 24, auc(a), auc(a, 24 * 112), summ(a, b)))
  cat(sprintf("Same amt  : Cmax %.3f ug/mL           | AUC0-300d %.1f | AUC0-16wk %.0f | %s\n",
              max(a2) * ug_per_mL, auc(a2), auc(a2, 24 * 112), summ(a2, b)))
  cat(sprintf("Choi 2025 : Cmax %.3f ug/mL day %.1f | AUC0-300d %.1f | AUC0-16wk %.0f\n\n",
              max(b) * ug_per_mL, tt[which.max(b)] / 24, auc(b), auc(b, 24 * 112)))
}
cat("Prolia label, 60 mg SC: Cmax 6.75 ug/mL, Tmax 10 days, AUC0-16wk 316 ug*day/mL\n")
cat("Compare every number with the EXPECTED RESULTS block at the top of this file.\n")
