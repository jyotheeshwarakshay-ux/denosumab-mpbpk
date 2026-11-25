# ===========================================================================
# PAPER 2 - STEP 4a: total plasma volume sensitivity.  Read this block first.
#
#   ONE CHANGE vs R/14 (Step 2a) and R/16 (Step 3): the central plasma
#   compartment is 3.126 - sum(V_V) = 1.4121 L instead of 3.126 L, so that
#   central plasma + organ vascular plasma = 3.126 L in total.
#
#   WHY: the Jones 2019 code (sources/PSP4-8-738-s002.txt, line 300) sets
#   V_Plasma = 3.126 L for the central compartment and adds the organ
#   vascular volumes V_V (lines 233-249, sum 1.714 L) as separate
#   compartments: total vascular plasma 4.84 L. The Shah-lab human
#   parameter table in Chang et al. 2022 (J Pharmacokinet Pharmacodyn,
#   doi 10.1007/s10928-022-09824-w, Table 1) totals 3.135 L of plasma
#   (tissue vascular 1.986 L + arterial/venous 1.149 L). This script tests
#   the "total plasma = 3.126 L" reading. It is a SENSITIVITY ANALYSIS, not
#   a claim that the Jones 2019 implementation is wrong.
#
#   Everything else as in Step 2a (RANKL QSS in plasma, R0 15.23 nM as a
#   CONCENTRATION; SC F 0.64, ka 0.00883 /h, CL 3.06 mL/h) and Step 3
#   (Choi 2025 reference model, doses 6/60/210 mg). Free concentrations are
#   reported because Sutjandra 2011 measured serum FREE denosumab.
#
#   EXPECTED RESULTS (independent Python, LSODA, rtol 1e-10, atol 1e-14):
#     central plasma 1.4121 L; Vss 4.151 L (Gibiansky 3.96 L/66 kg)
#     linear terminal half-life (IV, no target) 41.1 days
#     60 mg SC: mass balance ~8e-12; total Cmax 7.06 ug/mL at day 11.0;
#       free Cmax 6.99; t1/2 (day 30-120) 29.4 days; via RANKL by day 300
#       28.1%; AUC0-16wk total 374, free 363 ug*day/mL
#     PBPK/Choi median ratio (total, free): 6 mg 0.20, 0.11;
#       60 mg 1.08, 1.08; 210 mg 1.20, 1.21
#     Cmax PBPK vs Choi (total): 6 mg 0.578 vs 0.651; 60 mg 7.061 vs 5.631;
#       210 mg 26.451 vs 20.627
#   Prolia label 60 mg SC: Cmax 6.75 ug/mL, Tmax 10 d, t1/2 25.4 d,
#   AUC0-16wk 316 ug*day/mL.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/05_full_15_tissues.R")

cat("\n\n=======================================================================\n")
cat("=== R/17_total_plasma_volume.R - Paper 2, Step 4a ===\n")
cat("=======================================================================\n\n")

MW <- 147000; ug_per_mL <- MW * 1e-6
CL_lin <- 3.06 / 1000; F_sc <- 0.64; ka <- 0.00883            # Sutjandra 2011
Kss <- 1.56; R0 <- 15.23; kint <- 0.022; ksyn <- 0.01          # Choi 2025
kdeg <- ksyn / R0

V_c_new <- 3.126 - sum(V_V[all_organs])
p_vol <- organ_params; p_vol["V_central"] <- V_c_new
cat(sprintf("Central plasma = 3.126 - %.4f = %.4f L   [Python: 1.4121]\n", sum(V_V[all_organs]), V_c_new))

# --- Vss from the no-elimination equilibrium (file 05 model, new volume) ---
eq <- rxSolve(mod, p_vol, et(amt = 1, cmt = "Central", time = 0) |> et(c(0, 20000)))
eq_last <- eq[nrow(eq), ]
Vss <- sum(unlist(eq_last[state_cols])) / (eq_last$Central / V_c_new)
cat(sprintf("Vss = %.3f L   [Python: 4.151; Gibiansky 3.96 L/66 kg]\n", Vss))

# --- Linear terminal half-life (IV, no target) ---
i_ce <- grep("^d/dt\\(Central\\)", lines)
l_lin <- lines; l_lin[i_ce] <- paste(lines[i_ce], "- CL_lin * C_central")
mod_lin <- rxode2(paste(l_lin, collapse = "\n"))
rl <- rxSolve(mod_lin, c(p_vol, CL_lin = CL_lin),
              et(amt = 1, cmt = "Central", time = 0) |> et(seq(0, 20000, by = 10)),
              atol = 1e-14, rtol = 1e-10)
tl <- rl$time >= 10000
sl_lin <- coef(lm(log(rl$Central[tl] / V_c_new) ~ rl$time[tl]))[2]
cat(sprintf("Linear terminal half-life: %.1f days   [Python: 41.1]\n\n", log(2) / -unname(sl_lin) / 24))

# --- Step 2a model (RANKL QSS in plasma) with the new volume ---
i_cc <- grep("^C_central <-", lines)
nl <- lines
nl[i_cc] <- paste("Ctot <- Central / V_central", "bq <- Ctot - Rtot - Kss",
                  "C_central <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot))",
                  "RC <- Ctot - C_central", sep = "\n")
nl[i_ce] <- paste(lines[i_ce], "- CL_lin * C_central + ka * Depot - kint * RC * V_central")
nl <- c("Rtot(0) <- R0", "f(Depot) <- F_sc", nl, "d/dt(Depot) <- -ka * Depot",
        "d/dt(Rtot) <- ksyn - kdeg * Rtot - (kint - kdeg) * RC",
        "d/dt(Eliminated) <- CL_lin * C_central",
        "d/dt(TMDD_elim) <- kint * RC * V_central")
mod4 <- rxode2(paste(nl, collapse = "\n"))
p4 <- c(p_vol, CL_lin = CL_lin, F_sc = F_sc, ka = ka, Kss = Kss, R0 = R0,
        kint = kint, ksyn = ksyn, kdeg = kdeg)

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

tt <- seq(0, 24 * 300, by = 1)
qss_free <- function(Ct, R) { b <- Ct - R - Kss; 0.5 * (b + sqrt(b^2 + 4 * Kss * Ct)) }
run4 <- function(dose_mg) {
  d <- dose_mg * 1e-3 / MW * 1e9
  r <- rxSolve(mod4, p4, et(amt = d, cmt = "Depot", time = 0) |> et(tt), atol = 1e-14, rtol = 1e-10)
  Ct <- r$Central / V_c_new
  list(Ct = Ct, Cf = qss_free(Ct, r$Rtot), r = r, d = d)
}
runC <- function(dose_mg) {
  r <- rxSolve(choi, p_choi, et(amt = dose_mg * 1e-3 / MW * 1e9, cmt = "Depot", time = 0) |> et(tt),
               atol = 1e-14, rtol = 1e-10)
  list(Ct = r$Ctot, Cf = qss_free(r$Ctot, r$Rtot))
}
auc <- function(c, upto) { k <- tt <= upto; sum(diff(tt[k]) * (head(c[k], -1) + tail(c[k], -1)) / 2) / 24 * ug_per_mL }

o <- run4(60)
drug_cols <- c(state_cols, "Depot", "Eliminated", "TMDD_elim")
mb <- max(abs(rowSums(o$r[, drug_cols]) - F_sc * o$d)) / (F_sc * o$d)
k <- which.max(o$Ct); w <- tt >= 24 * 30 & tt <= 24 * 120
sl <- coef(lm(log(o$Ct[w]) ~ tt[w]))[2]
cat("=== 60 mg SC (RANKL in plasma, new volume) ===\n")
cat(sprintf("Mass balance rel. error %.1e   [Python: ~8e-12]\n", mb))
cat(sprintf("Total Cmax %.2f ug/mL at day %.1f | free Cmax %.2f   [Python: 7.06, day 11.0, 6.99]\n",
            o$Ct[k] * ug_per_mL, tt[k] / 24, max(o$Cf) * ug_per_mL))
cat(sprintf("Half-life (day 30-120) %.1f days | via RANKL by day 300 %.1f%%   [Python: 29.4, 28.1]\n",
            log(2) / -unname(sl) / 24, 100 * o$r$TMDD_elim[nrow(o$r)] / (F_sc * o$d)))
cat(sprintf("AUC0-16wk total %.0f, free %.0f ug*day/mL   [Python: 374, 363; label 316]\n\n",
            auc(o$Ct, 24 * 112), auc(o$Cf, 24 * 112)))

cat("=== PBPK (new volume) vs Choi 2025 ===\n")
for (dm in c(6, 60, 210)) {
  a <- run4(dm); b <- runC(dm)
  ww <- tt >= 24 & b$Ct * ug_per_mL > 0.01
  cat(sprintf("%3d mg: median ratio total %.2f, free %.2f | Cmax %.3f vs %.3f ug/mL\n", dm,
              median(a$Ct[ww] / b$Ct[ww]), median(a$Cf[ww] / b$Cf[ww]),
              max(a$Ct) * ug_per_mL, max(b$Ct) * ug_per_mL))
}
cat("[Python: 6 mg 0.20, 0.11, 0.578 vs 0.651 | 60 mg 1.08, 1.08, 7.061 vs 5.631 | 210 mg 1.20, 1.21, 26.451 vs 20.627]\n")
