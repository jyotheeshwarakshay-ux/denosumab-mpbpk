# ===========================================================================
# PAPER 2 - STEP 2a: subcutaneous dosing + RANKL binding in PLASMA (QSS).
# Read this block first.
#
#   CHANGES vs R/13_whole_body_elimination.R:
#   1) SC depot: 60 mg (408.2 nmol, MW 147 kDa) into "Depot", bioavailability
#      F = 0.64, first-order absorption ka = 0.00883 /h into Central.
#      Source: Sutjandra et al. 2011 (Clin Pharmacokinet 50:793-807).
#   2) RANKL target in plasma with the quasi-steady-state (QSS) TMDD
#      approximation. "Central" now holds TOTAL drug in plasma (free +
#      RANKL-bound). Free concentration:
#        C = 0.5*[(Ctot - Rtot - Kss) + sqrt((Ctot - Rtot - Kss)^2 + 4*Kss*Ctot)]
#      Only FREE drug leaves plasma (to tissues, and by linear CL).
#      Target-mediated elimination: kint * RC * V_central, RC = Ctot - C.
#      Target: d(Rtot)/dt = ksyn - kdeg*Rtot - (kint - kdeg)*RC, Rtot(0) = R0.
#      Target parameters from Choi et al. 2025 (Front Pharmacol 16:1631034,
#      Table 3, CC-BY): Kss 1.56 nM, R0 15.23 nM, kint 0.022 /h,
#      ksyn 0.01 nM/h, kdeg = ksyn/R0 (derived, as in Paper 1).
#      Why Choi's target parameters may be used with Sutjandra's CL: Choi's
#      CL and volumes are apparent (divided by F), but concentrations, and
#      therefore Kss, R0, kint, kdeg, are not distorted by F.
#
#   ASSUMPTIONS (stated, to be tested later, not hidden):
#   - RANKL sits in plasma only (version 2a). Version 2b will test bone.
#   - Plasma RANKL concentration equals Choi's R0, although our plasma
#     volume (3.126 L) differs from Choi's central volume. A "same target
#     AMOUNT" alternative is a later sensitivity analysis.
#   - The drug-RANKL complex stays in plasma.
#   - Linear CL acts on FREE drug (as in the popPK models).
#
#   EXPECTED RESULTS (independent Python implementation, LSODA, rtol 1e-10):
#     Check 0  no drug: Rtot stays at R0 (drift ~0)
#     Check 1  drug mass balance (depot + body + linear elim + RANKL elim
#              = F*dose): rel. error ~1e-11
#     Check 2  target switched off: AUC*CL/(F*dose) = 1.0000
#     Total Cmax 4.48 ug/mL at day 12.5; free Cmax 4.40 ug/mL at day 12.7
#     Apparent half-life (total, day 30-120): 24.1 days
#     Total conc (ug/mL): day 1 1.662, day 7 4.045, day 10 4.403,
#                         day 28 3.643, day 84 1.061, day 168 0.001
#     Rtot (nM): day 1 12.54, day 7 1.19, day 28 0.48, day 84 0.55, day 168 7.23
#     By day 400: 54.1% eliminated via RANKL, 45.9% via linear CL
#   Prolia label, 60 mg SC: Cmax 6.75 ug/mL, Tmax 10 days, t1/2 25.4 days.
#   If your R output differs from the Python values beyond rounding, STOP.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/05_full_15_tissues.R")   # builds `lines`, `organ_params`, `state_cols`

cat("\n\n=======================================================================\n")
cat("=== R/14_sc_dosing_rankl_plasma.R - Paper 2, Step 2a ===\n")
cat("=======================================================================\n\n")

# ---------------------------------------------------------------------------
# Parameters (units: nmol, L, h)
# ---------------------------------------------------------------------------

MW        <- 147000                 # g/mol, denosumab (as in Paper 1)
dose_nmol <- 60e-3 / MW * 1e9       # 60 mg -> 408.2 nmol
CL_lin    <- 3.06 / 1000            # L/h   Sutjandra 2011
F_sc      <- 0.64                   # -     Sutjandra 2011
ka        <- 0.00883                # 1/h   Sutjandra 2011
Kss       <- 1.56                   # nM    Choi 2025 Table 3
R0        <- 15.23                  # nM    Choi 2025 Table 3
kint      <- 0.022                  # 1/h   Choi 2025 Table 3
ksyn      <- 0.01                   # nM/h  Choi 2025 Table 3
kdeg      <- ksyn / R0              # 1/h   derived (steady state)

ug_per_mL <- MW * 1e-6              # 1 nmol/L = 0.147 ug/mL

# ---------------------------------------------------------------------------
# Build the model text from file 05's `lines`.
# C_central is redefined as the FREE plasma concentration, so every flow
# that leaves plasma (to tissues) automatically uses free drug.
# ---------------------------------------------------------------------------

i_cc <- grep("^C_central <-", lines)
i_ce <- grep("^d/dt\\(Central\\)", lines)
stopifnot(length(i_cc) == 1, length(i_ce) == 1)

new_lines <- lines
new_lines[i_cc] <- paste(
  "Ctot <- Central / V_central",
  "bq <- Ctot - Rtot - Kss",
  "C_central <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot))",
  "RC <- Ctot - C_central",
  sep = "\n")
new_lines[i_ce] <- paste(lines[i_ce],
  "- CL_lin * C_central + ka * Depot - kint * RC * V_central")
new_lines <- c(
  "Rtot(0) <- R0",
  "f(Depot) <- F_sc",
  new_lines,
  "d/dt(Depot) <- -ka * Depot",
  "d/dt(Rtot) <- ksyn - kdeg * Rtot - (kint - kdeg) * RC",
  "d/dt(Eliminated) <- CL_lin * C_central",
  "d/dt(TMDD_elim) <- kint * RC * V_central")

mod2 <- rxode2(paste(new_lines, collapse = "\n"))

p_on  <- c(organ_params, CL_lin = CL_lin, F_sc = F_sc, ka = ka, Kss = Kss,
           R0 = R0, kint = kint, ksyn = ksyn, kdeg = kdeg)
p_off <- p_on; p_off["R0"] <- 0; p_off["ksyn"] <- 0   # target switched off

drug_cols <- c(state_cols, "Depot", "Eliminated", "TMDD_elim")

# ---------------------------------------------------------------------------
# Check 0: no drug -> Rtot stays at R0
# ---------------------------------------------------------------------------

ev0 <- et(seq(0, 5000, by = 50))
r0  <- rxSolve(mod2, p_on, ev0, atol = 1e-12, rtol = 1e-10)
cat("=== Check 0: no drug, Rtot drift ===\n")
cat("Max |Rtot - R0|:", signif(max(abs(r0$Rtot - R0)), 3), "nM (must be ~0)\n\n")

# ---------------------------------------------------------------------------
# Main run: 60 mg SC
# ---------------------------------------------------------------------------

ev <- et(amt = dose_nmol, cmt = "Depot", time = 0) |> et(seq(0, 24 * 400, by = 1))
r  <- rxSolve(mod2, p_on, ev, atol = 1e-12, rtol = 1e-10)

drug_total <- rowSums(r[, drug_cols])
mb_err <- max(abs(drug_total - F_sc * dose_nmol)) / (F_sc * dose_nmol)
cat("=== Check 1: drug mass balance (= F * dose) ===\n")
cat("Max relative error:", signif(mb_err, 3), "  Passed (< 1e-8)?", mb_err < 1e-8, "\n\n")

# ---------------------------------------------------------------------------
# Check 2: target switched off -> linear SC; AUC*CL/(F*dose) must be 1
# ---------------------------------------------------------------------------

ev_l <- et(amt = dose_nmol, cmt = "Depot", time = 0) |> et(seq(0, 24 * 2000, by = 1))
rl   <- rxSolve(mod2, p_off, ev_l, atol = 1e-12, rtol = 1e-10)
Cl   <- rl$Central / V_central
auc_l <- sum(diff(rl$time) * (head(Cl, -1) + tail(Cl, -1)) / 2)
cat("=== Check 2: target off, AUC*CL/(F*dose) ===\n")
cat("Ratio:", round(auc_l * CL_lin / (F_sc * dose_nmol), 4), "  [Python: 1.0000]\n\n")

# ---------------------------------------------------------------------------
# Results vs label
# ---------------------------------------------------------------------------

tt    <- r$time
Ctot_r <- r$Central / V_central
bq_r   <- Ctot_r - r$Rtot - Kss
Cfree  <- 0.5 * (bq_r + sqrt(bq_r^2 + 4 * Kss * Ctot_r))

cat("=== Results (60 mg SC) ===\n")
i1 <- which.max(Ctot_r); i2 <- which.max(Cfree)
cat(sprintf("Total Cmax %.2f ug/mL at day %.1f   [Python: 4.48, day 12.5]\n",
            Ctot_r[i1] * ug_per_mL, tt[i1] / 24))
cat(sprintf("Free  Cmax %.2f ug/mL at day %.1f   [Python: 4.40, day 12.7]\n",
            Cfree[i2] * ug_per_mL, tt[i2] / 24))

m <- tt >= 24 * 30 & tt <= 24 * 120
sl <- coef(lm(log(Ctot_r[m]) ~ tt[m]))[2]
cat(sprintf("Apparent half-life (total, day 30-120): %.1f days   [Python: 24.1]\n",
            log(2) / -unname(sl) / 24))
cat("Prolia label, 60 mg SC: Cmax 6.75 ug/mL, Tmax 10 days, t1/2 25.4 days\n\n")

for (d in c(1, 7, 10, 28, 84, 168)) {
  k <- which(tt == d * 24)
  cat(sprintf("Day %3d: total %.3f ug/mL  free %.4f  Rtot %.2f nM\n",
              d, Ctot_r[k] * ug_per_mL, Cfree[k] * ug_per_mL, r$Rtot[k]))
}
last2 <- r[nrow(r), ]
cat(sprintf("\nBy day 400: eliminated via RANKL %.1f%%, linear %.1f%%   [Python: 54.1 / 45.9]\n",
            100 * last2$TMDD_elim / (F_sc * dose_nmol),
            100 * last2$Eliminated / (F_sc * dose_nmol)))
