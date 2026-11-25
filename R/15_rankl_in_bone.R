# ===========================================================================
# PAPER 2 - STEP 2b: RANKL in BONE interstitium instead of plasma.
# Read this block first.
#
#   QUESTION: if the same total amount of RANKL sits in bone interstitial
#   fluid instead of plasma, does the plasma denosumab profile change? If
#   not, plasma data cannot tell where the target is. If yes, plasma data
#   carry information about target location.
#
#   CHANGES vs R/14_sc_dosing_rankl_plasma.R:
#   - Plasma has NO target: Central is linear again (CL_lin on plasma conc).
#   - Target sits in BoneIS, QSS: BoneIS holds total drug (free + bound);
#     free bone conc C_BoneIS from the QSS formula; only free drug leaves
#     BoneIS by lymph; target-mediated loss kint * RC * Vinter_Bone.
#   - Same total target AMOUNT as 2a: R0_bone = R0 * V_plasma / V_BoneIS
#     = 15.23 * 3.126 / 1.891 = 25.18 nM; ksyn_bone = kdeg * R0_bone.
#     Kss, kint, kdeg unchanged (Choi 2025).
#   - SC dosing unchanged (F 0.64, ka 0.00883 /h, Sutjandra 2011).
#
#   EXPECTED RESULTS (independent Python implementation, LSODA, rtol 1e-10):
#     mass balance rel. error ~2e-11
#     plasma total Cmax 5.63 ug/mL at day 12.4
#     apparent half-life (day 30-120): 41.8 days
#     eliminated via RANKL by day 400: 20.1%
#     plasma total conc (ug/mL): day 7 5.146, day 28 4.680, day 84 1.854,
#                                day 140 0.732
#     free bone conc (ug/mL): day 7 0.0053, day 28 0.0078, day 84 0.0037
#     free RANKL in bone (fraction of baseline): day 7 0.942, day 28 0.709,
#                                                day 84 0.598, day 140 0.722
#   For comparison, 2a (RANKL in plasma): Cmax 4.48, half-life 24.1 days,
#   54.1% via RANKL. Prolia label: Cmax 6.75 ug/mL, t1/2 25.4 days.
#
#   INTERPRETATION LIMITS (do not overclaim):
#   - Bone delivery in this skeleton is known to be low (Bone ABC 5.70% vs
#     7.27% benchmark; no FcRn, no endothelial compartment), so weak bone
#     RANKL suppression here may partly be a delivery artefact.
#   - No published human time course of total RANKL after denosumab was
#     found; serum FREE RANKL is reported to fall (D-BIOMARK, Breast Cancer
#     Res 2025, doi 10.1186/s13058-025-01996-w). Choi 2025 states drug
#     concentrations alone could not estimate kdeg and kint.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/05_full_15_tissues.R")

cat("\n\n=======================================================================\n")
cat("=== R/15_rankl_in_bone.R - Paper 2, Step 2b ===\n")
cat("=======================================================================\n\n")

MW        <- 147000
dose_nmol <- 60e-3 / MW * 1e9
CL_lin    <- 3.06 / 1000            # Sutjandra 2011
F_sc      <- 0.64                   # Sutjandra 2011
ka        <- 0.00883                # Sutjandra 2011
Kss       <- 1.56                   # Choi 2025
kint      <- 0.022                  # Choi 2025
kdeg      <- 0.01 / 15.23           # Choi 2025, derived ksyn/R0
V_BIS     <- unname(organ_params["Vinter_Bone"])
R0_bone   <- 15.23 * V_central / V_BIS
ksyn_bone <- kdeg * R0_bone
ug_per_mL <- MW * 1e-6

cat(sprintf("R0 in bone interstitium: %.2f nM   [Python: 25.18]\n\n", R0_bone))

i_cb <- grep("^C_BoneIS", lines)
i_bi <- grep("^d/dt\\(BoneIS\\)", lines)
i_ce <- grep("^d/dt\\(Central\\)", lines)
stopifnot(length(i_cb) == 1, length(i_bi) == 1, length(i_ce) == 1)

nl <- lines
nl[i_cb] <- paste(
  "CtotB <- BoneIS / Vinter_Bone",
  "bqB <- CtotB - RtotB - Kss",
  "C_BoneIS <- 0.5 * (bqB + sqrt(bqB * bqB + 4 * Kss * CtotB))",
  "RCB <- CtotB - C_BoneIS",
  sep = "\n")
nl[i_bi] <- paste(lines[i_bi], "- kint * RCB * Vinter_Bone")
nl[i_ce] <- paste(lines[i_ce], "- CL_lin * C_central + ka * Depot")
nl <- c("RtotB(0) <- R0_bone", "f(Depot) <- F_sc", nl,
        "d/dt(Depot) <- -ka * Depot",
        "d/dt(RtotB) <- ksyn_bone - kdeg * RtotB - (kint - kdeg) * RCB",
        "d/dt(Eliminated) <- CL_lin * C_central",
        "d/dt(TMDD_elim) <- kint * RCB * Vinter_Bone")

mod3 <- rxode2(paste(nl, collapse = "\n"))
p3 <- c(organ_params, CL_lin = CL_lin, F_sc = F_sc, ka = ka, Kss = Kss,
        kint = kint, kdeg = kdeg, R0_bone = R0_bone, ksyn_bone = ksyn_bone)

ev <- et(amt = dose_nmol, cmt = "Depot", time = 0) |> et(seq(0, 24 * 400, by = 1))
r  <- rxSolve(mod3, p3, ev, atol = 1e-12, rtol = 1e-10)

drug_cols <- c(state_cols, "Depot", "Eliminated", "TMDD_elim")
mb <- max(abs(rowSums(r[, drug_cols]) - F_sc * dose_nmol)) / (F_sc * dose_nmol)
cat("=== Check 1: drug mass balance ===\n")
cat("Max relative error:", signif(mb, 3), "  Passed (< 1e-8)?", mb < 1e-8, "\n\n")

tt   <- r$time
Cpl  <- r$Central / V_central
CtB  <- r$BoneIS / V_BIS
bq   <- CtB - r$RtotB - Kss
CfB  <- 0.5 * (bq + sqrt(bq^2 + 4 * Kss * CtB))
fR   <- (r$RtotB - r$RtotB * CfB / (Kss + CfB)) / R0_bone

k1 <- which.max(Cpl)
m  <- tt >= 24 * 30 & tt <= 24 * 120
sl <- coef(lm(log(Cpl[m]) ~ tt[m]))[2]
last3 <- r[nrow(r), ]
cat("=== Results (60 mg SC, RANKL in bone) ===\n")
cat(sprintf("Plasma Cmax %.2f ug/mL at day %.1f   [Python: 5.63, day 12.4]\n",
            Cpl[k1] * ug_per_mL, tt[k1] / 24))
cat(sprintf("Apparent half-life (day 30-120): %.1f days   [Python: 41.8]\n",
            log(2) / -unname(sl) / 24))
cat(sprintf("Eliminated via RANKL by day 400: %.1f%%   [Python: 20.1]\n\n",
            100 * last3$TMDD_elim / (F_sc * dose_nmol)))
for (d in c(7, 28, 84, 140)) {
  k <- which(tt == d * 24)
  cat(sprintf("Day %3d: plasma %.3f ug/mL | free bone %.4f ug/mL | free RANKL in bone %.3f of baseline\n",
              d, Cpl[k] * ug_per_mL, CfB[k] * ug_per_mL, fR[k]))
}
cat("[Python: plasma 5.146 / 4.680 / 1.854 / 0.732; free bone 0.0053 / 0.0078 / 0.0037 / -;",
    "free RANKL 0.942 / 0.709 / 0.598 / 0.722]\n")
cat("Compare 2a (RANKL in plasma): Cmax 4.48, half-life 24.1 days, 54.1% via RANKL.\n")
