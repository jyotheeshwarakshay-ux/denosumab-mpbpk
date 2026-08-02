# ===========================================================================
# FINDING - read this block first. Diagnostic calculation only: no model
# build, no new ODEs. Builds on R/11_bone_fcrn_shared_pool.R's steady-state
# output.
#
#   QUESTION TESTED: could bone's missing bulk CELLULAR volume explain the
#   ~1.4 percentage-point gap between the current model's Bone ABC (5.86%,
#   file 11) and the published benchmark (7.27%)? File 11 counts bone drug
#   in vascular + interstitial + endosomal/membrane compartments only - it
#   has no bulk cellular compartment holding drug, even though bone's
#   cellular volume (V_C["Bone"] = 7.817 L) is over 4x its interstitial
#   volume (V_IntS["Bone"] = 1.891 L, both from R/01_parameters.R).
#
#   BOUNDING CALCULATION: added a hypothetical amount of drug to bone's
#   cellular volume under two assumptions, and recomputed ABC each time -
#     (a) LOWER bound - cellular concentration = interstitial concentration
#         (added drug = C_BoneIS * V_C[Bone])       -> ABC = 20.31%
#     (b) UPPER bound - cellular concentration = plasma concentration
#         (added drug = C_central * V_C[Bone])      -> ABC = 82.76%
#   Both BOUNDS massively overshoot the 7.27% benchmark - the lower bound
#   alone overshoots by ~3x.
#
#   CONCLUSION: the gap is NOT a missing-volume problem. Bone's cellular
#   volume is more than large enough to close a 1.4-point gap; the actual
#   limiting factor is that real intracellular antibody concentration is
#   kinetically limited (receptor/pinocytosis-driven uptake), far below
#   bulk interstitial or plasma levels - filling the cellular volume at
#   either of those concentrations is not physiologically plausible. This
#   is consistent with, not a coincidence relative to, how the source
#   itself is built: it models "cellular" content specifically AS the
#   FcRn endosomal/membrane machinery (VM, E7, E6a, E7b, ISM) - exactly
#   the compartments already being built in this project - rather than as
#   a separate bulk-equilibrium pool.
#
#   THEREFORE: the residual gap to 7.27% is an FcRn-KINETICS question
#   (total FcRn amount, binding/recycling rates, multi-site cooperative
#   binding, endogenous IgG competition), not a volume question.
#   Hypothesis "missing bulk cellular volume explains the gap" is RULED
#   OUT by this calculation.
# ===========================================================================
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.

source("R/11_bone_fcrn_shared_pool.R")

last <- res[nrow(res), ]
C_central_last <- last$Central / V_central
C_BoneIS_last  <- last$BoneIS / unname(V_IntS["Bone"])

bone_drug_current <- last$BoneV + last$BoneIS + last$BoneVM + last$BoneE7 +
                      last$BoneE6a + last$BoneE7b + last$BoneISM + last$BoneBound_E6a

abc_current_pct <- 100 * (bone_drug_current / unname(V_Organ["Bone"])) / C_central_last

V_C_Bone    <- unname(V_C["Bone"])
V_IntS_Bone <- unname(V_IntS["Bone"])

cat("\n\n=== Diagnostic: cellular-volume bounding calculation ===\n\n")
cat("1. Bone total drug (current, counted compartments):", format(bone_drug_current, digits = 6), "\n")
cat("   Bone ABC (current):", round(abc_current_pct, 4), "%\n\n")

cat("2. C_BoneIS (steady state):", format(C_BoneIS_last, digits = 6), "uM\n")
cat("   C_central (steady state):", format(C_central_last, digits = 6), "uM\n\n")

cat("3. V_C[Bone] (cellular volume):", V_C_Bone, "L\n")
cat("   V_IntS[Bone] (interstitial volume):", V_IntS_Bone, "L\n\n")

# (a) LOWER bound: cellular concentration = interstitial concentration
added_drug_lower <- C_BoneIS_last * V_C_Bone
bone_drug_lower  <- bone_drug_current + added_drug_lower
abc_lower_pct    <- 100 * (bone_drug_lower / unname(V_Organ["Bone"])) / C_central_last

# (b) UPPER bound: cellular concentration = plasma concentration
added_drug_upper <- C_central_last * V_C_Bone
bone_drug_upper  <- bone_drug_current + added_drug_upper
abc_upper_pct    <- 100 * (bone_drug_upper / unname(V_Organ["Bone"])) / C_central_last

cat("=== Bounding estimates ===\n")
cat(sprintf("(a) LOWER bound (cellular conc = interstitial conc):\n"))
cat(sprintf("    Added drug = C_BoneIS * V_C[Bone] = %.6g * %.6g = %.6g\n",
            C_BoneIS_last, V_C_Bone, added_drug_lower))
cat(sprintf("    New total bone drug = %.6g\n", bone_drug_lower))
cat(sprintf("    ABC (lower bound) = %.4f %%\n\n", abc_lower_pct))

cat(sprintf("(b) UPPER bound (cellular conc = plasma conc):\n"))
cat(sprintf("    Added drug = C_central * V_C[Bone] = %.6g * %.6g = %.6g\n",
            C_central_last, V_C_Bone, added_drug_upper))
cat(sprintf("    New total bone drug = %.6g\n", bone_drug_upper))
cat(sprintf("    ABC (upper bound) = %.4f %%\n\n", abc_upper_pct))

cat("=== Comparison ===\n")
cat(sprintf("Current (no cellular compartment):    %.4f %%\n", abc_current_pct))
cat(sprintf("Lower bound (cellular = IS conc):     %.4f %%\n", abc_lower_pct))
cat(sprintf("Upper bound (cellular = plasma conc): %.4f %%\n", abc_upper_pct))
cat(sprintf("Benchmark:                              7.27 %%\n\n"))

cat("=== Conclusion ===\n")
cat("Both bounds overshoot the 7.27% benchmark (lower bound overshoots by\n")
cat("~3x). Bulk cellular volume is RULED OUT as the gap explanation - the\n")
cat("volume is plenty large enough; what's missing is a physiological\n")
cat("reason for drug to reach anywhere near bulk concentration inside\n")
cat("cells, and that reason is FcRn kinetics, not raw space.\n")
