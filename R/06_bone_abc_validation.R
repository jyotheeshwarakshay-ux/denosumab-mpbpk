# Antibody biodistribution coefficient (ABC) validation, built on top of
# R/05_full_15_tissues.R (sourced directly below - reuses its 15-tissue,
# liver-portal-routed model, its steady-state solve, and all its objects:
# res, V_central, organ_params, all_organs, conc_last(), last).
#
# ABC definition, matching the source's own readout (PSP4-8-738-s002.txt
# line 982):
#   mAb_Biodistr_Coeff_in_tissue[i] = (Total_EXG_intra_Tissue[i]/V_Organ[i])
#                                      / C_EXG_Plasma
# i.e. total intra-tissue antibody amount, divided by that tissue's whole
# organ volume, divided by plasma concentration.
#
# KNOWN SIMPLIFICATION: the source's Total_EXG_intra_Tissue sums vascular +
# vascular membrane + THREE endosomal sub-compartments + interstitial
# membrane + interstitial space (source lines 972-980) - i.e. it includes
# the full endothelial/endosomal (cellular) machinery, because that model
# has FcRn. Our distribution skeleton (05) has no cellular/endothelial
# compartment at all - only vascular and interstitial. So here,
# "total intra-tissue amount" = vascular amount + interstitial amount ONLY.
# This is expected to sit BELOW the full model's ABC, not a bug: we are
# missing whatever antibody the full model holds in the endothelial/
# endosomal space. The point of this file is to see how much of the 7.27%
# benchmark a pure convective distribution skeleton can explain on its own,
# before FcRn is added.
#
# Distribution skeleton only - NO FcRn, NO cellular compartments added here.
#
# Run from the project root (denosumab-mpbpk/) so source() calls resolve.

source("R/05_full_15_tissues.R")

cat("\n\n")
cat("=======================================================================\n")
cat("=== R/06_bone_abc_validation.R - ABC computation (builds on above) ===\n")
cat("=======================================================================\n\n")

# ---------------------------------------------------------------------------
# ABC[i] = (VascularAmt_i + InterstitialAmt_i) / V_Organ[i] / C_central
# at the final (steady-state) timepoint. V_Organ comes from
# R/01_parameters.R (already in scope via 05's source() chain).
# ---------------------------------------------------------------------------

abc <- sapply(all_organs, function(o) {
  intra_tissue_amt <- last[[paste0(o, "V")]] + last[[paste0(o, "IS")]]
  (intra_tissue_amt / unname(V_Organ[o])) / C_central_last
})

abc_pct <- abc * 100
abc_pct_sorted <- sort(abc_pct)

# ---------------------------------------------------------------------------
# Step 3: all 15 tissue ABCs, sorted, as % of plasma
# ---------------------------------------------------------------------------

cat("=== Step 3: ABC per tissue (% of plasma concentration), sorted, at t =",
    last$time, "h ===\n")
print(round(abc_pct_sorted, 3))

# ---------------------------------------------------------------------------
# Step 4: Bone ABC vs the Shah & Betts 2013 benchmark (7.27%)
# ---------------------------------------------------------------------------

bone_abc_pct <- unname(abc_pct["Bone"])
benchmark_pct <- 7.27
fold_diff <- bone_abc_pct / benchmark_pct

cat("\n=== Step 4: Bone ABC vs benchmark ===\n")
cat("Bone ABC (this distribution-only skeleton): ", round(bone_abc_pct, 3), "%\n")
cat("Published benchmark (Shah & Betts 2013):    ", benchmark_pct, "%\n")
cat("Ratio (this model / benchmark):             ", round(fold_diff, 3), "\n")
if (fold_diff < 1) {
  cat("-> This model's Bone ABC is", round(1 / fold_diff, 2),
      "-fold LOWER than the benchmark.\n")
} else {
  cat("-> This model's Bone ABC is", round(fold_diff, 2),
      "-fold HIGHER than the benchmark.\n")
}

# ---------------------------------------------------------------------------
# Step 5: interpretation - is Bone within a plausible (~2-fold) range of
# the benchmark? 2-fold is the accuracy band Shah & Betts themselves report
# for the validation of this ABC approach across species.
# ---------------------------------------------------------------------------

lower_bound <- benchmark_pct / 2
upper_bound <- benchmark_pct * 2
within_2fold <- bone_abc_pct >= lower_bound && bone_abc_pct <= upper_bound

cat("\n=== Step 5: interpretation ===\n")
cat("Bone ABC (distribution skeleton only):", round(bone_abc_pct, 3), "%\n")
cat("Benchmark (Shah & Betts 2013):        ", benchmark_pct, "%\n")
cat("2-fold plausibility band:             [", round(lower_bound, 3), "%, ",
    round(upper_bound, 3), "%]\n", sep = "")
cat("Within 2-fold of benchmark?", within_2fold, "\n")
if (within_2fold) {
  cat("-> Distribution (vascular + interstitial) alone already accounts for\n",
      "   Bone's ABC within the accuracy band the source itself reports.\n",
      "   Any remaining gap to the full model's number is attributable to\n",
      "   the missing endothelial/endosomal (cellular) compartment and FcRn,\n",
      "   which this skeleton deliberately does not include yet.\n", sep = "")
} else if (bone_abc_pct < lower_bound) {
  cat("-> Distribution alone falls short of the benchmark by more than\n",
      "   2-fold. Some of the gap is expected (missing cellular fraction),\n",
      "   but a shortfall this large is worth re-checking once FcRn is\n",
      "   added, in case a parameter or term is off rather than just\n",
      "   'the cellular compartment isn't here yet'.\n", sep = "")
} else {
  cat("-> Distribution alone already EXCEEDS the benchmark by more than\n",
      "   2-fold, which is unexpected for a skeleton that is missing a\n",
      "   compartment (it should if anything run low, not high) - worth\n",
      "   double-checking the Bone-specific parameters/terms before adding\n",
      "   FcRn on top of this.\n", sep = "")
}
