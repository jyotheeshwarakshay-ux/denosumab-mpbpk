# ===========================================================================
# SUMMARY - read this block first.
#
#   ONE VARIABLE CHANGED vs R/10_bone_fcrn_real_dose.R: the FcRn pool
#   denominator. 10 computed FcRn_Conc using Bone's own endosomal volume
#   alone (V_endosomal_Bone) - flagged there as overstating local FcRn
#   availability. This file instead sums V_endosomal[i] across all 15
#   tissues (matching the source's ARRAYSUM(V_endosomal[*]), source line
#   911) and uses THAT as the denominator - the whole-body-consistent
#   FcRn concentration. Bone-only scoping overstated local FcRn by
#   6.77-fold (FcRn_Conc 298 uM -> 44 uM once corrected). Everything else
#   - real 60 mg dose, IV-equivalent bolus, E6a-only rescue, single
#   binding, drug only, no endogenous IgG competition - is UNCHANGED from
#   10. The other 14 tissues still do not physically compete for FcRn in
#   this model (there's no FcRn state anywhere but Bone's E6a); only the
#   CONCENTRATION value now reflects what it would be if they did.
#
#   RESULT: Bone ABC drops from 7.28% (10, bone-only pool) to 5.86% (this
#   file) - landing just above the no-FcRn distribution baseline (5.70%).
#   This CONFIRMS the 7.28% near-match to the 7.27% benchmark in 10 was a
#   compensating-error artifact of the overstated bone-only FcRn pool,
#   exactly as 10's own header predicted before this file was run.
#
#   FINDING: at correctly-scarce whole-body FcRn, this E6a-only,
#   single-binding, drug-only rescue contributes only ~0.16 percentage
#   points to Bone ABC (5.86% vs 5.70% no-FcRn). This minimal-scope FcRn
#   mechanism does NOT explain the distribution model's residual gap to
#   the 7.27% benchmark - most of that gap remains unaccounted for.
#
#   5.86% IS AN UPPER BOUND for this scope, not a settled number: adding
#   the still-missing endogenous IgG competition would push it LOWER
#   (real IgG competes for the same now-properly-scarce FcRn pool, which
#   this file still lets the drug access uncontested).
#
#   STILL CARRIES: E6a-only rescue, single binding, drug only / no EDG
#   competition, IV-equivalent bolus (SC absorption deferred), single
#   dose (Q6-month regimen deferred).
#
#   Endothelial_Cell_Frac for all 15 tissues (source lines 217-231) is
#   needed to compute each tissue's endosomal volume - transcribed below,
#   verified to sum to 1.0000 across all 15 (a whole-body fractional
#   split, confirming these are the right 15 numbers to sum over).
# ===========================================================================
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.

library(rxode2)

source("R/01_parameters.R")
source("R/07_fcrn_parameters.R")

FR <- 0.715   # source line 140 [-]. Missing from 07_fcrn_parameters.R (as in 08/09/10).

# ---------------------------------------------------------------------------
# Endothelial_Cell_Frac, all 15 tissues (source lines 217-231), mapped to
# organ names via the index key confirmed earlier (Lung=1 ... Other=15).
# ---------------------------------------------------------------------------

Endothelial_Cell_Frac_15 <- c(
  Lung     = 0.0834,  # source line 217
  Liver    = 0.1877,  # source line 218
  Heart    = 0.0011,  # source line 219
  Muscle   = 0.1928,  # source line 220
  Skin     = 0.0819,  # source line 221
  Adipose  = 0.0999,  # source line 222
  Bone     = 0.1478,  # source line 223
  Brain    = 0.0115,  # source line 224
  Kidney   = 0.0157,  # source line 225
  SI       = 0.0121,  # source line 226
  LI       = 0.0209,  # source line 227
  Pancreas = 0.0001,  # source line 228
  Thymus   = 0.0001,  # source line 229
  Spleen   = 0.0499,  # source line 230
  Other    = 0.0951   # source line 231
)

stopifnot(abs(sum(Endothelial_Cell_Frac_15) - 1) < 1e-9)

# ---------------------------------------------------------------------------
# V_endosomal[i] for all 15 tissues (source lines 908, 910), and their sum
# ---------------------------------------------------------------------------

Organ_Endothelial_Cell_15 <- Total_Endothelial_Cell * Endothelial_Cell_Frac_15 * Scale_Factor
V_endosomal_15 <- CL_Up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell_15 * 1e-6 * 1e-9

cat("=== V_endosomal[i], all 15 tissues ===\n")
print(format(V_endosomal_15, digits = 6, scientific = FALSE))
cat("Sum:", format(sum(V_endosomal_15), digits = 6), "\n\n")

stopifnot(abs(unname(V_endosomal_15["Bone"]) - V_endosomal_Bone) < 1e-12)  # sanity vs 07

# ---------------------------------------------------------------------------
# FcRn_Conc: old (bone-only denominator, from 07/09/10) vs new (whole-body)
# ---------------------------------------------------------------------------

FcRn_Conc_BoneOnly <- FcRn_Availability * total_FcRn_in_nmole * 1e-3 / V_endosomal_Bone
FcRn_Conc_Shared    <- FcRn_Availability * total_FcRn_in_nmole * 1e-3 / sum(V_endosomal_15)

cat("=== FcRn_Conc: bone-only denominator (10) vs whole-body denominator (this file) ===\n")
cat("FcRn_Conc, bone-only (old):  ", format(FcRn_Conc_BoneOnly, digits = 6), "uM\n")
cat("FcRn_Conc, whole-body (new): ", format(FcRn_Conc_Shared, digits = 6), "uM\n")
cat("Ratio (old / new):           ", format(FcRn_Conc_BoneOnly / FcRn_Conc_Shared, digits = 6),
    "-> bone-only overstated local FcRn by this factor\n\n")

FcRn_E6a_init_amount <- FcRn_Conc_Shared * V_E6a_Bone   # only this line differs functionally from 10

# ---------------------------------------------------------------------------
# Real dose conversion - unchanged from 10
# ---------------------------------------------------------------------------

dose_mg      <- 60
MW_denosumab <- 147000
V_plasma_L   <- unname(V_Plasma)

C0_plasma_uM <- dose_mg / (V_plasma_L * 1000) / MW_denosumab * 1e6

# ---------------------------------------------------------------------------
# FcRn binding constants for E6a - unchanged from 10
# ---------------------------------------------------------------------------

k_off_6_EXG <- (KD6 / 1000) * k_on_6_EXG
stopifnot(abs(k_off_6_EXG - 56.42) < 1e-6)

# ---------------------------------------------------------------------------
# Organ groups (identical to 05/08/09/10)
# ---------------------------------------------------------------------------

generic_organs <- c("Lung", "Heart", "Muscle", "Skin", "Adipose", "Bone",
                     "Brain", "Kidney", "SI", "LI", "Pancreas", "Thymus",
                     "Spleen", "Other")
gut_organs   <- c("Spleen", "Pancreas", "SI", "LI")
plain_organs <- setdiff(generic_organs, gut_organs)
all_organs   <- c(generic_organs, "Liver")
non_bone_plain_organs <- setdiff(plain_organs, "Bone")

V_central   <- unname(V_Plasma)
V_ln        <- unname(V_LN)
L_LymphNode <- sum(LF[all_organs])

organ_params <- c(V_central = V_central, V_ln = V_ln, L_LymphNode = L_LymphNode,
                   FR = FR, Prob_deg = Prob_deg,
                   CL_UP_Bone = unname(CL_UP_Bone),
                   V_VM_Bone  = unname(V_VM_Bone),
                   V_E7_Bone  = unname(V_E7_Bone),
                   V_E6a_Bone = unname(V_E6a_Bone),
                   V_E7b_Bone = unname(V_E7b_Bone),
                   V_ISM_Bone = unname(V_ISM_Bone),
                   k_on_6_EXG  = unname(k_on_6_EXG),
                   k_off_6_EXG = unname(k_off_6_EXG))

for (o in all_organs) {
  organ_params[paste0("PLQ_", o)]     <- unname(PLQ[o])
  organ_params[paste0("LF_", o)]      <- unname(LF[o])
  organ_params[paste0("sigmaV_", o)]  <- unname(sigma_V[o])
  organ_params[paste0("sigmaIS_", o)] <- unname(sigma_IS[o])
  organ_params[paste0("Vvasc_", o)]   <- unname(V_V[o])
  organ_params[paste0("Vinter_", o)]  <- unname(V_IntS[o])
}

# ---------------------------------------------------------------------------
# Build the rxode2 model text - byte-for-byte identical to 10's model
# ---------------------------------------------------------------------------

lines <- c("C_central <- Central / V_central", "C_ln <- LymphNode / V_ln")

for (o in all_organs) {
  lines <- c(lines,
    sprintf("C_%sV  <- %sV  / Vvasc_%s",  o, o, o),
    sprintf("C_%sIS <- %sIS / Vinter_%s", o, o, o))
}

lines <- c(lines,
  "C_BoneVM     <- BoneVM     / V_VM_Bone",
  "C_BoneE7     <- BoneE7     / V_E7_Bone",
  "C_BoneE6a    <- BoneE6a    / V_E6a_Bone",
  "C_BoneE7b    <- BoneE7b    / V_E7b_Bone",
  "C_BoneISM    <- BoneISM    / V_ISM_Bone",
  "C_BoneFcRn   <- BoneFcRn_E6a  / V_E6a_Bone",
  "C_BoneBound  <- BoneBound_E6a / V_E6a_Bone",
  "on_rate_E6a  <- k_on_6_EXG  * C_BoneE6a * C_BoneFcRn * V_E6a_Bone",
  "off_rate_E6a <- k_off_6_EXG * C_BoneBound * V_E6a_Bone")

central_terms <- character(0)
for (o in plain_organs) {
  central_terms <- c(central_terms,
    sprintf("- PLQ_%s * C_central + (PLQ_%s - LF_%s) * C_%sV", o, o, o, o))
}
for (g in gut_organs) {
  central_terms <- c(central_terms,
    sprintf("- PLQ_%s * C_central + (PLQ_%s - LF_%s) * C_LiverV", g, g, g))
}
central_terms <- c(central_terms,
  "- PLQ_Liver * C_central + (PLQ_Liver - LF_Liver) * C_LiverV")
lines <- c(lines, sprintf("d/dt(Central) <- %s + L_LymphNode * C_ln",
                           paste(central_terms, collapse = " ")))

ln_terms <- character(0)
for (o in all_organs) {
  ln_terms <- c(ln_terms, sprintf("+ LF_%s * (1 - sigmaIS_%s) * C_%sIS", o, o, o))
}
lines <- c(lines, sprintf("d/dt(LymphNode) <- %s - L_LymphNode * C_ln",
                           paste(ln_terms, collapse = " ")))

for (o in non_bone_plain_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            o, o, o, o, o, o, o, o),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            o, o, o, o, o, o, o))
}

for (g in gut_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            g, g, g, g, g, g, g, g),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            g, g, g, g, g, g, g))
}

portal_inflow_terms <- sapply(gut_organs, function(g)
  sprintf("+ (PLQ_%s - LF_%s) * C_%sV", g, g, g))
outlet_terms <- c("(PLQ_Liver - LF_Liver)",
                   sapply(gut_organs, function(g) sprintf("(PLQ_%s - LF_%s)", g, g)))
lines <- c(lines,
  sprintf("d/dt(LiverV) <- PLQ_Liver * C_central %s - C_LiverV * (%s) - LF_Liver * (1 - sigmaV_Liver) * C_LiverV",
          paste(portal_inflow_terms, collapse = " "), paste(outlet_terms, collapse = " + ")),
  "d/dt(LiverIS) <- LF_Liver * (1 - sigmaV_Liver) * C_LiverV - LF_Liver * (1 - sigmaIS_Liver) * C_LiverIS")

lines <- c(lines,
  "d/dt(BoneV)  <- PLQ_Bone * C_central - (PLQ_Bone - LF_Bone) * C_BoneV - LF_Bone * (1 - sigmaV_Bone) * C_BoneV - CL_UP_Bone * FR * C_BoneV + CL_UP_Bone * FR * C_BoneVM",
  "d/dt(BoneIS) <- LF_Bone * (1 - sigmaV_Bone) * C_BoneV - LF_Bone * (1 - sigmaIS_Bone) * C_BoneIS - CL_UP_Bone * (1 - FR) * C_BoneIS + CL_UP_Bone * (1 - FR) * C_BoneISM"
)

lines <- c(lines,
  "d/dt(BoneVM) <- CL_UP_Bone * FR * C_BoneV - CL_UP_Bone * FR * C_BoneVM - CL_UP_Bone * FR * C_BoneVM + CL_UP_Bone * FR * C_BoneE7b",
  "d/dt(BoneE7) <- CL_UP_Bone * FR * C_BoneVM + CL_UP_Bone * (1 - FR) * C_BoneISM - CL_UP_Bone * C_BoneE7",
  "d/dt(BoneE6a) <- CL_UP_Bone * C_BoneE7 - CL_UP_Bone * C_BoneE6a - on_rate_E6a + off_rate_E6a",
  "d/dt(BoneFcRn_E6a) <- - on_rate_E6a + off_rate_E6a",
  "d/dt(BoneBound_E6a) <- on_rate_E6a - off_rate_E6a - CL_UP_Bone * C_BoneBound * V_E6a_Bone",
  "d/dt(BoneE7b) <- CL_UP_Bone * (1 - Prob_deg) * C_BoneE6a - CL_UP_Bone * C_BoneE7b + CL_UP_Bone * C_BoneBound * V_E6a_Bone",
  "d/dt(BoneISM) <- CL_UP_Bone * (1 - FR) * C_BoneIS - CL_UP_Bone * (1 - FR) * C_BoneISM + CL_UP_Bone * (1 - FR) * C_BoneE7b - CL_UP_Bone * (1 - FR) * C_BoneISM",
  "d/dt(BoneDegraded) <- CL_UP_Bone * Prob_deg * C_BoneE6a"
)

model_txt <- paste(lines, collapse = "\n")

mod <- rxode2(model_txt)

stopifnot(length(mod$state) == 40)

# ---------------------------------------------------------------------------
# Dose: real IV-equivalent bolus, unchanged from 10. Simulate 0-20000 h.
# ---------------------------------------------------------------------------

dose_amt <- C0_plasma_uM * V_central

ev <- et(amt = dose_amt, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 10))

res <- rxSolve(mod, organ_params, ev,
                inits = c(BoneFcRn_E6a = unname(FcRn_E6a_init_amount)))

# ---------------------------------------------------------------------------
# Conservation check
# ---------------------------------------------------------------------------

live_cols <- c("Central", "LymphNode",
               as.vector(rbind(paste0(all_organs, "V"), paste0(all_organs, "IS"))),
               "BoneVM", "BoneE7", "BoneE6a", "BoneE7b", "BoneISM",
               "BoneBound_E6a")

total_live <- rowSums(res[, live_cols])
total_all  <- total_live + res$BoneDegraded

cat("=== Conservation check: total_live + BoneDegraded vs dose_amt ===\n")
cat("dose_amt:", format(dose_amt, digits = 6), "\n")
cat("min: ", format(min(total_all), digits = 6), "\n")
cat("max: ", format(max(total_all), digits = 6), "\n")
cat("Relative range:", (max(total_all) - min(total_all)) / mean(total_all), "\n")
cat("Conserved to dose_amt throughout (rel. range < 1e-8)?",
    (max(total_all) - min(total_all)) / mean(total_all) < 1e-8, "\n\n")

# ---------------------------------------------------------------------------
# Key result: four-way ABC comparison
# ---------------------------------------------------------------------------

last <- res[nrow(res), ]
C_central_last <- last$Central / V_central

bone_drug <- last$BoneV + last$BoneIS + last$BoneVM + last$BoneE7 +
             last$BoneE6a + last$BoneE7b + last$BoneISM + last$BoneBound_E6a

abc_shared_pool_pct <- 100 * (bone_drug / unname(V_Organ["Bone"])) / C_central_last

abc_bone_only_pct <- 7.28
no_fcrn_pct       <- 5.70
benchmark_pct     <- 7.27

cat("=== Key result: Bone ABC, shared (whole-body) FcRn pool vs references ===\n")
cat("Bone ABC, shared FcRn pool (this file, 11):   ", round(abc_shared_pool_pct, 4), "%\n")
cat("Bone ABC, bone-only FcRn pool (10):           ", abc_bone_only_pct, "%\n")
cat("No-FcRn distribution-only (06):               ", no_fcrn_pct, "%\n")
cat("Published benchmark (Shah & Betts 2013):      ", benchmark_pct, "%\n")
