# ===========================================================================
# SUMMARY - read this block first.
#
#   ONE VARIABLE CHANGED vs R/09_bone_fcrn_rescue.R: the dimensionless
#   "1 unit" dose is replaced with a real 60 mg denosumab dose
#   (C0_plasma_uM = 0.1306 uM), given as an IV-equivalent bolus. Nothing
#   else in the model differs from 09.
#
#   RESULT: Bone ABC moves from 6.36% (09, arbitrary dose) to 7.28% (this
#   file) from dose realism ALONE - a ~14% shift. This shows FcRn rescue
#   is dose-sensitive at clinical denosumab concentrations: binding
#   saturation is non-negligible around ~0.13 uM, not just a theoretical
#   concern.
#
#   THE 7.28% NEAR-MATCH TO THE 7.27% BENCHMARK IS NOT VALIDATION. Two
#   known simplifications both bias the ABC UPWARD: (1) the bone-scoped
#   FcRn pool overstates local FcRn availability (see simplification 3
#   below and 09's header), and (2) there is no endogenous IgG (EDG)
#   competing for that same FcRn pool, so the drug gets uncontested
#   access to receptors a real cell would have to share. Both push the
#   number up. A correctly-built model (whole-body FcRn competition
#   across all 15 tissues, plus real EDG competition) would most likely
#   land LOWER than 7.28% - the closeness to 7.27% here is best read as
#   compensating errors landing in the same place, not quantitative
#   correctness.
#
#   STILL CARRIES every R/09 simplification, unchanged: E6a-only rescue,
#   single binding, drug only (no EDG), bone-scoped FcRn pool. Plus two
#   new ones specific to this file: IV-equivalent bolus (SC absorption
#   deferred) and single dose (Q6-month regimen deferred). See the full
#   list below.
# ===========================================================================
#
# Everything else - E6a-only rescue, single binding, drug only, no
# endogenous IgG competition, bone-scoped FcRn pool, all of 09's other
# documented simplifications - is UNCHANGED. See 09's header for the full
# list of what those mean; not repeated here except where the real dose
# changes the interpretation.
#
# DOCUMENTED SIMPLIFICATIONS specific to this file:
#
#   1. IV-EQUIVALENT BOLUS. Denosumab is clinically given SUBCUTANEOUSLY
#      (~60% bioavailability, first-order absorption from a depot - the
#      source itself models this via `ka`/`F`, source lines 37-38). This
#      file doses the full amount directly into Central at t=0, as if IV.
#      SC absorption (depot compartment, ka, F=0.60) is DEFERRED as a
#      separate future step, not implemented here.
#
#   2. SINGLE DOSE ONLY. Real denosumab dosing is Q6-month (repeat SC
#      injections). This file simulates one dose to steady-state decay,
#      not a multi-dose regimen. That's also deferred future work.
#
#   3. This resolves (not just documents) the "absolute dose scale"
#      caveat flagged in 09: because FcRn binding is saturable, dosing at
#      an arbitrary "1 unit" and dosing at a real ~0.13 uM plasma
#      concentration are NOT interchangeable - the whole point of this
#      file is to put the dose on the same real, physical (uM) scale that
#      FcRn_Conc_Bone (derived from a real 1022 nmol FcRn abundance,
#      R/07_fcrn_parameters.R) already uses, so the two are finally
#      comparable in absolute terms.
#
# Dose parameters (denosumab-specific, verified against FDA label - NOT
# the source .txt's generic values):
#   dose_mg       = 60      [mg]     denosumab clinical dose (flat, not
#                                     mg/kg - denosumab is dosed flat, so
#                                     no BW scaling is applied, unlike the
#                                     source's Dose_in_mgkg*BW formula)
#   MW_denosumab  = 147000  [g/mol]  denosumab-specific MW (FDA label);
#                                     source line 111 uses a generic
#                                     MW_EDG = 150000 for both species -
#                                     deliberately NOT reused here
#   V_plasma      = 3.126   [L]      already in R/01_parameters.R
#
# Conversion (same structural form as source line 11, adapted for a flat
# dose instead of Dose_in_mgkg*BW):
#   C0_plasma_uM = dose_mg / (V_plasma * 1000) / MW_denosumab * 1e6
#
# Central is amount-based (as in every file in this project), so the
# equivalent dosed AMOUNT is C0_plasma_uM * V_central - this also puts
# the dose on the same absolute uM/umol scale as FcRn_Conc_Bone, so the
# binding kinetics now see a physically meaningful ratio of drug to FcRn,
# not an arbitrary one.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.

library(rxode2)

source("R/01_parameters.R")
source("R/07_fcrn_parameters.R")

FR <- 0.715   # source line 140 [-]. Missing from 07_fcrn_parameters.R (as in 08/09).

# ---------------------------------------------------------------------------
# Real dose conversion
# ---------------------------------------------------------------------------

dose_mg      <- 60        # [mg] denosumab clinical dose, flat (FDA label)
MW_denosumab <- 147000    # [g/mol] denosumab-specific MW (FDA label)
V_plasma_L   <- unname(V_Plasma)   # [L], from 01_parameters.R (3.126 L, source line 300)

C0_plasma_uM <- dose_mg / (V_plasma_L * 1000) / MW_denosumab * 1e6

cat("=== Real dose conversion ===\n")
cat("C0_plasma_uM = dose_mg / (V_plasma * 1000) / MW_denosumab * 1e6\n")
cat("             =", dose_mg, "/ (", V_plasma_L, "* 1000 ) /", MW_denosumab, "* 1e6\n")
cat("             =", format(C0_plasma_uM, digits = 6), "uM\n\n")

# ---------------------------------------------------------------------------
# FcRn binding constants for E6a (drug, single site only) - unchanged from 09
# ---------------------------------------------------------------------------

k_off_6_EXG <- (KD6 / 1000) * k_on_6_EXG   # source line 201 formula; = 56.42 /h
stopifnot(abs(k_off_6_EXG - 56.42) < 1e-6)

FcRn_Conc_Bone <- FcRn_Availability * total_FcRn_in_nmole * 1e-3 / V_endosomal_Bone
FcRn_E6a_init_amount <- FcRn_Conc_Bone * V_E6a_Bone

# ---------------------------------------------------------------------------
# Organ groups (identical to 05/08/09)
# ---------------------------------------------------------------------------

generic_organs <- c("Lung", "Heart", "Muscle", "Skin", "Adipose", "Bone",
                     "Brain", "Kidney", "SI", "LI", "Pancreas", "Thymus",
                     "Spleen", "Other")
gut_organs   <- c("Spleen", "Pancreas", "SI", "LI")
plain_organs <- setdiff(generic_organs, gut_organs)
all_organs   <- c(generic_organs, "Liver")
non_bone_plain_organs <- setdiff(plain_organs, "Bone")

# ---------------------------------------------------------------------------
# Central / lymph node volumes, ARRAYSUM(LF[*]) over all 15 - unchanged
# ---------------------------------------------------------------------------

V_central   <- unname(V_Plasma)
V_ln        <- unname(V_LN)
L_LymphNode <- sum(LF[all_organs])

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

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
# Build the rxode2 model text - byte-for-byte identical to 09's model
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
# Dose: real IV-equivalent bolus. amt = C0_plasma_uM * V_central, so the
# initial Central concentration equals C0_plasma_uM exactly - putting the
# dose on the same absolute scale as FcRn_Conc_Bone. Simulate 0-20000 h.
# ---------------------------------------------------------------------------

dose_amt <- C0_plasma_uM * V_central

ev <- et(amt = dose_amt, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 10))

res <- rxSolve(mod, organ_params, ev,
                inits = c(BoneFcRn_E6a = unname(FcRn_E6a_init_amount)))

# ---------------------------------------------------------------------------
# Conservation check: total = total_live (free drug everywhere, incl.
# BoneBound_E6a) + BoneDegraded, should equal dose_amt throughout.
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
# Key result: Bone ABC at the real dose, vs 09's arbitrary-dose result, vs
# the no-FcRn distribution value, vs the published benchmark
# ---------------------------------------------------------------------------

last <- res[nrow(res), ]
C_central_last <- last$Central / V_central

bone_drug_real_dose <- last$BoneV + last$BoneIS + last$BoneVM + last$BoneE7 +
                       last$BoneE6a + last$BoneE7b + last$BoneISM + last$BoneBound_E6a

abc_real_dose_pct <- 100 * (bone_drug_real_dose / unname(V_Organ["Bone"])) / C_central_last

abc_arbitrary_dose_pct <- 6.36
no_fcrn_pct             <- 5.70
benchmark_pct           <- 7.27

cat("=== Key result: Bone ABC, real dose vs arbitrary dose vs references ===\n")
cat("Bone ABC at REAL dose (this file, 60 mg):     ", round(abc_real_dose_pct, 4), "%\n")
cat("Bone ABC at arbitrary '1 unit' dose (09):     ", abc_arbitrary_dose_pct, "%\n")
cat("No-FcRn distribution-only (06):               ", no_fcrn_pct, "%\n")
cat("Published benchmark (Shah & Betts 2013):      ", benchmark_pct, "%\n")
