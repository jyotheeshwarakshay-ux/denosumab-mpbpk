# ===========================================================================
# PROOF-OF-CONCEPT - NOT A CALIBRATED PREDICTION. Read this block first.
#
#   - Scope: FcRn rescue in E6a ONLY, single binding, DRUG ONLY. No
#     endogenous IgG (EDG) competition, no double (2nd-site) binding, no
#     binding anywhere else (VM/E7/E7b/ISM are trafficking-only, exactly
#     as in R/08_bone_fcrn_trafficking.R).
#   - FcRn pool is BONE-SCOPED: total_FcRn_in_nmole / V_endosomal_Bone,
#     NOT divided across all 15 tissues' endosomal volume the way the
#     source does (ARRAYSUM, source line 911). The other 14 tissues never
#     compete for this pool here. This OVERSTATES local FcRn availability
#     - the true rescue, if all 15 tissues competed for one shared pool,
#     is almost certainly SMALLER than what this file shows.
#   - Dose is a dimensionless "1 unit", NOT a real dose-scaled value.
#     FcRn binding is SATURABLE (unlike every prior file's linear
#     distribution equations), so the resulting ABC is dose-dependent
#     here - a different absolute dose would give a different rescued
#     ABC. 6.36% is illustrative of the mechanism's DIRECTION, not a
#     dose-calibrated bone ABC.
#   - Result: adding the rescue moves Bone ABC from 5.70% (no FcRn,
#     R/06) to 6.36% here - toward the 7.27% published benchmark -
#     confirming the FcRn rescue mechanism works in the right direction.
#     Nothing more than that should be read into the specific number.
# ===========================================================================
#
# 40-compartment rxode2 model: R/08_bone_fcrn_trafficking.R (the Bone
# trafficking belt, code-reviewed and verified) plus FcRn RESCUE in E6a
# ONLY - single binding, DRUG ONLY, no endogenous IgG competition, no
# double binding. A documented proof-of-concept scope, not a full FcRn
# reproduction.
#
# "BoneE6a" (state name, carried over unchanged from 08) is the free-drug
# pool in E6a - equivalent to what the task spec calls "BoneE6a_free".
# Kept the same name as 08 rather than renaming, since we're building on
# 08's compartment directly; "BoneBound_E6a" is the new FcRn-bound pool.
#
# DOCUMENTED SIMPLIFICATIONS (this file only, on top of 08's and earlier
# files' already-documented ones):
#
#   1. FcRn_Conc is scoped to BONE ONLY: FcRn_Availability *
#      total_FcRn_in_nmole * 1e-3 / V_endosomal_Bone. The source computes
#      this once, GLOBALLY, as total FcRn divided by
#      ARRAYSUM(V_endosomal[*]) across all 15 tissues (source line 911).
#      Using Bone's own endosomal volume as the sole denominator here
#      substantially inflates the local FcRn concentration relative to
#      what Bone would actually see if all 15 tissues competed for the
#      same total_FcRn_in_nmole pool - Bone's endosomal volume is a small
#      fraction of the body total. This is a deliberate proof-of-concept
#      simplification (matches the task's explicit instruction), not a
#      bug, but it means the rescue effect below is almost certainly
#      LARGER than a body-wide-consistent FcRn allocation would produce.
#
#   2. BoneFcRn_E6a has NO trafficking in/out (unlike the source, where
#      free FcRn moves between E7/E6a/E7b via CL_UP, source lines
#      779-781). The receptor pool is held local to E6a for this scope.
#
#   3. Bound drug leaving E6a arrives in E7b as FREE drug
#      (d/dt(BoneE7b) gains the same term BoneBound_E6a loses to
#      exocytosis). The source would keep it FcRn-bound
#      (Bound_E7b) through E7b's own binding kinetics; since this scope
#      has no binding in E7b, the complex is treated as releasing on
#      arrival.
#
#   4. ABSOLUTE DOSE SCALE (flagged, not resolved here): every prior file
#      in this project used a dimensionless "1 unit" bolus, which was
#      harmless for files 02-06 because those models are entirely linear
#      (ABC is a ratio that's scale-invariant). FcRn binding is NOT
#      linear - it depends on the absolute amount of drug relative to the
#      absolute amount of FcRn. FcRn_Conc here is derived from a real
#      physical quantity (total_FcRn_in_nmole = 1022 nmol, source line
#      119), while the "1 unit" dose is not on a real mg/kg -> uM scale
#      (source line 11 converts a real Dose_in_mgkg to real uM; this
#      project has never done that conversion). The rescue percentage
#      reported below should be read as a directional/mechanistic
#      demonstration (does adding FcRn binding measurably change Bone's
#      ABC), not as a dose-calibrated quantitative prediction.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.

library(rxode2)

source("R/01_parameters.R")
source("R/07_fcrn_parameters.R")

FR <- 0.715   # source line 140 [-]. Missing from 07_fcrn_parameters.R (as in 08).

# ---------------------------------------------------------------------------
# FcRn binding constants for E6a (drug, single site only)
# ---------------------------------------------------------------------------

k_off_6_EXG <- (KD6 / 1000) * k_on_6_EXG   # source line 201 formula; = 0.7 * 80.6 = 56.42 /h
stopifnot(abs(k_off_6_EXG - 56.42) < 1e-6)

# Bone-only-scoped FcRn concentration - see DOCUMENTED SIMPLIFICATION 1 above.
FcRn_Conc_Bone <- FcRn_Availability * total_FcRn_in_nmole * 1e-3 / V_endosomal_Bone

FcRn_E6a_init_amount <- FcRn_Conc_Bone * V_E6a_Bone   # INIT = FcRn_Conc * V_E6a_Bone, as an amount

# ---------------------------------------------------------------------------
# Organ groups (identical to 05/08)
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
# Build the rxode2 model text
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

# --- Central plasma (unchanged) -------------------------------------------

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

# --- Shared lymph node (unchanged) ----------------------------------------

ln_terms <- character(0)
for (o in all_organs) {
  ln_terms <- c(ln_terms, sprintf("+ LF_%s * (1 - sigmaIS_%s) * C_%sIS", o, o, o))
}
lines <- c(lines, sprintf("d/dt(LymphNode) <- %s - L_LymphNode * C_ln",
                           paste(ln_terms, collapse = " ")))

# --- Non-Bone plain tissues (9): unchanged --------------------------------

for (o in non_bone_plain_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            o, o, o, o, o, o, o, o),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            o, o, o, o, o, o, o))
}

# --- Gut tissues (4): unchanged --------------------------------------------

for (g in gut_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            g, g, g, g, g, g, g, g),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            g, g, g, g, g, g, g))
}

# --- Liver: unchanged --------------------------------------------------

portal_inflow_terms <- sapply(gut_organs, function(g)
  sprintf("+ (PLQ_%s - LF_%s) * C_%sV", g, g, g))
outlet_terms <- c("(PLQ_Liver - LF_Liver)",
                   sapply(gut_organs, function(g) sprintf("(PLQ_%s - LF_%s)", g, g)))
lines <- c(lines,
  sprintf("d/dt(LiverV) <- PLQ_Liver * C_central %s - C_LiverV * (%s) - LF_Liver * (1 - sigmaV_Liver) * C_LiverV",
          paste(portal_inflow_terms, collapse = " "), paste(outlet_terms, collapse = " + ")),
  "d/dt(LiverIS) <- LF_Liver * (1 - sigmaV_Liver) * C_LiverV - LF_Liver * (1 - sigmaIS_Liver) * C_LiverIS")

# --- BONE: vascular + interstitial (unchanged from 08) ---------------------

lines <- c(lines,
  "d/dt(BoneV)  <- PLQ_Bone * C_central - (PLQ_Bone - LF_Bone) * C_BoneV - LF_Bone * (1 - sigmaV_Bone) * C_BoneV - CL_UP_Bone * FR * C_BoneV + CL_UP_Bone * FR * C_BoneVM",
  "d/dt(BoneIS) <- LF_Bone * (1 - sigmaV_Bone) * C_BoneV - LF_Bone * (1 - sigmaIS_Bone) * C_BoneIS - CL_UP_Bone * (1 - FR) * C_BoneIS + CL_UP_Bone * (1 - FR) * C_BoneISM"
)

# --- BONE: trafficking belt, with FcRn rescue added to E6a/E7b ------------

lines <- c(lines,
  # VM: unchanged from 08
  "d/dt(BoneVM) <- CL_UP_Bone * FR * C_BoneV - CL_UP_Bone * FR * C_BoneVM - CL_UP_Bone * FR * C_BoneVM + CL_UP_Bone * FR * C_BoneE7b",

  # E7: unchanged from 08
  "d/dt(BoneE7) <- CL_UP_Bone * FR * C_BoneVM + CL_UP_Bone * (1 - FR) * C_BoneISM - CL_UP_Bone * C_BoneE7",

  # E6a (free): 08's belt terms UNCHANGED (full CL_UP outflow, feeding both
  # E7b's (1-Prob_deg) rescue-by-trafficking and Degraded's Prob_deg loss),
  # PLUS the new binding on/off terms
  "d/dt(BoneE6a) <- CL_UP_Bone * C_BoneE7 - CL_UP_Bone * C_BoneE6a - on_rate_E6a + off_rate_E6a",

  # BoneFcRn_E6a: consumed by binding, freed by unbinding, no trafficking (simplification 2)
  "d/dt(BoneFcRn_E6a) <- - on_rate_E6a + off_rate_E6a",

  # BoneBound_E6a: formed/released by binding, exocytosed forward to E7b
  "d/dt(BoneBound_E6a) <- on_rate_E6a - off_rate_E6a - CL_UP_Bone * C_BoneBound * V_E6a_Bone",

  # E7b: 08's terms UNCHANGED, PLUS rescued (bound-turned-free) drug arriving from E6a
  "d/dt(BoneE7b) <- CL_UP_Bone * (1 - Prob_deg) * C_BoneE6a - CL_UP_Bone * C_BoneE7b + CL_UP_Bone * C_BoneBound * V_E6a_Bone",

  # ISM: unchanged from 08
  "d/dt(BoneISM) <- CL_UP_Bone * (1 - FR) * C_BoneIS - CL_UP_Bone * (1 - FR) * C_BoneISM + CL_UP_Bone * (1 - FR) * C_BoneE7b - CL_UP_Bone * (1 - FR) * C_BoneISM",

  # Degraded: unchanged from 08 - the only sink, acts on free E6a only
  "d/dt(BoneDegraded) <- CL_UP_Bone * Prob_deg * C_BoneE6a"
)

model_txt <- paste(lines, collapse = "\n")

mod <- rxode2(model_txt)

stopifnot(length(mod$state) == 40)  # 38 (file 08) + BoneFcRn_E6a + BoneBound_E6a

# ---------------------------------------------------------------------------
# Single IV bolus of 1 unit into Central at t=0; BoneFcRn_E6a starts
# pre-loaded with the local free-FcRn pool. Simulate 0-20000 h.
# ---------------------------------------------------------------------------

ev <- et(amt = 1, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 10))

res <- rxSolve(mod, organ_params, ev,
                inits = c(BoneFcRn_E6a = unname(FcRn_E6a_init_amount)))

# ---------------------------------------------------------------------------
# Conservation check: total = total_live (free drug everywhere, incl.
# BoneBound_E6a) + BoneDegraded. Should be 1 throughout. BoneFcRn_E6a is
# NOT drug and is excluded.
# ---------------------------------------------------------------------------

live_cols <- c("Central", "LymphNode",
               as.vector(rbind(paste0(all_organs, "V"), paste0(all_organs, "IS"))),
               "BoneVM", "BoneE7", "BoneE6a", "BoneE7b", "BoneISM",
               "BoneBound_E6a")

total_live <- rowSums(res[, live_cols])
total_all  <- total_live + res$BoneDegraded

cat("=== Conservation check: total_live + BoneDegraded ===\n")
cat("min: ", min(total_all), "\n")
cat("max: ", max(total_all), "\n")
cat("Relative range:", (max(total_all) - min(total_all)) / mean(total_all), "\n")
cat("Total = 1 throughout (rel. range < 1e-8)?",
    (max(total_all) - min(total_all)) / mean(total_all) < 1e-8, "\n\n")

cat("total_live at t=0:     ", total_live[1], "\n")
cat("total_live at t=20000: ", tail(total_live, 1), "\n")
cat("BoneDegraded at t=20000:", tail(res$BoneDegraded, 1), "\n\n")

# ---------------------------------------------------------------------------
# Key result: Bone ABC with rescue, vs no-FcRn (5.70%) vs benchmark (7.27%)
# ---------------------------------------------------------------------------

last <- res[nrow(res), ]
C_central_last <- last$Central / V_central

bone_drug_with_rescue <- last$BoneV + last$BoneIS + last$BoneVM + last$BoneE7 +
                          last$BoneE6a + last$BoneE7b + last$BoneISM + last$BoneBound_E6a

abc_with_rescue_pct <- 100 * (bone_drug_with_rescue / unname(V_Organ["Bone"])) / C_central_last

no_fcrn_pct   <- 5.70
benchmark_pct <- 7.27

cat("=== Key result: Bone ABC with E6a rescue, at t =", last$time, "h ===\n")
cat("Bone ABC with rescue:            ", round(abc_with_rescue_pct, 4), "%\n")
cat("No-FcRn distribution-only (06):  ", no_fcrn_pct, "%\n")
cat("Published benchmark (2013):      ", benchmark_pct, "%\n")
cat("Rescue vs no-FcRn, absolute diff:", round(abc_with_rescue_pct - no_fcrn_pct, 4), "percentage points\n")
cat("Rescue vs no-FcRn, fold change:  ", round(abc_with_rescue_pct / no_fcrn_pct, 4), "\n")
cat("Rescue vs benchmark, ratio:      ", round(abc_with_rescue_pct / benchmark_pct, 4), "\n")
