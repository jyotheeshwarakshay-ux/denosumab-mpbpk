# 38-compartment rxode2 model: the full 15-tissue distribution skeleton
# (R/05_full_15_tissues.R) plus the FcRn TRAFFICKING belt (pinocytosis /
# endosomal transit only - no FcRn binding, no endogenous IgG) added to
# BONE ONLY. The other 14 tissues are byte-for-byte the same equations as
# 05. This is the first non-conservative piece added to this project:
# BoneDegraded is a genuine sink (see source lines 454-467 and the
# Prob_deg mechanism documented earlier), not a simplification.
#
# New Bone compartments (free EXG only, trafficking terms only):
#   BoneVM, BoneE7, BoneE6a, BoneE7b, BoneISM   (source lines 424-480,
#     CL_up/CL_UP terms only - all K_on/K_off FcRn-binding terms and all
#     k_on_PS polyspecificity terms omitted, to be added later)
#   BoneDegraded                                 (accumulator, no outflow)
#
# Bone's existing BoneV and BoneIS equations (from 05) gain the matching
# pinocytosis exchange terms with BoneVM and BoneISM respectively - every
# source vascular/interstitial equation already has these terms (source
# lines 361-364, 472-475), 05 just had nowhere for them to go before now.
#
# NOTE: FR = 0.715 (source line 140) is used throughout the trafficking
# equations but was not transcribed into R/07_fcrn_parameters.R - that
# file's parameter list (as specified) omitted it. Defined locally below
# instead of silently editing 07.
#
# All states remain AMOUNT-based (matching every prior file in this
# project), NOT concentration-based like the source's Berkeley Madonna
# equations. Each term borrowed from the source is the pre-division
# (i.e. already amount-rate) content of that equation's bracket - no
# trailing "/V[...]" is applied here, unlike the source, because our
# state variables are amounts, not concentrations.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.

library(rxode2)

source("R/01_parameters.R")
source("R/07_fcrn_parameters.R")

FR <- 0.715   # source line 140 [-]. Missing from 07_fcrn_parameters.R; defined here instead.

# ---------------------------------------------------------------------------
# Organ groups (identical to 05)
# ---------------------------------------------------------------------------

generic_organs <- c("Lung", "Heart", "Muscle", "Skin", "Adipose", "Bone",
                     "Brain", "Kidney", "SI", "LI", "Pancreas", "Thymus",
                     "Spleen", "Other")
gut_organs   <- c("Spleen", "Pancreas", "SI", "LI")
plain_organs <- setdiff(generic_organs, gut_organs)
all_organs   <- c(generic_organs, "Liver")

non_bone_plain_organs <- setdiff(plain_organs, "Bone")   # 9 tissues, fully unchanged

stopifnot(length(all_organs) == 15)
stopifnot(length(non_bone_plain_organs) == 9)

# ---------------------------------------------------------------------------
# Central / lymph node volumes, ARRAYSUM(LF[*]) over all 15 - unchanged
# ---------------------------------------------------------------------------

V_central   <- unname(V_Plasma)
V_ln        <- unname(V_LN)
L_LymphNode <- sum(LF[all_organs])

# ---------------------------------------------------------------------------
# Per-organ parameters (all 15 tissues, including Liver) - unchanged from 05
# ---------------------------------------------------------------------------

organ_params <- c(V_central = V_central, V_ln = V_ln, L_LymphNode = L_LymphNode,
                   FR = FR, Prob_deg = Prob_deg,
                   CL_UP_Bone = unname(CL_UP_Bone),
                   V_VM_Bone  = unname(V_VM_Bone),
                   V_E7_Bone  = unname(V_E7_Bone),
                   V_E6a_Bone = unname(V_E6a_Bone),
                   V_E7b_Bone = unname(V_E7b_Bone),
                   V_ISM_Bone = unname(V_ISM_Bone))

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
  "C_BoneVM  <- BoneVM  / V_VM_Bone",
  "C_BoneE7  <- BoneE7  / V_E7_Bone",
  "C_BoneE6a <- BoneE6a / V_E6a_Bone",
  "C_BoneE7b <- BoneE7b / V_E7b_Bone",
  "C_BoneISM <- BoneISM / V_ISM_Bone")

# --- Central plasma (unchanged from 05 - Bone's return credit is still via
#     C_BoneV, untouched by the new membrane compartments) -----------------

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

# --- Shared lymph node (unchanged from 05) --------------------------------

ln_terms <- character(0)
for (o in all_organs) {
  ln_terms <- c(ln_terms, sprintf("+ LF_%s * (1 - sigmaIS_%s) * C_%sIS", o, o, o))
}
lines <- c(lines, sprintf("d/dt(LymphNode) <- %s - L_LymphNode * C_ln",
                           paste(ln_terms, collapse = " ")))

# --- Non-Bone plain tissues (9): byte-for-byte identical to 05 -----------

for (o in non_bone_plain_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            o, o, o, o, o, o, o, o),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            o, o, o, o, o, o, o))
}

# --- Gut tissues (4): unchanged from 05 -----------------------------------

for (g in gut_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            g, g, g, g, g, g, g, g),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            g, g, g, g, g, g, g))
}

# --- Liver: unchanged from 05 ----------------------------------------------

portal_inflow_terms <- sapply(gut_organs, function(g)
  sprintf("+ (PLQ_%s - LF_%s) * C_%sV", g, g, g))
outlet_terms <- c("(PLQ_Liver - LF_Liver)",
                   sapply(gut_organs, function(g) sprintf("(PLQ_%s - LF_%s)", g, g)))
lines <- c(lines,
  sprintf("d/dt(LiverV) <- PLQ_Liver * C_central %s - C_LiverV * (%s) - LF_Liver * (1 - sigmaV_Liver) * C_LiverV",
          paste(portal_inflow_terms, collapse = " "), paste(outlet_terms, collapse = " + ")),
  "d/dt(LiverIS) <- LF_Liver * (1 - sigmaV_Liver) * C_LiverV - LF_Liver * (1 - sigmaIS_Liver) * C_LiverIS")

# --- BONE: vascular + interstitial, now with membrane exchange terms -----

lines <- c(lines,
  "d/dt(BoneV)  <- PLQ_Bone * C_central - (PLQ_Bone - LF_Bone) * C_BoneV - LF_Bone * (1 - sigmaV_Bone) * C_BoneV - CL_UP_Bone * FR * C_BoneV + CL_UP_Bone * FR * C_BoneVM",
  "d/dt(BoneIS) <- LF_Bone * (1 - sigmaV_Bone) * C_BoneV - LF_Bone * (1 - sigmaIS_Bone) * C_BoneIS - CL_UP_Bone * (1 - FR) * C_BoneIS + CL_UP_Bone * (1 - FR) * C_BoneISM"
)

# --- BONE: new FcRn trafficking belt (free drug only, no binding) --------

lines <- c(lines,
  # VM: in from V, out to V, out to E7, in from E7b
  "d/dt(BoneVM) <- CL_UP_Bone * FR * C_BoneV - CL_UP_Bone * FR * C_BoneVM - CL_UP_Bone * FR * C_BoneVM + CL_UP_Bone * FR * C_BoneE7b",

  # E7: in from VM, in from ISM, out to E6a
  "d/dt(BoneE7) <- CL_UP_Bone * FR * C_BoneVM + CL_UP_Bone * (1 - FR) * C_BoneISM - CL_UP_Bone * C_BoneE7",

  # E6a: in from E7, out (full CL_UP) - splits downstream into E7b and Degraded
  "d/dt(BoneE6a) <- CL_UP_Bone * C_BoneE7 - CL_UP_Bone * C_BoneE6a",

  # E7b: in from E6a ((1-Prob_deg) fraction only), out to VM+ISM (full)
  "d/dt(BoneE7b) <- CL_UP_Bone * (1 - Prob_deg) * C_BoneE6a - CL_UP_Bone * C_BoneE7b",

  # ISM: in from IS, out to E7, in from E7b, out to IS
  "d/dt(BoneISM) <- CL_UP_Bone * (1 - FR) * C_BoneIS - CL_UP_Bone * (1 - FR) * C_BoneISM + CL_UP_Bone * (1 - FR) * C_BoneE7b - CL_UP_Bone * (1 - FR) * C_BoneISM",

  # Degraded: accumulator, no outflow - the ONLY sink in the whole model
  "d/dt(BoneDegraded) <- CL_UP_Bone * Prob_deg * C_BoneE6a"
)

model_txt <- paste(lines, collapse = "\n")

mod <- rxode2(model_txt)

stopifnot(length(mod$state) == 38)  # 2 (Central, LymphNode) + 15*2 (organs) + 5 (Bone membranes) + 1 (Degraded)

# ---------------------------------------------------------------------------
# Single IV bolus of 1 unit into Central at t=0; simulate 0-20000 h
# ---------------------------------------------------------------------------

ev <- et(amt = 1, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 10))

res <- rxSolve(mod, organ_params, ev)

# ---------------------------------------------------------------------------
# Conservation check: total_live (everything except BoneDegraded) should
# DECREASE over time (degradation is now active); total_live + BoneDegraded
# should stay exactly 1 (the dose) throughout, to floating-point precision.
# ---------------------------------------------------------------------------

live_cols <- c("Central", "LymphNode",
               as.vector(rbind(paste0(all_organs, "V"), paste0(all_organs, "IS"))),
               "BoneVM", "BoneE7", "BoneE6a", "BoneE7b", "BoneISM")

total_live <- rowSums(res[, live_cols])
total_all  <- total_live + res$BoneDegraded

cat("=== Check: total_live (excludes BoneDegraded) ===\n")
cat("total_live at t=0:     ", total_live[1], "\n")
cat("total_live at t=20000: ", tail(total_live, 1), "\n")
cat("Decreased?", tail(total_live, 1) < total_live[1], "\n\n")

cat("=== Check: total_live + BoneDegraded (should be 1 throughout) ===\n")
cat("min: ", min(total_all), "\n")
cat("max: ", max(total_all), "\n")
cat("Relative range:", (max(total_all) - min(total_all)) / mean(total_all), "\n")
cat("Conserved to floating-point precision (rel. range < 1e-8)?",
    (max(total_all) - min(total_all)) / mean(total_all) < 1e-8, "\n\n")

cat("BoneDegraded at t=20000:", tail(res$BoneDegraded, 1), "\n")
