# 32-compartment rxode2 prototype: all 15 Shah & Betts tissues, extending
# R/04_generic_14_tissues.R by adding LIVER with correct portal-vein
# routing. This is the ONE new structural feature; the 10 "plain" generic
# tissues (Lung, Heart, Muscle, Skin, Adipose, Bone, Brain, Kidney, Thymus,
# Other) are untouched from file 04. Distribution skeleton ONLY - no FcRn,
# no clearance/elimination anywhere.
#
# Portal routing (verified verbatim against source lines 375-387, 855-872):
#
#   LIVER vascular receives:
#     - arterial inflow from lung: PLQ[Liver] * C_central     (like every
#       tissue)
#     - portal inflow from Spleen, Pancreas, SI, LI: for each gut organ g,
#       (PLQ[g]-LF[g]) * C_vasc_g
#
#   LIVER vascular drains at the COMBINED flow, at LIVER's own
#   concentration:
#     outlet = C_LiverV * ( (PLQ[Liver]-LF[Liver])
#                          + (PLQ[Spleen]-LF[Spleen])
#                          + (PLQ[Pancreas]-LF[Pancreas])
#                          + (PLQ[SI]-LF[SI])
#                          + (PLQ[LI]-LF[LI]) )
#   plus the normal vascular->interstitial lymph term
#   (1-sigma_V[Liver])*LF[Liver]*C_LiverV, and LiverIS drains to the shared
#   node exactly like every other tissue.
#
#   The four gut organs (Spleen, Pancreas, SI, LI):
#     - vascular->interstitial and interstitial->node terms UNCHANGED
#       (same generic pattern as file 04)
#     - arterial inflow from central UNCHANGED: + PLQ[g] * C_central
#     - but the venous return leg is redirected: gutV loses
#       (PLQ[g]-LF[g])*C_gutV to LIVER's portal inflow (not back to
#       Central at its own concentration). Central's credit for this flow
#       is folded into the liver-return terms below, at C_LiverV.
#
#   Central plasma: for the 10 plain tissues, unchanged
#   (-PLQ*C_central + PLQ*C_vasc). For Liver and the 4 gut organs, Central
#   loses the usual arterial term (-PLQ*C_central) but its ONLY return
#   credit from this group is the combined liver-outlet flow, all at
#   C_LiverV:
#     + (PLQ[Liver]-LF[Liver])  * C_LiverV
#     + (PLQ[Spleen]-LF[Spleen])* C_LiverV
#     + (PLQ[Pancreas]-LF[Pancreas]) * C_LiverV
#     + (PLQ[SI]-LF[SI])        * C_LiverV
#     + (PLQ[LI]-LF[LI])        * C_LiverV
#   which exactly matches Liver's combined outlet term - see the term-by-
#   term conservation argument in the project notes.
#
#   L_LymphNode is now ARRAYSUM(LF[*]) over all 15 tissues (14 generic +
#   Liver), not just the 14 from file 04.
#
# Built programmatically (loops over organ name vectors), matching file 04.
#
# Run from the project root (denosumab-mpbpk/) so the relative source()
# path below resolves.

library(rxode2)

source("R/01_parameters.R")

# ---------------------------------------------------------------------------
# Organ groups
# ---------------------------------------------------------------------------

generic_organs <- c("Lung", "Heart", "Muscle", "Skin", "Adipose", "Bone",
                     "Brain", "Kidney", "SI", "LI", "Pancreas", "Thymus",
                     "Spleen", "Other")                    # same 14 as file 04

gut_organs   <- c("Spleen", "Pancreas", "SI", "LI")         # portal-draining subset
plain_organs <- setdiff(generic_organs, gut_organs)         # 10, fully unchanged from file 04

all_organs <- c(generic_organs, "Liver")                    # all 15

stopifnot(length(generic_organs) == 14)
stopifnot(length(gut_organs) == 4)
stopifnot(length(plain_organs) == 10)
stopifnot(length(all_organs) == 15)
stopifnot(all(all_organs %in% organs))                      # sanity vs 01_parameters.R

# ---------------------------------------------------------------------------
# Central / lymph node volumes, and ARRAYSUM(LF[*]) over all 15 tissues
# ---------------------------------------------------------------------------

V_central   <- unname(V_Plasma)
V_ln        <- unname(V_LN)
L_LymphNode <- sum(LF[all_organs])

# ---------------------------------------------------------------------------
# Per-organ parameters (all 15 tissues, including Liver)
# ---------------------------------------------------------------------------

organ_params <- c(V_central = V_central, V_ln = V_ln, L_LymphNode = L_LymphNode)

for (o in all_organs) {
  organ_params[paste0("PLQ_", o)]      <- unname(PLQ[o])
  organ_params[paste0("LF_", o)]       <- unname(LF[o])
  organ_params[paste0("sigmaV_", o)]   <- unname(sigma_V[o])
  organ_params[paste0("sigmaIS_", o)]  <- unname(sigma_IS[o])
  organ_params[paste0("Vvasc_", o)]    <- unname(V_V[o])
  organ_params[paste0("Vinter_", o)]   <- unname(V_IntS[o])
}

# ---------------------------------------------------------------------------
# Build the rxode2 model text by looping over organ groups.
# Compartments (32 total): Central, LymphNode, and per organ (all 15)
# <Organ>V / <Organ>IS (vascular / interstitial amounts).
# ---------------------------------------------------------------------------

lines <- c(
  "C_central <- Central / V_central",
  "C_ln <- LymphNode / V_ln"
)

for (o in all_organs) {
  lines <- c(lines,
    sprintf("C_%sV  <- %sV  / Vvasc_%s",  o, o, o),
    sprintf("C_%sIS <- %sIS / Vinter_%s", o, o, o)
  )
}

# --- Central plasma -----------------------------------------------------

central_terms <- character(0)

for (o in plain_organs) {
  central_terms <- c(central_terms,
    sprintf("- PLQ_%s * C_central + PLQ_%s * C_%sV", o, o, o))
}

for (g in gut_organs) {
  central_terms <- c(central_terms,
    sprintf("- PLQ_%s * C_central + (PLQ_%s - LF_%s) * C_LiverV", g, g, g))
}

central_terms <- c(central_terms,
  "- PLQ_Liver * C_central + (PLQ_Liver - LF_Liver) * C_LiverV")

lines <- c(lines,
  sprintf("d/dt(Central) <- %s + L_LymphNode * C_ln",
          paste(central_terms, collapse = " ")))

# --- Shared lymph node ----------------------------------------------------

ln_terms <- character(0)
for (o in all_organs) {
  ln_terms <- c(ln_terms,
    sprintf("+ LF_%s * (1 - sigmaIS_%s) * C_%sIS", o, o, o))
}

lines <- c(lines,
  sprintf("d/dt(LymphNode) <- %s - L_LymphNode * C_ln",
          paste(ln_terms, collapse = " ")))

# --- Plain tissues (10): unchanged generic pattern from file 04 -----------

for (o in plain_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - PLQ_%s * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            o, o, o, o, o, o, o),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            o, o, o, o, o, o, o)
  )
}

# --- Gut tissues (4): interstitial side unchanged; vascular return
#     redirected to the portal path instead of back to Central -----------

for (g in gut_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - (PLQ_%s - LF_%s) * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            g, g, g, g, g, g, g, g),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            g, g, g, g, g, g, g)
  )
}

# --- Liver: arterial + portal inflow, combined outlet ---------------------

portal_inflow_terms <- sapply(gut_organs, function(g)
  sprintf("+ (PLQ_%s - LF_%s) * C_%sV", g, g, g))

outlet_terms <- c("(PLQ_Liver - LF_Liver)",
                   sapply(gut_organs, function(g) sprintf("(PLQ_%s - LF_%s)", g, g)))

liver_V_eq <- sprintf(
  "d/dt(LiverV) <- PLQ_Liver * C_central %s - C_LiverV * (%s) - LF_Liver * (1 - sigmaV_Liver) * C_LiverV",
  paste(portal_inflow_terms, collapse = " "),
  paste(outlet_terms, collapse = " + ")
)

liver_IS_eq <- "d/dt(LiverIS) <- LF_Liver * (1 - sigmaV_Liver) * C_LiverV - LF_Liver * (1 - sigmaIS_Liver) * C_LiverIS"

lines <- c(lines, liver_V_eq, liver_IS_eq)

model_txt <- paste(lines, collapse = "\n")

mod <- rxode2(model_txt)

stopifnot(length(mod$state) == 32)   # Central + LymphNode + 15*(V + IS)

# ---------------------------------------------------------------------------
# Single IV bolus of 1 unit into Central at t=0; simulate 0-20000 h (same
# long window as file 04, for the same reason: Bone's slow interstitial
# time constant).
# ---------------------------------------------------------------------------

ev <- et(amt = 1, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 10))

res <- rxSolve(mod, organ_params, ev)

# ---------------------------------------------------------------------------
# Check 1: mass conservation across all 32 compartments, at every timepoint.
# ---------------------------------------------------------------------------

state_cols <- c("Central", "LymphNode",
                as.vector(rbind(paste0(all_organs, "V"),
                                 paste0(all_organs, "IS"))))

total_mass <- rowSums(res[, state_cols])

mass_min <- min(total_mass)
mass_max <- max(total_mass)
mass_range <- mass_max - mass_min
mass_rel_range <- mass_range / mean(total_mass)
mass_conserved <- mass_rel_range < 1e-8

cat("=== Check 1: mass conservation (32 compartments) ===\n")
cat("Total mass min:      ", mass_min, "\n")
cat("Total mass max:      ", mass_max, "\n")
cat("Absolute range:      ", mass_range, "\n")
cat("Relative range:      ", mass_rel_range, "\n")
cat("Mass conserved (rel. range < 1e-8)?", mass_conserved, "\n\n")

last <- res[nrow(res), ]

C_central_last <- last$Central / V_central
C_ln_last      <- last$LymphNode / V_ln

conc_last <- function(o, side) {
  vol_name <- if (side == "V") paste0("Vvasc_", o) else paste0("Vinter_", o)
  last[[paste0(o, side)]] / organ_params[vol_name]
}

# ---------------------------------------------------------------------------
# Check 2: lymph node inflow/outflow balance at the final timepoint
# (now summed over all 15 tissues).
# ---------------------------------------------------------------------------

ln_inflow <- 0
for (o in all_organs) {
  C_is_last <- conc_last(o, "IS")
  ln_inflow <- ln_inflow +
    organ_params[paste0("LF_", o)] * (1 - organ_params[paste0("sigmaIS_", o)]) * C_is_last
}
ln_outflow <- L_LymphNode * C_ln_last

ln_diff <- unname(ln_inflow - ln_outflow)
ln_rel_diff <- abs(ln_diff) / unname(ln_outflow)
ln_balanced <- ln_rel_diff < 1e-6

cat("=== Check 2: lymph node balance at t =", last$time, "h ===\n")
cat("Inflow  (sum of 15 tissues' filtered interstitial outflow):", unname(ln_inflow), "\n")
cat("Outflow (L_LymphNode * C_ln):                                ", unname(ln_outflow), "\n")
cat("Difference:                                                  ", ln_diff, "\n")
cat("Relative difference:                                         ", ln_rel_diff, "\n")
cat("Lymph node balanced (rel. diff < 1e-6)?", ln_balanced, "\n\n")

# ---------------------------------------------------------------------------
# Check 3 (NEW): portal consistency - total flow INTO liver's vascular
# (arterial + 4 portal) vs total flow OUT (combined outlet + liver's own
# lymph term), at the final timepoint. If liver's vascular compartment is
# neither accumulating nor leaking blood, these should match.
# ---------------------------------------------------------------------------

C_LiverV_last <- conc_last("Liver", "V")

liver_inflow <- organ_params["PLQ_Liver"] * C_central_last
for (g in gut_organs) {
  C_gV_last <- conc_last(g, "V")
  liver_inflow <- liver_inflow +
    (organ_params[paste0("PLQ_", g)] - organ_params[paste0("LF_", g)]) * C_gV_last
}

combined_outlet_rate <- organ_params["PLQ_Liver"] - organ_params["LF_Liver"]
for (g in gut_organs) {
  combined_outlet_rate <- combined_outlet_rate +
    (organ_params[paste0("PLQ_", g)] - organ_params[paste0("LF_", g)])
}

liver_outflow <- C_LiverV_last * combined_outlet_rate +
  organ_params["LF_Liver"] * (1 - organ_params["sigmaV_Liver"]) * C_LiverV_last

liver_diff <- unname(liver_inflow - liver_outflow)
liver_rel_diff <- abs(liver_diff) / unname(liver_outflow)
liver_balanced <- liver_rel_diff < 1e-6

cat("=== Check 3: liver vascular portal consistency at t =", last$time, "h ===\n")
cat("Inflow  (arterial + 4 portal):                 ", unname(liver_inflow), "\n")
cat("Outflow (combined outlet + liver lymph term):  ", unname(liver_outflow), "\n")
cat("Difference:                                    ", liver_diff, "\n")
cat("Relative difference:                           ", liver_rel_diff, "\n")
cat("Liver vascular balanced (rel. diff < 1e-6)?", liver_balanced, "\n\n")

# ---------------------------------------------------------------------------
# Check 4: interstitial:plasma ratio per tissue (all 15), sorted, at the
# final timepoint. The 14 generic tissues should be unchanged from file 04;
# Liver's ratio is new.
# ---------------------------------------------------------------------------

ratios <- sapply(all_organs, function(o) {
  C_is_last <- conc_last(o, "IS")
  unname(C_is_last / C_central_last)
})

ratios_sorted <- sort(ratios)

cat("=== Check 4: C_interstitial / C_central per tissue, at t =", last$time, "h (sorted) ===\n")
print(round(ratios_sorted, 4))

cat("\n=== Comparison to file 04 (14 generic tissues only) ===\n")
file04_ratios <- c(Brain = 0.0125, Muscle = 0.0625, Other = 0.0625, Skin = 0.0625,
                    Adipose = 0.0625, Heart = 0.0625, LI = 0.0625, Lung = 0.0625,
                    Pancreas = 0.1250, SI = 0.1250, Thymus = 0.1250, Kidney = 0.1250,
                    Spleen = 0.1874, Bone = 0.1874)
comparison <- data.frame(
  tissue      = names(file04_ratios),
  file04      = unname(file04_ratios),
  file05      = round(unname(ratios[names(file04_ratios)]), 4)
)
comparison$diff <- comparison$file05 - comparison$file04
print(comparison)

cat("\nLiver interstitial:plasma ratio (new):", round(unname(ratios["Liver"]), 4), "\n")
