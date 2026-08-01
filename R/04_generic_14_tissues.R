# 30-compartment rxode2 prototype: Central plasma -> {14 generic tissues}
# vascular -> {14 generic tissues} interstitial -> ONE shared Lymph node ->
# back to Central. Distribution skeleton ONLY - no FcRn, no
# clearance/elimination anywhere.
#
# This extends R/03_two_tissue_prototype.R (Bone + Muscle) from 2 tissues to
# all 14 GENERIC tissues - all 15 Shah & Betts organs EXCEPT Liver. Liver is
# held back deliberately: its portal-vein inflow (it receives outflow from
# Spleen, Pancreas, SI and LI before returning to Central - see source
# lines 375-387) makes its vascular equation structurally different from
# every other organ. The 14 tissues here all share the identical per-organ
# pattern already verified twice:
#   Central <-> tissue_i vascular      : PLQ[i], both directions
#   tissue_i vascular -> interstitial  : LF[i]*(1-sigma_V[i])*C_vasc_i
#   tissue_i interstitial -> LymphNode : LF[i]*(1-sigma_IS[i])*C_int_i
#
# Per the task, this is built PROGRAMMATICALLY: the R code below loops over
# the 14 organ names and assembles the rxode2 model as a text string, rather
# than hand-writing 14 copies of the same block. Adding Liver later means
# adding its special-case vascular equation once, not restructuring this
# loop.
#
# Run from the project root (denosumab-mpbpk/) so the relative source()
# path below resolves.

library(rxode2)

source("R/01_parameters.R")

# ---------------------------------------------------------------------------
# The 14 generic tissues (all 15 Shah & Betts organs except Liver)
# ---------------------------------------------------------------------------

generic_organs <- c("Lung", "Heart", "Muscle", "Skin", "Adipose", "Bone",
                     "Brain", "Kidney", "SI", "LI", "Pancreas", "Thymus",
                     "Spleen", "Other")

stopifnot(all(generic_organs %in% organs))          # sanity vs 01_parameters.R
stopifnot(length(generic_organs) == 14)
stopifnot(!"Liver" %in% generic_organs)

# ---------------------------------------------------------------------------
# Central / lymph node volumes, and the scoped ARRAYSUM(LF[*]): sum of the
# 14 generic tissues' lymph flow into the shared node. Liver's LF is not
# included yet - it gets added when Liver joins the framework.
# ---------------------------------------------------------------------------

V_central   <- unname(V_Plasma)
V_ln        <- unname(V_LN)
L_LymphNode <- sum(LF[generic_organs])

# ---------------------------------------------------------------------------
# Per-organ parameters, named <param>_<Organ> so they can be referenced by
# name inside the generated model text.
# ---------------------------------------------------------------------------

organ_params <- c(V_central = V_central, V_ln = V_ln, L_LymphNode = L_LymphNode)

for (o in generic_organs) {
  organ_params[paste0("PLQ_", o)]      <- unname(PLQ[o])
  organ_params[paste0("LF_", o)]       <- unname(LF[o])
  organ_params[paste0("sigmaV_", o)]   <- unname(sigma_V[o])
  organ_params[paste0("sigmaIS_", o)]  <- unname(sigma_IS[o])
  organ_params[paste0("Vvasc_", o)]    <- unname(V_V[o])
  organ_params[paste0("Vinter_", o)]   <- unname(V_IntS[o])
}

# ---------------------------------------------------------------------------
# Build the rxode2 model text by looping over the 14 organs.
# Compartments (30 total): Central, LymphNode, and per organ <Organ>V /
# <Organ>IS (vascular / interstitial amounts).
# ---------------------------------------------------------------------------

lines <- c(
  "C_central <- Central / V_central",
  "C_ln <- LymphNode / V_ln"
)

for (o in generic_organs) {
  lines <- c(lines,
    sprintf("C_%sV  <- %sV  / Vvasc_%s",  o, o, o),
    sprintf("C_%sIS <- %sIS / Vinter_%s", o, o, o)
  )
}

central_terms <- character(0)
ln_terms      <- character(0)
for (o in generic_organs) {
  central_terms <- c(central_terms,
    sprintf("- PLQ_%s * C_central + PLQ_%s * C_%sV", o, o, o))
  ln_terms <- c(ln_terms,
    sprintf("+ LF_%s * (1 - sigmaIS_%s) * C_%sIS", o, o, o))
}

lines <- c(lines,
  sprintf("d/dt(Central) <- %s + L_LymphNode * C_ln",
          paste(central_terms, collapse = " ")),
  sprintf("d/dt(LymphNode) <- %s - L_LymphNode * C_ln",
          paste(ln_terms, collapse = " "))
)

for (o in generic_organs) {
  lines <- c(lines,
    sprintf("d/dt(%sV)  <- PLQ_%s * C_central - PLQ_%s * C_%sV - LF_%s * (1 - sigmaV_%s) * C_%sV",
            o, o, o, o, o, o, o),
    sprintf("d/dt(%sIS) <- LF_%s * (1 - sigmaV_%s) * C_%sV - LF_%s * (1 - sigmaIS_%s) * C_%sIS",
            o, o, o, o, o, o, o)
  )
}

model_txt <- paste(lines, collapse = "\n")

mod <- rxode2(model_txt)

stopifnot(length(mod$state) == 30)   # Central + LymphNode + 14*(V + IS)

# ---------------------------------------------------------------------------
# Single IV bolus of 1 unit into Central at t=0; simulate 0-20000 h.
# A long window is used deliberately: the two-tissue prototype
# (03_two_tissue_prototype.R) showed Bone's interstitial time constant is
# ~456 h, so 2000 h (~4.4 time constants) was not enough to fully
# equilibrate the lymph node balance, even though mass conservation held
# throughout. 20000 h (~44 time constants) was confirmed there to close the
# balance to machine precision.
# ---------------------------------------------------------------------------

ev <- et(amt = 1, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 10))

res <- rxSolve(mod, organ_params, ev)

# ---------------------------------------------------------------------------
# Check 1: mass conservation across all 30 compartments, at every timepoint.
# ---------------------------------------------------------------------------

state_cols <- c("Central", "LymphNode",
                as.vector(rbind(paste0(generic_organs, "V"),
                                 paste0(generic_organs, "IS"))))

total_mass <- rowSums(res[, state_cols])

mass_min <- min(total_mass)
mass_max <- max(total_mass)
mass_range <- mass_max - mass_min
mass_rel_range <- mass_range / mean(total_mass)
mass_conserved <- mass_rel_range < 1e-8

cat("=== Check 1: mass conservation (30 compartments) ===\n")
cat("Total mass min:      ", mass_min, "\n")
cat("Total mass max:      ", mass_max, "\n")
cat("Absolute range:      ", mass_range, "\n")
cat("Relative range:      ", mass_rel_range, "\n")
cat("Mass conserved (rel. range < 1e-8)?", mass_conserved, "\n\n")

# ---------------------------------------------------------------------------
# Check 2: lymph node inflow/outflow balance at the final timepoint.
# ---------------------------------------------------------------------------

last <- res[nrow(res), ]

C_central_last <- last$Central / V_central
C_ln_last      <- last$LymphNode / V_ln

ln_inflow <- 0
for (o in generic_organs) {
  C_is_last <- last[[paste0(o, "IS")]] / organ_params[paste0("Vinter_", o)]
  ln_inflow <- ln_inflow +
    organ_params[paste0("LF_", o)] * (1 - organ_params[paste0("sigmaIS_", o)]) * C_is_last
}
ln_outflow <- L_LymphNode * C_ln_last

ln_diff <- unname(ln_inflow - ln_outflow)
ln_rel_diff <- abs(ln_diff) / unname(ln_outflow)
ln_balanced <- ln_rel_diff < 1e-6

cat("=== Check 2: lymph node balance at t =", last$time, "h ===\n")
cat("Inflow  (sum of 14 tissues' filtered interstitial outflow):", unname(ln_inflow), "\n")
cat("Outflow (L_LymphNode * C_ln):                                ", unname(ln_outflow), "\n")
cat("Difference:                                                  ", ln_diff, "\n")
cat("Relative difference:                                         ", ln_rel_diff, "\n")
cat("Lymph node balanced (rel. diff < 1e-6)?", ln_balanced, "\n\n")

# ---------------------------------------------------------------------------
# Check 3: sanity readout - interstitial-to-plasma concentration ratio per
# tissue at the final timepoint (a preview of the biodistribution
# coefficient). All should be well below 1. At full equilibrium, with
# sigma_IS uniform (0.2) across tissues and PLQ >> LF so C_vasc ~= C_central,
# the ratio reduces to ~(1-sigma_V)/(1-sigma_IS): tissues with HIGH sigma_V
# (tighter vasculature, e.g. Brain at 0.99) should sit LOWEST, and tissues
# with LOW sigma_V (leakier vasculature, e.g. Bone and Spleen at 0.85 - the
# lowest values in the table) should sit HIGHEST.
# ---------------------------------------------------------------------------

ratios <- sapply(generic_organs, function(o) {
  C_is_last <- last[[paste0(o, "IS")]] / organ_params[paste0("Vinter_", o)]
  unname(C_is_last / C_central_last)
})

ratios_sorted <- sort(ratios)

cat("=== Check 3: C_interstitial / C_central per tissue, at t =", last$time, "h (sorted) ===\n")
print(round(ratios_sorted, 4))
