# 6-compartment rxode2 prototype: Central plasma -> {Bone, Muscle} vascular
# -> {Bone, Muscle} interstitial -> ONE shared Lymph node -> back to
# Central. Distribution skeleton ONLY - no FcRn, no clearance/elimination
# anywhere. Extends R/02_single_tissue_prototype.R from one tissue to two,
# to verify that a shared lymph node correctly pools outflow from multiple
# tissues before the mechanism is scaled to all 15.
#
# Per-tissue structure (identical for Bone and Muscle, verified already in
# 02_single_tissue_prototype.R):
#   Central <-> tissue vascular       : PLQ, both directions
#   tissue vascular -> interstitial   : LF*(1-sigma_V)*C_vasc
#   tissue interstitial -> LymphNode  : LF*(1-sigma_IS)*C_int
#
# New in this file: the lymph node is shared. It receives the filtered
# interstitial outflow from BOTH tissues, and its outflow rate is now the
# real ARRAYSUM(LF[*]) - scoped to the two tissues present here:
#   L_LymphNode = LF[Bone] + LF[Muscle]
# (in 02_single_tissue_prototype.R this was a single-tissue placeholder,
# L_LymphNode = LF_bone, because Bone was the only contributor).
#
# Run from the project root (denosumab-mpbpk/) so the relative source()
# path below resolves.

library(rxode2)

source("R/01_parameters.R")

# ---------------------------------------------------------------------------
# Pull out Bone and Muscle parameters (plus central plasma and lymph node
# volumes)
# ---------------------------------------------------------------------------

V_central <- unname(V_Plasma)
V_ln      <- unname(V_LN)

V_bv  <- unname(V_V["Bone"])
V_bis <- unname(V_IntS["Bone"])
V_mv  <- unname(V_V["Muscle"])
V_mis <- unname(V_IntS["Muscle"])

PLQ_bone    <- unname(PLQ["Bone"])
PLQ_muscle  <- unname(PLQ["Muscle"])
LF_bone     <- unname(LF["Bone"])
LF_muscle   <- unname(LF["Muscle"])

sigma_V_bone    <- unname(sigma_V["Bone"])
sigma_V_muscle  <- unname(sigma_V["Muscle"])
sigma_IS_bone   <- unname(sigma_IS["Bone"])
sigma_IS_muscle <- unname(sigma_IS["Muscle"])

# Scoped ARRAYSUM(LF[*]): sum of both tissues' lymph flow into the shared node
L_LymphNode <- LF_bone + LF_muscle

params <- c(
  V_central       = V_central,
  V_ln            = V_ln,
  V_bv            = V_bv,
  V_bis           = V_bis,
  V_mv            = V_mv,
  V_mis           = V_mis,
  PLQ_bone        = PLQ_bone,
  PLQ_muscle      = PLQ_muscle,
  LF_bone         = LF_bone,
  LF_muscle       = LF_muscle,
  sigma_V_bone    = sigma_V_bone,
  sigma_V_muscle  = sigma_V_muscle,
  sigma_IS_bone   = sigma_IS_bone,
  sigma_IS_muscle = sigma_IS_muscle,
  L_LymphNode     = L_LymphNode
)

# ---------------------------------------------------------------------------
# Model: six amount-based compartments (Central, BoneV, BoneIS, MuscleV,
# MuscleIS, LymphNode). Every transfer term that leaves one compartment
# appears with the identical term (same rate, opposite sign) entering
# another - a closed loop with no sink, so total mass is conserved by
# construction.
# ---------------------------------------------------------------------------

mod <- rxode2({
  C_central <- Central   / V_central
  C_bv      <- BoneV     / V_bv
  C_bis     <- BoneIS    / V_bis
  C_mv      <- MuscleV   / V_mv
  C_mis     <- MuscleIS  / V_mis
  C_ln      <- LymphNode / V_ln

  d/dt(Central) <- -PLQ_bone   * C_central -
                     PLQ_muscle * C_central +
                     PLQ_bone   * C_bv +
                     PLQ_muscle * C_mv +
                     L_LymphNode * C_ln

  d/dt(BoneV)   <-  PLQ_bone * C_central -
                     PLQ_bone * C_bv -
                     LF_bone * (1 - sigma_V_bone) * C_bv

  d/dt(BoneIS)  <-  LF_bone * (1 - sigma_V_bone) * C_bv -
                     LF_bone * (1 - sigma_IS_bone) * C_bis

  d/dt(MuscleV) <-  PLQ_muscle * C_central -
                     PLQ_muscle * C_mv -
                     LF_muscle * (1 - sigma_V_muscle) * C_mv

  d/dt(MuscleIS) <- LF_muscle * (1 - sigma_V_muscle) * C_mv -
                      LF_muscle * (1 - sigma_IS_muscle) * C_mis

  # Lymph node receives the filtered interstitial outflow from BOTH
  # tissues, and drains back to Central at the pooled rate L_LymphNode
  d/dt(LymphNode) <- LF_bone   * (1 - sigma_IS_bone)   * C_bis +
                      LF_muscle * (1 - sigma_IS_muscle) * C_mis -
                      L_LymphNode * C_ln
})

# ---------------------------------------------------------------------------
# Single IV bolus of 1 unit into Central at t=0; simulate 0-2000 h
# ---------------------------------------------------------------------------

ev <- et(amt = 1, cmt = "Central", time = 0) |>
  et(seq(0, 2000, by = 1))

res <- rxSolve(mod, params, ev)

# ---------------------------------------------------------------------------
# Check 1: mass conservation - total antibody amount across all 6
# compartments at every timepoint.
# ---------------------------------------------------------------------------

total_mass <- res$Central + res$BoneV + res$BoneIS +
              res$MuscleV + res$MuscleIS + res$LymphNode

mass_min <- min(total_mass)
mass_max <- max(total_mass)
mass_range <- mass_max - mass_min
mass_rel_range <- mass_range / mean(total_mass)
mass_conserved <- mass_rel_range < 1e-8

cat("=== Check 1: mass conservation ===\n")
cat("Total mass min:      ", mass_min, "\n")
cat("Total mass max:      ", mass_max, "\n")
cat("Absolute range:      ", mass_range, "\n")
cat("Relative range:      ", mass_rel_range, "\n")
cat("Mass conserved (rel. range < 1e-8)?", mass_conserved, "\n\n")

# ---------------------------------------------------------------------------
# Check 2: lymph node inflow/outflow balance at steady state (last
# timepoint, t=2000 h). Recompute the two flux terms directly from the
# solved concentrations - if the node is neither accumulating nor leaking,
# these should match.
# ---------------------------------------------------------------------------

last <- res[nrow(res), ]

C_bis_last <- last$BoneIS   / V_bis
C_mis_last <- last$MuscleIS / V_mis
C_ln_last  <- last$LymphNode / V_ln

ln_inflow  <- LF_bone   * (1 - sigma_IS_bone)   * C_bis_last +
              LF_muscle * (1 - sigma_IS_muscle) * C_mis_last
ln_outflow <- L_LymphNode * C_ln_last

ln_diff <- ln_inflow - ln_outflow
ln_rel_diff <- abs(ln_diff) / ln_outflow
ln_balanced <- ln_rel_diff < 1e-6

cat("=== Check 2: lymph node balance at t =", last$time, "h ===\n")
cat("Inflow  (Bone + Muscle filtered interstitial outflow):", ln_inflow, "\n")
cat("Outflow (L_LymphNode * C_ln):                          ", ln_outflow, "\n")
cat("Difference:                                            ", ln_diff, "\n")
cat("Relative difference:                                   ", ln_rel_diff, "\n")
cat("Lymph node balanced (rel. diff < 1e-6)?", ln_balanced, "\n")
