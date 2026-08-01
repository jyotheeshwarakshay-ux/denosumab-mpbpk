# Minimal 4-compartment rxode2 prototype: Central plasma -> Bone vascular
# -> Bone interstitial -> Lymph node -> back to Central. Distribution
# skeleton ONLY - no FcRn, no clearance/elimination anywhere. Purpose:
# verify the convective transport mechanism conserves total antibody mass
# before this mechanism is scaled up to all 15 tissues.
#
# v2: adds the Lymph Node as its own compartment, matching the real source
# topology (PSP4-8-738-s002.txt lines 487-493, 833-897) instead of the v1
# shortcut that returned interstitial lymph flow directly to Central. In
# the source, ALL 15 tissues drain into ONE shared lymph node compartment,
# whose total outflow rate is L_LymphNode = ARRAYSUM(LF[*]) (source line
# 208) - the sum of every tissue's lymph flow, not just one tissue's. With
# only Bone in this prototype, Bone is the pool's only contributor, so the
# scoped equivalent used here is L_LymphNode = LF_bone. That keeps the
# system closed for this test, but it is a placeholder: once the other 14
# tissues are added, L_LymphNode must become the sum of all their LF
# values, not just Bone's.
#
# Run from the project root (denosumab-mpbpk/) so the relative source()
# path below resolves.

library(rxode2)

source("R/01_parameters.R")

# ---------------------------------------------------------------------------
# Pull out only the bone-specific parameters (plus central plasma volume
# and lymph node volume)
# ---------------------------------------------------------------------------

V_central     <- unname(V_Plasma)
V_bv          <- unname(V_V["Bone"])
V_bis         <- unname(V_IntS["Bone"])
V_ln          <- unname(V_LN)
PLQ_bone      <- unname(PLQ["Bone"])
LF_bone       <- unname(LF["Bone"])
sigma_V_bone  <- unname(sigma_V["Bone"])
sigma_IS_bone <- unname(sigma_IS["Bone"])

# Placeholder for ARRAYSUM(LF[*]) from the source - see note above. Equals
# LF_bone here because Bone is the only tissue in this prototype.
L_LymphNode <- LF_bone

params <- c(
  V_central     = V_central,
  V_bv          = V_bv,
  V_bis         = V_bis,
  V_ln          = V_ln,
  PLQ_bone      = PLQ_bone,
  LF_bone       = LF_bone,
  sigma_V_bone  = sigma_V_bone,
  sigma_IS_bone = sigma_IS_bone,
  L_LymphNode   = L_LymphNode
)

# ---------------------------------------------------------------------------
# Model: four amount-based compartments (Central, BoneV, BoneIS, LymphNode).
# Every transfer term that leaves one compartment appears with the identical
# term (same rate, opposite sign) entering another - a closed loop with no
# sink, so total mass is conserved by construction. That is the thing being
# tested below.
# ---------------------------------------------------------------------------

mod <- rxode2({
  C_central <- Central   / V_central
  C_bv      <- BoneV     / V_bv
  C_bis     <- BoneIS    / V_bis
  C_ln      <- LymphNode / V_ln

  # Central <-> Bone vascular: plasma flow delivers, plasma flow returns
  # Bone vascular -> Bone interstitial: lymph-driven convection, filtered
  #   by the vascular reflection coefficient (sigma_V)
  # Bone interstitial -> Lymph node: lymph-driven convection, filtered by
  #   the interstitial reflection coefficient (sigma_IS) - matches source
  #   line 492 exactly, unlike v1 which omitted this filter
  # Lymph node -> Central: total lymph flow carries pooled lymph back to
  #   plasma (source lines 872, 897)

  d/dt(Central) <- -PLQ_bone * C_central +
                     PLQ_bone * C_bv +
                     L_LymphNode * C_ln

  d/dt(BoneV)   <-  PLQ_bone * C_central -
                     PLQ_bone * C_bv -
                     LF_bone * (1 - sigma_V_bone) * C_bv

  d/dt(BoneIS)  <-  LF_bone * (1 - sigma_V_bone) * C_bv -
                     LF_bone * (1 - sigma_IS_bone) * C_bis

  d/dt(LymphNode) <- LF_bone * (1 - sigma_IS_bone) * C_bis -
                       L_LymphNode * C_ln
})

# ---------------------------------------------------------------------------
# Single IV bolus of 1 unit into Central at t=0; simulate 0-2000 h
# ---------------------------------------------------------------------------

ev <- et(amt = 1, cmt = "Central", time = 0) |>
  et(seq(0, 2000, by = 1))

res <- rxSolve(mod, params, ev)

# ---------------------------------------------------------------------------
# Mass conservation check: total antibody amount (Central + BoneV + BoneIS +
# LymphNode) at every timepoint. All four states are amount-based (not
# concentrations), so summing them directly gives total mass.
# ---------------------------------------------------------------------------

total_mass <- res$Central + res$BoneV + res$BoneIS + res$LymphNode

mass_min <- min(total_mass)
mass_max <- max(total_mass)
mass_range <- mass_max - mass_min
mass_rel_range <- mass_range / mean(total_mass)

conserved <- mass_rel_range < 1e-8

cat("Total mass min:      ", mass_min, "\n")
cat("Total mass max:      ", mass_max, "\n")
cat("Absolute range:      ", mass_range, "\n")
cat("Relative range:      ", mass_rel_range, "\n")
cat("Mass conserved (rel. range < 1e-8)?", conserved, "\n")
