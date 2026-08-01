# Physiological parameters for the Shah & Betts (2012) platform mAb PBPK model.
# Source: Jones HM et al. 2019, CPT: Pharmacometrics & Systems Pharmacology,
# supplementary file PSP4-8-738-s002.txt (Berkeley Madonna code), lines 64-335.
# Values transcribed here are the HUMAN parameter set (as used for human
# simulations in the source file, i.e. Table S2 equivalent).
#
# No ODEs / model logic here - parameters only.

# ---------------------------------------------------------------------------
# Organ name-to-index key (source lines 64-78)
# Kept for traceability back to the source file's array indexing; all
# vectors below are named by organ so this lookup isn't required downstream.
# ---------------------------------------------------------------------------

organ_index <- c(
  Lung     = 1,
  Liver    = 2,
  Heart    = 3,
  Muscle   = 4,
  Skin     = 5,
  Adipose  = 6,
  Bone     = 7,
  Brain    = 8,
  Kidney   = 9,
  SI       = 10,
  LI       = 11,
  Pancreas = 12,
  Thymus   = 13,
  Spleen   = 14,
  Other    = 15
)

organs <- names(organ_index)

# ---------------------------------------------------------------------------
# Reflection coefficients (source lines 142-158) [dimensionless]
# ---------------------------------------------------------------------------

sigma_V <- c(
  Lung     = 0.95,
  Liver    = 0.85,
  Heart    = 0.95,
  Muscle   = 0.95,
  Skin     = 0.95,
  Adipose  = 0.95,
  Bone     = 0.85,
  Brain    = 0.99,
  Kidney   = 0.90,
  SI       = 0.90,
  LI       = 0.95,
  Pancreas = 0.90,
  Thymus   = 0.90,
  Spleen   = 0.85,
  Other    = 0.95
)

# sigma_IS is a single value (0.2) applied to all 15 organs in the source.
sigma_IS <- setNames(rep(0.2, 15), organs)

# ---------------------------------------------------------------------------
# Tissue sub-compartment volumes (source lines 233-297) [L]
# Original source values are in mL (divided by 1000 in-line); kept as
# literal mL/1000 expressions here so each number can be checked directly
# against Table S2's mL figures.
# ---------------------------------------------------------------------------

# Vascular volume per organ
V_V <- c(
  Lung     = 55/1000,
  Liver    = 183/1000,
  Heart    = 13.1/1000,
  Muscle   = 662/1000,
  Skin     = 127/1000,
  Adipose  = 148/1000,
  Bone     = 224/1000,
  Brain    = 31.9/1000,
  Kidney   = 18.2/1000,
  SI       = 6.15/1000,
  LI       = 8.74/1000,
  Pancreas = 5.7/1000,
  Thymus   = 0.353/1000,
  Spleen   = 26.8/1000,
  Other    = 204/1000
)

# Interstitial volume per organ
V_IntS <- c(
  Lung     = 300/1000,
  Liver    = 429/1000,
  Heart    = 48.8/1000,
  Muscle   = 3910/1000,
  Skin     = 1125/1000,
  Adipose  = 2289/1000,
  Bone     = 1891/1000,
  Brain    = 261/1000,
  Kidney   = 49.8/1000,
  SI       = 67.1/1000,
  LI       = 95.3/1000,
  Pancreas = 18/1000,
  Thymus   = 1.09/1000,
  Spleen   = 44.3/1000,
  Other    = 831/1000
)

# Cellular volume per organ
V_C <- c(
  Lung     = 595/1000,
  Liver    = 1371/1000,
  Heart    = 267/1000,
  Muscle   = 24815/1000,
  Skin     = 2035/1000,
  Adipose  = 10840/1000,
  Bone     = 7817/1000,
  Brain    = 1124/1000,
  Kidney   = 247/1000,
  SI       = 305/1000,
  LI       = 434/1000,
  Pancreas = 74.7/1000,
  Thymus   = 4.65/1000,
  Spleen   = 127/1000,
  Other    = 3626/1000
)

# Total organ volume
V_Organ <- c(
  Lung     = 1000/1000,
  Liver    = 2143/1000,
  Heart    = 341/1000,
  Muscle   = 30078/1000,
  Skin     = 3408/1000,
  Adipose  = 13465/1000,
  Bone     = 10165/1000,
  Brain    = 1450/1000,
  Kidney   = 332/1000,
  SI       = 385/1000,
  LI       = 548/1000,
  Pancreas = 104/1000,
  Thymus   = 6.41/1000,
  Spleen   = 221/1000,
  Other    = 4852/1000
)

# ---------------------------------------------------------------------------
# Plasma and lymph node volumes (source lines 300-301) [L]
# Given directly in L in the source (human values).
# ---------------------------------------------------------------------------

V_Plasma <- 3.126
V_LN     <- 0.274

# ---------------------------------------------------------------------------
# Organ plasma flow rates (source lines 305-319) [L/h]
# Original source values are in mL/h (divided by 1000 in-line).
# ---------------------------------------------------------------------------

PLQ <- c(
  Lung     = 181913/1000,
  Liver    = 13210/1000,
  Heart    = 7752/1000,
  Muscle   = 33469/1000,
  Skin     = 11626/1000,
  Adipose  = 11233/1000,
  Bone     = 2591/1000,
  Brain    = 21453/1000,
  Kidney   = 36402/1000,
  SI       = 12368/1000,
  LI       = 12867/1000,
  Pancreas = 3056/1000,
  Thymus   = 353/1000,
  Spleen   = 6343/1000,
  # NOTE: source line 319 hardcodes this as a literal 9190/1000, in the same
  # style as every other PLQ[i] entry - it is NOT written as a formula in
  # the source. However, 9190 = PLQ[Lung] (181913, i.e. cardiac output)
  # minus the sum of all other 13 organ flows (172723), so it is exactly
  # the residual that makes peripheral flow balance to cardiac output.
  # Table S2 (Jones et al. 2019) lists Other plasma flow as 5521 mL/h,
  # which does NOT match this 9190 figure - the discrepancy is real and
  # unresolved; 9190 is what the model code actually runs with, so that is
  # what is transcribed here.
  Other    = 9190/1000
)

# ---------------------------------------------------------------------------
# Lymph flow rates (source lines 321-335) [L/h]
# LF = PLQ * 0.002 for every organ EXCEPT Liver, which sums the portal
# venous inputs (Liver + SI + LI + Spleen + Pancreas, net of their own
# lymph flow) before applying the 0.002 factor - reproduced exactly as in
# the source rather than simplified, since it is a direct physiological
# input value, not model logic.
# ---------------------------------------------------------------------------

LF <- setNames(PLQ * 0.002, organs)

LF["Liver"] <- (PLQ["Liver"] +
                  (PLQ["SI"]       - LF["SI"]) +
                  (PLQ["LI"]       - LF["LI"]) +
                  (PLQ["Spleen"]   - LF["Spleen"]) +
                  (PLQ["Pancreas"] - LF["Pancreas"])) * 0.002
