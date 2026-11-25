# ===========================================================================
# PAPER 2 - STEP 1: whole-body linear elimination.  Read this block first.
#
#   ONE CHANGE vs R/05_full_15_tissues.R: linear elimination from CENTRAL
#   PLASMA, tracked in a new sink compartment "Eliminated".
#
#     d/dt(Central)    gains  - CL_lin * C_central
#     d/dt(Eliminated)  =     + CL_lin * C_central
#
#   CL_lin = 3.06 mL/h = 0.00306 L/h, the typical LINEAR clearance of
#   denosumab per 66 kg from Sutjandra et al. 2011, Clin Pharmacokinet
#   50:793-807 (doi 10.2165/11594240-000000000-00000).
#
#   ASSUMPTION (stated, not hidden): elimination happens from plasma only.
#   The Shah & Betts source eliminates antibody in the endosomes of every
#   tissue via FcRn; Fronton, Pilari & Huisinga 2014 (J PKPD 41:87-107)
#   also consider elimination from tissues. Plasma-only is the simplest
#   verifiable version; tissue elimination can be tested later as a
#   sensitivity analysis. No RANKL binding here - that is Step 2.
#
#   Dose: 60 mg as an IV-equivalent bolus into Central. Subcutaneous
#   absorption and bioavailability (64%, Sutjandra 2011) come in Step 2.
#
#   EXPECTED RESULTS (from an independent Python re-implementation of this
#   exact system, solved with a matrix exponential; the Python copy first
#   reproduced file 06's Bone ABC of 5.701% with CL = 0):
#     Check 1  mass balance (body + Eliminated = dose): rel. error < 1e-8
#     Check 2  terminal half-life ~ 56.6 days (1358 h)
#     Check 3  AUC_plasma * CL / dose = 1.000 (exact for plasma-only CL)
#     Plasma concentration: day 1 ~11.47, day 7 ~9.76, day 28 ~7.12,
#                           day 84 ~3.49, day 168 ~1.24 mg/L
#   If your R output differs from these by more than rounding, STOP and
#   report it - one of the two implementations has a bug.
#
#   WHAT THE HALF-LIFE MEANS: the Prolia label reports a mean half-life of
#   25.4 days after 60 mg SC. Our linear-only model gives ~57 days. That
#   gap is EXPECTED at this stage: the label value includes RANKL-mediated
#   (target-mediated) elimination, which this model does not have yet.
#   Step 2 tests whether adding RANKL binding closes the gap. Do not tune
#   CL_lin to force a match.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/05_full_15_tissues.R")   # builds `lines`, `organ_params`, `res`, `last` (CL = 0)

cat("\n\n=======================================================================\n")
cat("=== R/13_whole_body_elimination.R - Paper 2, Step 1 ===\n")
cat("=======================================================================\n\n")

# ---------------------------------------------------------------------------
# Steady-state volume of distribution from file 05's equilibrium (CL = 0):
# Vss = total amount in body / plasma concentration.
# ---------------------------------------------------------------------------

Vss <- sum(unlist(last[state_cols])) / (last$Central / V_central)
cat("Vss from file 05 equilibrium:", round(Vss, 3), "L   [Python: 5.865 L]\n\n")

# ---------------------------------------------------------------------------
# Add elimination: edit ONLY the Central equation, add the sink.
# ---------------------------------------------------------------------------

CL_lin <- 3.06 / 1000                      # L/h (Sutjandra 2011, per 66 kg)

i_central <- grep("^d/dt\\(Central\\)", lines)
stopifnot(length(i_central) == 1)

lines_el <- lines
lines_el[i_central] <- paste(lines[i_central], "- CL_lin * C_central")
lines_el <- c(lines_el, "d/dt(Eliminated) <- CL_lin * C_central")

mod_el <- rxode2(paste(lines_el, collapse = "\n"))
stopifnot(length(mod_el$state) == 33)      # 32 from file 05 + Eliminated

params_el <- c(organ_params, CL_lin = CL_lin)

dose <- 60                                  # mg
ev_el <- et(amt = dose, cmt = "Central", time = 0) |>
  et(seq(0, 20000, by = 1))

res_el <- rxSolve(mod_el, params_el, ev_el, atol = 1e-12, rtol = 1e-10)

# ---------------------------------------------------------------------------
# Check 1: mass balance - amount in body + amount eliminated = dose, always.
# ---------------------------------------------------------------------------

body_amt  <- rowSums(res_el[, state_cols])
total_amt <- body_amt + res_el$Eliminated
mass_rel_err <- max(abs(total_amt - dose)) / dose

cat("=== Check 1: mass balance (body + Eliminated = dose) ===\n")
cat("Max relative error:", signif(mass_rel_err, 3), "\n")
cat("Passed (< 1e-8)?   ", mass_rel_err < 1e-8, "\n\n")

# ---------------------------------------------------------------------------
# Check 2: terminal half-life from the log-linear slope of plasma
# concentration over the last part of the simulation (10000-20000 h).
# ---------------------------------------------------------------------------

C_pl <- res_el$Central / V_central
tail_idx <- res_el$time >= 10000
slope <- coef(lm(log(C_pl[tail_idx]) ~ res_el$time[tail_idx]))[2]
t_half_h <- log(2) / -unname(slope)

cat("=== Check 2: terminal half-life ===\n")
cat("Terminal half-life:", round(t_half_h), "h =", round(t_half_h / 24, 1),
    "days   [Python: 1358 h = 56.6 days]\n")
cat("Rough expectation ln2*Vss/CL:", round(log(2) * Vss / CL_lin / 24, 1), "days\n")
cat("Prolia label (60 mg SC, includes target-mediated elimination): 25.4 days\n\n")

# ---------------------------------------------------------------------------
# Check 3: AUC. For elimination from plasma only, AUC_plasma(0-inf) = dose/CL
# exactly. Trapezoid AUC to 20000 h plus a log-linear tail extrapolation.
# ---------------------------------------------------------------------------

tt <- res_el$time
auc_0_t <- sum(diff(tt) * (head(C_pl, -1) + tail(C_pl, -1)) / 2)
auc_inf <- auc_0_t + tail(C_pl, 1) / -unname(slope)
auc_ratio <- auc_inf * CL_lin / dose

cat("=== Check 3: AUC_plasma * CL / dose (must be ~1) ===\n")
cat("Ratio:", round(auc_ratio, 4), "   [Python: 1.000000]\n\n")

# ---------------------------------------------------------------------------
# Plasma concentrations at fixed days, for the comparison table.
# ---------------------------------------------------------------------------

cat("=== Plasma concentration (mg/L) ===\n")
for (d in c(1, 7, 28, 84, 168)) {
  cat(sprintf("Day %3d: %.3f   ", d, C_pl[tt == d * 24]))
}
cat("\n[Python: 11.474, 9.762, 7.118, 3.489, 1.244]\n")
