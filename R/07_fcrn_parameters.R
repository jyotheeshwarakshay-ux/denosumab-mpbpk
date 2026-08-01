# FcRn / endothelial parameters for BONE only, transcribed from
# sources/PSP4-8-738-s002.txt. Sourced ON TOP of R/01_parameters.R (which
# is NOT modified) - this file adds the parameters needed for FcRn
# trafficking that 01_parameters.R does not cover (that file only has
# tissue volumes/flows/reflection coefficients for the distribution
# skeleton).
#
# Parameters and derived volumes ONLY - no ODEs, no model equations.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.

source("R/01_parameters.R")

# ---------------------------------------------------------------------------
# Transcribed parameters (source line cited on each)
# ---------------------------------------------------------------------------

Total_Endothelial_Cell <- 1.422e9          # source line 125 [number]. Mouse total endothelial cell count.
human_endo_scale_factor <- 603.7           # source line 130 [-]. Scales mouse endothelial count to human.
Scale_Factor <- human_endo_scale_factor    # source line 132 [-]. "To simulate human."

Endothelial_Cell_Frac <- c(Bone = 0.1478)  # source line 223 (index 7 = Bone) [-]

CL_Up_in_nL_per_hour_per_million_cells <- 150  # source line 134 [nL/h per 1e6 cells]

pino_time <- 10.8 / 60                     # source line 136 [h]. (Defined as `Pino_time`
                                            # at line 136, referenced lowercase `pino_time`
                                            # in the V_endosomal formula at line 910 - same
                                            # variable, Berkeley Madonna is case-insensitive.)

total_FcRn_in_nmole <- 1022                # source line 119 (`total_FcRn_in_nmole_par`,
                                            # "1022 for human"), assigned to
                                            # `total_FcRn_in_nmole` at line 214 [nmol]

FcRn_Availability <- 1.0                   # source line 160 [-]

E6a_Vol_pct <- 0.33                        # source line 183 [-]
E7b_Vol_pct <- (1 - E6a_Vol_pct) / 2       # source line 184 [-]
E7_Vol_pct  <- (1 - E6a_Vol_pct) / 2       # source line 185 [-]

tau_VM  <- 1 / 60                          # source line 180 [h]
tau_ISM <- 1 / 60                          # source line 181 [h]

Prob_deg <- 0.95                           # source line 162 [-]
FcRn_recycle_fraction <- 0.99              # source line 138 [-]

# FcRn binding constants (1st-site, exogenous drug)
KD6 <- 700                                 # source line 50 (`KD6_WT`) [nM], pH 6.0
KD7 <- 154077                              # source line 51 (`KD7_WT`) [nM], pH 7.4
k_on_6_EXG <- 8.06e7 / 1e6                 # source line 169 [1/(uM*h)]
k_on_7_EXG <- 1.61e7 / 5 / 1e6             # source line 170 [1/(uM*h)]

on_rate_ratio_1st_vs_2nd <- 83.7           # source line 127 (`on_rate_ratio_1st_vs_2nd_binding`) [-]
kdeg_FcRn_Ab <- log(2) / 11.1              # source line 164 (`LOGN(2)/11.1`) [1/h]

# ---------------------------------------------------------------------------
# Derived quantities for Bone - formulas match source lines 908-923 exactly.
# Printed with the formula alongside each value so the arithmetic can be
# checked line by line.
# ---------------------------------------------------------------------------

# source line 908: Organ_Endothelial_Cell[i] = Total_Endothelial_Cell * Endothelial_Cell_Frac[i] * Scale_Factor
Organ_Endothelial_Cell_Bone <- Total_Endothelial_Cell * Endothelial_Cell_Frac["Bone"] * Scale_Factor

# source line 921: CL_UP[i] = CL_Up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell[i] * 1E-15
# NOTE - unit label question (report as instructed, not resolved here):
# Table S3 (PDF) labels CLup as [ml/hr]. The arithmetic below, combined
# with how CL_UP is used downstream (divided into V_v[j] terms which are
# in L, alongside PLQ/LF also in L/h), is only dimensionally consistent if
# CL_UP actually comes out in L/h, not mL/h as labeled. This is flagged,
# not resolved - verify against your own unit-consistency check once the
# FcRn ODEs are wired in, the same way the PLQ[Other] table/code
# discrepancy was flagged rather than silently picked one way.
CL_UP_Bone <- CL_Up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell_Bone * 1e-15

# source line 910: V_endosomal[i] = CL_Up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell[i] * 1E-6 * 1E-9
V_endosomal_Bone <- CL_Up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell_Bone * 1e-6 * 1e-9

# source lines 917-919: V_E7/V_E6a/V_E7b[i] = V_endosomal[i] * <vol_pct>
V_E7_Bone  <- V_endosomal_Bone * E7_Vol_pct
V_E6a_Bone <- V_endosomal_Bone * E6a_Vol_pct
V_E7b_Bone <- V_endosomal_Bone * E7b_Vol_pct

# source line 923-924: V_VM[i] = CL_UP[i]*tau_VM ; V_ISM[i] = CL_UP[i]*tau_ISM
V_VM_Bone  <- CL_UP_Bone * tau_VM
V_ISM_Bone <- CL_UP_Bone * tau_ISM

# ---------------------------------------------------------------------------
# Print every computed value with the formula used
# ---------------------------------------------------------------------------

cat("Organ_Endothelial_Cell[Bone] = Total_Endothelial_Cell * Endothelial_Cell_Frac[Bone] * Scale_Factor\n")
cat("                             =", Total_Endothelial_Cell, "*", unname(Endothelial_Cell_Frac["Bone"]), "*", Scale_Factor, "\n")
cat("                             =", format(unname(Organ_Endothelial_Cell_Bone), scientific = TRUE, digits = 6), "\n\n")

cat("CL_UP[Bone] = CL_Up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell[Bone] * 1e-15\n")
cat("            =", CL_Up_in_nL_per_hour_per_million_cells, "*", format(unname(Organ_Endothelial_Cell_Bone), scientific = TRUE, digits = 6), "* 1e-15\n")
cat("            =", format(unname(CL_UP_Bone), digits = 6), "  [unit-label question flagged above - see comment]\n\n")

cat("V_endosomal[Bone] = CL_Up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell[Bone] * 1e-6 * 1e-9\n")
cat("                  =", CL_Up_in_nL_per_hour_per_million_cells, "*", format(pino_time, digits = 6), "*",
    format(unname(Organ_Endothelial_Cell_Bone), scientific = TRUE, digits = 6), "* 1e-15\n")
cat("                  =", format(unname(V_endosomal_Bone), digits = 6), "\n\n")

cat("V_E7[Bone]  = V_endosomal[Bone] * E7_Vol_pct  =", format(unname(V_endosomal_Bone), digits = 6), "*", E7_Vol_pct,
    "=", format(unname(V_E7_Bone), digits = 6), "\n")
cat("V_E6a[Bone] = V_endosomal[Bone] * E6a_Vol_pct =", format(unname(V_endosomal_Bone), digits = 6), "*", E6a_Vol_pct,
    "=", format(unname(V_E6a_Bone), digits = 6), "\n")
cat("V_E7b[Bone] = V_endosomal[Bone] * E7b_Vol_pct =", format(unname(V_endosomal_Bone), digits = 6), "*", E7b_Vol_pct,
    "=", format(unname(V_E7b_Bone), digits = 6), "\n\n")

cat("V_VM[Bone]  = CL_UP[Bone] * tau_VM  =", format(unname(CL_UP_Bone), digits = 6), "*", tau_VM,
    "=", format(unname(V_VM_Bone), digits = 6), "\n")
cat("V_ISM[Bone] = CL_UP[Bone] * tau_ISM =", format(unname(CL_UP_Bone), digits = 6), "*", tau_ISM,
    "=", format(unname(V_ISM_Bone), digits = 6), "\n")
