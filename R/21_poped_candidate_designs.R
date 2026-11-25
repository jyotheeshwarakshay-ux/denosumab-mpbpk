# ===========================================================================
# PAPER 2 - STEP 5b: why late samples inform Kss, and practical designs
# that keep every sample above the assay LLOQ.  Read this block first.
#
#   PART 1 (mechanism test): typical total and FREE denosumab after a 60 mg
#   SC dose (Choi 2025 model, as R/20), at 2160, 2880, 3600, 4380 h, and
#   when free drug falls below Kss (1.56 nM) and total drug below the assay
#   LLOQ. LLOQ = 20 ng/mL (Choi et al. 2025, Section 2.2) = 0.1361 nM
#   with MW 147 kDa.
#   PART 2 (designs): PopED FIM (FO, reduced FIM; same model, IIV, error,
#   dosing and N = 80 as R/20) for candidate sampling schedules. Every
#   sample is at most 4200 h after a dose, i.e. before the TYPICAL
#   concentration crosses the LLOQ (~4318-4372 h). Individuals with faster
#   elimination will still go below LLOQ earlier; the FIM ignores this
#   censoring, so Step 6 (simulation-estimation with BLQ handling) must
#   confirm.
#
#   EXPECTED RESULTS (reference = CRAN PopED 0.7.0, deSolve 1.42, run on
#   26 Nov 2026; [brackets] = dev PopED 0.7.0.9002, deSolve 1.40, R 4.3.3):
#   A/B agree within 1 point across PopED versions; C 1.5, D 3.0; E–G differ widely (ill-conditioned FIM).
#   Part 1 (after dose 1; identical in both PopED versions):
#     +2160 h total 8.5498, free 8.1014 nM, free/Kss 5.19
#     +2880 h total 3.9967, free 3.5578, 2.28
#     +3600 h total 1.4099, free 1.0200, 0.65
#     +4380 h total 0.0947, free 0.0288, 0.02 (total below LLOQ)
#     free < Kss from +3404 h (dose 1), +3456.64 h (dose 2), +3458.28 h (dose 3)
#     total < LLOQ from +4318 h, +4370.64 h, +4372.28 h
#   Part 2 (RSE%):             KSS            KINT          R0            KSYN           CL           VC            d_KSS
#     A v2 (Paper 1 final)  24  78.9 [79.8]    8.3 [8.4]     14.1 [14.1]   6.7 [6.8]      3.3 [3.3]    9.2 [9.2]     722.7 [713.0]
#     B v2 + 3600           27  38.8 [38.7]    7.3 [7.2]     14.0 [13.9]   3.9 [3.8]      3.0 [3.0]    7.0 [7.0]     433.9 [438.6]
#     C v2 + 3000/3600/4200 33  29.9 [28.4]    7.1 [6.8]     13.8 [13.6]   3.8 [3.7]      3.0 [3.0]    6.8 [6.8]     268.0 [260.8]
#     D 8 samples, 3 late   24  33.6 [30.6]    14.0 [12.9]   14.7 [14.1]   4.2 [4.2]      3.0 [3.0]    7.0 [7.0]     298.0 [285.3]
#     E 4 samples incl late 12  384.6 [396.9]  79.7 [105.7]  68.4 [122.3]  106.2 [172.1]  18.0 [14.3]  83.1 [124.9]  1937.3 [1837.3]
#     F 4 samples, no late  12  362.0 [433.4]  35.2 [35.9]   19.1 [20.1]   21.3 [24.3]    7.9 [9.4]    31.8 [38.0]   3911.8 [4281.0]
#     G 6 late-weighted     18  126.8 [157.0]  56.8 [57.8]   35.6 [53.8]   45.7 [68.8]    3.3 [3.3]    42.3 [51.8]   1541.7 [1336.5]
#
# Run from the project root (denosumab-mpbpk/).
# ===========================================================================

suppressMessages({ library(PopED); library(deSolve) })
cat("PopED", as.character(packageVersion("PopED")), "| deSolve", as.character(packageVersion("deSolve")), "\n\n")

src <- readLines("R/20_poped_baseline.R")
eval(parse(text = src[grep("^MW <-", src):(grep("^for \\(nm in", src) - 1)]))   # model, ff, sfg, make_db

# --- Part 1 ----------------------------------------------------------------
Kss_t <- 1.56; LLOQ <- 20 / 147000 * 1e3
p <- c(KA = 0.0078, VC = 1.58, VP = 6.06, CL = 0.006, Q = 0.20, KINT = 0.022, Kss = Kss_t, KSYN = 0.01, KDEG = 0.01 / 15.23)
times <- sort(unique(c(seq(0, 18 * month_h, by = 1), dose_times)))
evd <- data.frame(var = "depot", time = dose_times, value = dose_nmol, method = "add")
out <- as.data.frame(suppressWarnings(ode(c(depot = 0, Ctot = 0, Cp = 0, Rtot = 15.23), times, ode_qss, p,
                     events = list(data = evd), method = "lsoda", rtol = 1e-8, atol = 1e-10)))
bq <- out$Ctot - out$Rtot - Kss_t; out$Cfree <- 0.5 * (bq + sqrt(bq^2 + 4 * Kss_t * out$Ctot))
cat(sprintf("=== Part 1: typical total and free drug (LLOQ = %.4f nM) ===\n", LLOQ))
for (o in c(2160, 2880, 3600, 4380)) {
  k <- which(out$time == o)
  cat(sprintf("  +%4d h: total %.4f, free %.4f nM | free/Kss %.2f | total >= LLOQ: %s\n",
              o, out$Ctot[k], out$Cfree[k], out$Cfree[k] / Kss_t, out$Ctot[k] >= LLOQ))
}
for (d in seq_along(dose_times)) {
  w <- out$time > dose_times[d] + 300 & out$time < (if (d < 3) dose_times[d + 1] else Inf)
  i1 <- which(w & out$Cfree < Kss_t)[1]; i2 <- which(w & out$Ctot < LLOQ)[1]
  cat(sprintf("  dose %d: free < Kss from +%.2f h; total < LLOQ from +%.2f h\n",
              d, out$time[i1] - dose_times[d], out$time[i2] - dose_times[d]))
}

# --- Part 2 ----------------------------------------------------------------
designs <- list(
  "A v2 (Paper 1 final)"       = c(12, 24, 72, 168, 336, 720, 1440, 2160),
  "B v2 + 3600"                = c(12, 24, 72, 168, 336, 720, 1440, 2160, 3600),
  "C v2 + 3000/3600/4200"      = c(12, 24, 72, 168, 336, 720, 1440, 2160, 3000, 3600, 4200),
  "D 8 samples, 3 moved late"  = c(12, 24, 168, 720, 2160, 3000, 3600, 4200),
  "E 4 samples incl. late"     = c(168, 720, 3600, 4200),
  "F 4 samples, none late"     = c(24, 168, 720, 2160),
  "G 6 samples, late-weighted" = c(168, 2160, 3000, 3400, 3800, 4200))
cat("\n=== Part 2: expected RSE% per design (80 subjects, 3 doses) ===\n")
for (nm in names(designs)) {
  ev <- evaluate_design(make_db(designs[[nm]]))
  cat(sprintf("%-28s n=%2d | KSS %5.1f | KINT %5.1f | R0 %5.1f | KSYN %5.1f | CL %4.1f | VC %5.1f | d_KSS %6.1f\n",
      nm, 3 * length(designs[[nm]]), ev$rse["KSS"], ev$rse["KINT"], ev$rse["R0"], ev$rse["KSYN"],
      ev$rse["CL"], ev$rse["VC"], ev$rse["d_KSS"]))
}
cat("Compare with the EXPECTED RESULTS block in the header.\n")
