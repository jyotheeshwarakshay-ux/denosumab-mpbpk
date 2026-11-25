# ===========================================================================
# PAPER 2 - carry-overs before Step 5.  Read this block first.
#
#   No new model. This script prints, from R, three sets of numbers that
#   were so far computed only in Python, so they can enter the manuscript:
#   1. RANKL-mediated vs linear elimination at ONE common day (day 400,
#      60 mg SC, % of absorbed dose) for: full model with plasma 3.126 L
#      (Step 2a) and 1.4121 L (Step 4a), their lumped 2- and 3-compartment
#      versions (Step 4b), and RANKL in bone (Step 2b).
#   2. Step 2b (RANKL in bone): bone total/plasma ratio, share of bone drug
#      bound to RANKL, and total bone RANKL at days 7, 28, 84, 140.
#   3. Step 4b: when the lumped/full FREE ratio is lowest at 210 mg, the
#      full-model concentration there, and the share of time points off
#      by more than 10% (window: day >= 1 while full free > 0.01 ug/mL).
#
#   EXPECTED RESULTS (independent Python, LSODA, rtol 1e-10, atol 1e-14):
#   1. via RANKL / linear (%):
#      2a full plasma 3.126   54.1 / 45.9 | lumped 2-cmt 54.0 / 46.0 | 3-cmt 54.1 / 45.9
#      4a full plasma 1.4121  28.1 / 71.9 | lumped 2-cmt 27.9 / 72.1 | 3-cmt 28.1 / 71.9
#      2b full RANKL in bone  20.1 / 79.7  (about 0.2% still in the body)
#   2. day 7: ratio 0.017, bound 93.8%, RANKL 24.25 nM
#      day 28: 0.021, 92.0%, 18.45 | day 84: 0.021, 90.6%, 15.30
#      day 140: 0.021, 92.1%, 18.28
#   3. plasma 3.126:  2-cmt min 0.798 at day 235.3, full 0.010 ug/mL, 3.1% off
#                     3-cmt min 0.987 at day 235.3, 0.0% off
#      plasma 1.4121: 2-cmt min 0.549 at day 239.2, full 0.010 ug/mL, 12.0% off
#                     3-cmt min 0.977 at day 239.2, 0.0% off
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/05_full_15_tissues.R")

cat("\n\n=======================================================================\n")
cat("=== R/19_carry_overs.R - Paper 2 carry-overs ===\n")
cat("=======================================================================\n\n")

MW <- 147000; ug_per_mL <- MW * 1e-6
CL_lin <- 3.06 / 1000; F_sc <- 0.64; ka <- 0.00883            # Sutjandra 2011
Kss <- 1.56; R0 <- 15.23; kint <- 0.022; ksyn <- 0.01          # Choi 2025
kdeg <- ksyn / R0
d60 <- 60e-3 / MW * 1e9
fast_org <- c("Kidney", "Lung", "Thymus", "SI", "Pancreas", "Heart", "Spleen", "LI", "Liver", "Brain")
mid_org  <- c("Other", "Skin", "Muscle", "Adipose")
qss_free <- function(Ct, R) { b <- Ct - R - Kss; 0.5 * (b + sqrt(b^2 + 4 * Kss * Ct)) }
base_p <- c(CL_lin = CL_lin, F_sc = F_sc, ka = ka, Kss = Kss, kint = kint, kdeg = kdeg)

# --- same matrix build and lumping as R/18 --------------------------------
build_A <- function(V_pl) {
  nm <- c("Central", "LymphNode", as.vector(rbind(paste0(organs, "V"), paste0(organs, "IS"))))
  vol <- setNames(numeric(length(nm)), nm)
  vol["Central"] <- V_pl; vol["LymphNode"] <- V_LN
  for (o in organs) { vol[paste0(o, "V")] <- V_V[o]; vol[paste0(o, "IS")] <- V_IntS[o] }
  A <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  flow <- function(rate, from, to) {
    k <- rate / vol[from]; A[from, from] <<- A[from, from] - k; A[to, from] <<- A[to, from] + k
  }
  gut <- c("Spleen", "Pancreas", "SI", "LI")
  for (o in organs) {
    flow(PLQ[o], "Central", paste0(o, "V"))
    flow(LF[o] * (1 - sigma_V[o]), paste0(o, "V"), paste0(o, "IS"))
    flow(LF[o] * (1 - sigma_IS[o]), paste0(o, "IS"), "LymphNode")
  }
  for (o in setdiff(organs, c(gut, "Liver"))) flow(PLQ[o] - LF[o], paste0(o, "V"), "Central")
  for (g in gut) flow(PLQ[g] - LF[g], paste0(g, "V"), "LiverV")
  flow((PLQ["Liver"] - LF["Liver"]) + sum(PLQ[gut] - LF[gut]), "LiverV", "Central")
  flow(sum(LF[organs]), "LymphNode", "Central")
  A
}
equilibrium <- function(A) { s <- svd(A); x <- abs(s$v[, which.min(s$d)]); names(x) <- rownames(A); x / sum(x) }
lump <- function(A, x, V_pl, groups) {
  Cpe <- unname(x["Central"] / V_pl); G <- names(groups)
  L <- matrix(0, length(G), length(G), dimnames = list(G, G))
  V <- vapply(groups, function(g) sum(x[g]) / Cpe, numeric(1))
  for (a in G) for (b in G) {
    w <- x[groups[[b]]] / sum(x[groups[[b]]])
    L[a, b] <- sum(A[groups[[a]], groups[[b]], drop = FALSE] %*% w)
  }
  list(L = L, V = V)
}

# --- full model with elimination tracking ---------------------------------
i_cc <- grep("^C_central <-", lines); i_ce <- grep("^d/dt\\(Central\\)", lines)
nl <- lines
nl[i_cc] <- paste("Ctot <- Central / V_central", "bq <- Ctot - Rtot - Kss",
                  "C_central <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot))",
                  "RC <- Ctot - C_central", sep = "\n")
nl[i_ce] <- paste(lines[i_ce], "- CL_lin * C_central + ka * Depot - kint * RC * V_central")
nl <- c("Rtot(0) <- R0", "f(Depot) <- F_sc", nl, "d/dt(Depot) <- -ka * Depot",
        "d/dt(Rtot) <- ksyn - kdeg * Rtot - (kint - kdeg) * RC",
        "d/dt(TMDD_elim) <- kint * RC * V_central", "d/dt(Lin_elim) <- CL_lin * C_central")
mod_full <- rxode2(paste(nl, collapse = "\n"))

mod_l2 <- rxode2({
  Ctot <- Ac / Vc; bq <- Ctot - Rtot - Kss
  C <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot)); RC <- Ctot - C; Acf <- C * Vc
  d/dt(Depot) <- -ka * Depot
  d/dt(Ac) <- L11 * Acf + L12 * Ap + ka * Depot - kint * RC * Vc - CL_lin * C
  d/dt(Ap) <- L21 * Acf + L22 * Ap
  d/dt(Rtot) <- ksyn_a - kdeg * Rtot - (kint - kdeg) * RC
  d/dt(TMDD_elim) <- kint * RC * Vc
  d/dt(Lin_elim) <- CL_lin * C
  Rtot(0) <- R0_a
  f(Depot) <- F_sc
})
mod_l3 <- rxode2({
  Ctot <- Ac / Vc; bq <- Ctot - Rtot - Kss
  C <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot)); RC <- Ctot - C; Acf <- C * Vc
  d/dt(Depot) <- -ka * Depot
  d/dt(Ac) <- L11 * Acf + L12 * Ai + L13 * Ab + ka * Depot - kint * RC * Vc - CL_lin * C
  d/dt(Ai) <- L21 * Acf + L22 * Ai + L23 * Ab
  d/dt(Ab) <- L31 * Acf + L32 * Ai + L33 * Ab
  d/dt(Rtot) <- ksyn_a - kdeg * Rtot - (kint - kdeg) * RC
  d/dt(TMDD_elim) <- kint * RC * Vc
  d/dt(Lin_elim) <- CL_lin * C
  Rtot(0) <- R0_a
  f(Depot) <- F_sc
})

lump_params <- function(V_pl, groups) {
  A <- build_A(V_pl); x <- equilibrium(A); lp <- lump(A, x, V_pl, groups)
  Vc <- unname(lp$V["C"])
  p <- c(base_p, Vc = Vc, R0_a = R0 * V_pl / Vc, ksyn_a = ksyn * V_pl / Vc)
  for (i in seq_len(nrow(lp$L))) for (j in seq_len(ncol(lp$L))) p[sprintf("L%d%d", i, j)] <- lp$L[i, j]
  p
}
central <- c("Central", "LymphNode", paste0(organs, "V"), paste0(fast_org, "IS"))
g2 <- list(C = central, P = c(paste0(mid_org, "IS"), "BoneIS"))
g3 <- list(C = central, I = paste0(mid_org, "IS"), B = "BoneIS")
pct <- function(r) 100 * c(r$TMDD_elim[nrow(r)], r$Lin_elim[nrow(r)]) / (F_sc * d60)

# --- 1. Elimination at day 400 --------------------------------------------
ev400 <- et(amt = d60, cmt = "Depot", time = 0) |> et(c(0, 24 * 400))
cat("=== 1. Elimination by day 400, 60 mg SC: via RANKL / linear (% of absorbed dose) ===\n")
for (V_pl in c(3.126, 3.126 - sum(V_V[organs]))) {
  pf <- c(organ_params, base_p, R0 = R0, ksyn = ksyn); pf["V_central"] <- V_pl
  e <- pct(rxSolve(mod_full, pf, ev400, atol = 1e-14, rtol = 1e-10))
  cat(sprintf("  full, plasma %.4f L : %5.1f / %5.1f\n", V_pl, e[1], e[2]))
  e <- pct(rxSolve(mod_l2, lump_params(V_pl, g2), ev400, atol = 1e-14, rtol = 1e-10))
  cat(sprintf("     lumped 2-cmt       : %5.1f / %5.1f\n", e[1], e[2]))
  e <- pct(rxSolve(mod_l3, lump_params(V_pl, g3), ev400, atol = 1e-14, rtol = 1e-10))
  cat(sprintf("     lumped 3-cmt       : %5.1f / %5.1f\n", e[1], e[2]))
}

# --- Step 2b model (RANKL in bone), as R/15 --------------------------------
V_BIS <- unname(organ_params["Vinter_Bone"])
R0_bone <- 15.23 * V_central / V_BIS; ksyn_bone <- kdeg * R0_bone
i_cb <- grep("^C_BoneIS", lines); i_bi <- grep("^d/dt\\(BoneIS\\)", lines)
nb <- lines
nb[i_cb] <- paste("CtotB <- BoneIS / Vinter_Bone", "bqB <- CtotB - RtotB - Kss",
                  "C_BoneIS <- 0.5 * (bqB + sqrt(bqB * bqB + 4 * Kss * CtotB))",
                  "RCB <- CtotB - C_BoneIS", sep = "\n")
nb[i_bi] <- paste(lines[i_bi], "- kint * RCB * Vinter_Bone")
nb[i_ce] <- paste(lines[i_ce], "- CL_lin * C_central + ka * Depot")
nb <- c("RtotB(0) <- R0_bone", "f(Depot) <- F_sc", nb, "d/dt(Depot) <- -ka * Depot",
        "d/dt(RtotB) <- ksyn_bone - kdeg * RtotB - (kint - kdeg) * RCB",
        "d/dt(TMDD_elim) <- kint * RCB * Vinter_Bone", "d/dt(Lin_elim) <- CL_lin * C_central")
mod_bone <- rxode2(paste(nb, collapse = "\n"))
pb <- c(organ_params, base_p, R0_bone = R0_bone, ksyn_bone = ksyn_bone)
e <- pct(rxSolve(mod_bone, pb, ev400, atol = 1e-12, rtol = 1e-10))
cat(sprintf("  full, RANKL in bone  : %5.1f / %5.1f\n", e[1], e[2]))
cat("[Python: 54.1/45.9, 54.0/46.0, 54.1/45.9 | 28.1/71.9, 27.9/72.1, 28.1/71.9 | 20.1/79.7]\n\n")

# --- 2. Step 2b bone figures ------------------------------------------------
cat("=== 2. RANKL in bone: bone total/plasma, bound share in bone, total bone RANKL ===\n")
rb <- rxSolve(mod_bone, pb, et(amt = d60, cmt = "Depot", time = 0) |> et(c(7, 28, 84, 140) * 24),
              atol = 1e-12, rtol = 1e-10)
for (k in seq_len(nrow(rb))) {
  CtB <- rb$BoneIS[k] / V_BIS; CfB <- qss_free(CtB, rb$RtotB[k]); Cp <- rb$Central[k] / V_central
  cat(sprintf("  day %3d: ratio %.3f | bound %.1f%% | RANKL %.2f nM\n",
              rb$time[k] / 24, CtB / Cp, 100 * (CtB - CfB) / CtB, rb$RtotB[k]))
}
cat("[Python: 0.017/93.8/24.25 | 0.021/92.0/18.45 | 0.021/90.6/15.30 | 0.021/92.1/18.28]\n\n")

# --- 3. Timing of the lumped-model deviation at 210 mg ----------------------
cat("=== 3. Lumped/full FREE ratio at 210 mg: minimum, when, full conc there, % points off >10% ===\n")
tt <- seq(0, 24 * 300, by = 1)
ev210 <- et(amt = 210e-3 / MW * 1e9, cmt = "Depot", time = 0) |> et(tt)
for (V_pl in c(3.126, 3.126 - sum(V_V[organs]))) {
  pf <- c(organ_params, base_p, R0 = R0, ksyn = ksyn); pf["V_central"] <- V_pl
  rf <- rxSolve(mod_full, pf, ev210, atol = 1e-14, rtol = 1e-10)
  a <- qss_free(rf$Central / V_pl, rf$Rtot)
  for (lab in c("2-cmt", "3-cmt")) {
    p <- if (lab == "2-cmt") lump_params(V_pl, g2) else lump_params(V_pl, g3)
    m <- if (lab == "2-cmt") mod_l2 else mod_l3
    r <- rxSolve(m, p, ev210, atol = 1e-14, rtol = 1e-10)
    b <- qss_free(r$Ac / p["Vc"], r$Rtot)
    w <- tt >= 24 & a * ug_per_mL > 0.01
    q <- b[w] / a[w]; i <- which.min(q)
    cat(sprintf("  plasma %.4f %s: min %.3f at day %.1f (full %.3f ug/mL) | off >10%%: %.1f%%\n",
                V_pl, lab, q[i], tt[w][i] / 24, a[w][i] * ug_per_mL, 100 * mean(abs(q - 1) > 0.10)))
  }
}
cat("[Python: 0.798 d235.3 0.010 3.1% | 0.987 d235.3 0.0% | 0.549 d239.2 0.010 12.0% | 0.977 d239.2 0.0%]\n")
