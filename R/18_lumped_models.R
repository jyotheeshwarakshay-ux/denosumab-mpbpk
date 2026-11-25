# ===========================================================================
# PAPER 2 - STEP 4b: lumping the full PBPK-TMDD model into 2- and
# 3-compartment models.  Read this block first.
#
#   METHOD (Pilari & Huisinga 2010, J PKPD 37:365-405): compartments are
#   grouped by time-scale; inside a group, drug is split in the same
#   proportions as at equilibrium. This preserves Vss exactly.
#     Central  = plasma + organ vascular + lymph node + fast tissue fluid
#                (kidney, lung, thymus, SI, pancreas, heart, spleen, LI,
#                liver, brain; time constants < 10 h)
#     Intermediate = other, skin, muscle, adipose tissue fluid (56-127 h)
#     Bone     = bone tissue fluid (456 h)
#   2-cmt model: Central + (Intermediate + Bone). 3-cmt: all three.
#   The linear model is rebuilt here as a matrix (same equations as
#   R/05); Check A confirms it reproduces file 05's equilibrium.
#
#   RANKL (QSS, in plasma) carries over exactly into the lumped central
#   compartment with an APPARENT target concentration:
#     R0_apparent = R0 * V_plasma / Vc_lumped   (ksyn scaled the same way)
#   because the bound drug sits only in plasma but the lumped central
#   concentration is expressed per Vc_lumped. Kss, kint, kdeg unchanged.
#
#   Both plasma settings are run: 3.126 L (Jones 2019 code) and 1.4121 L
#   (Step 4a sensitivity, total plasma 3.126 L).
#   Comparison = FREE plasma concentration (Sutjandra 2011 measured free
#   denosumab), lumped / full, from day 1 while full > 0.01 ug/mL.
#
#   EXPECTED RESULTS (independent Python, LSODA, rtol 1e-10, atol 1e-14):
#   Plasma 3.126 L:
#     2-cmt VC 4.999, VP 0.866 L, Q 7.3 mL/h, Vss 5.865, R0_app 9.52 nM
#       6 mg median 0.999 (0.998-1.001); 60 mg 1.023 (0.979-1.095);
#       210 mg 1.008 (0.798-1.017); Cmax 0.034/4.317/18.204 vs full 0.034/4.403/18.575
#     3-cmt VC 4.999, VI 0.511, VB 0.355 L, Q_I 6.6, Q_B 0.8 mL/h
#       6 mg 1.000 (0.999-1.000); 60 mg 1.001 (0.995-1.001);
#       210 mg 1.000 (0.987-1.001); Cmax 0.034/4.397/18.548
#   Plasma 1.4121 L:
#     2-cmt VC 3.285, VP 0.866 L, Q 7.3 mL/h, Vss 4.151, R0_app 6.55 nM
#       6 mg 0.994 (0.983-1.067); 60 mg 1.022 (0.811-1.035);
#       210 mg 0.999 (0.549-1.029); Cmax 0.120/6.798/25.665 vs full 0.121/6.985/26.376
#     3-cmt: 6 mg 1.000 (0.998-1.004); 60 mg 1.000 (0.987-1.002);
#       210 mg 1.000 (0.977-1.002); Cmax 0.121/6.970/26.319
#   Published for context: Sutjandra Vc 2.49 L/66 kg, R0 614 ng/mL;
#   Gibiansky Vc 2.62, Vss 3.96 L/66 kg; Choi Vc/F 1.58, Vp/F 6.06,
#   Q/F 0.20 L/h, R0 15.23 nM.
#
# Run from the project root (denosumab-mpbpk/) so source() resolves.
# ===========================================================================

source("R/05_full_15_tissues.R")

cat("\n\n=======================================================================\n")
cat("=== R/18_lumped_models.R - Paper 2, Step 4b ===\n")
cat("=======================================================================\n\n")

MW <- 147000; ug_per_mL <- MW * 1e-6
CL_lin <- 3.06 / 1000; F_sc <- 0.64; ka <- 0.00883            # Sutjandra 2011
Kss <- 1.56; R0 <- 15.23; kint <- 0.022; ksyn <- 0.01          # Choi 2025
kdeg <- ksyn / R0

fast_org <- c("Kidney", "Lung", "Thymus", "SI", "Pancreas", "Heart", "Spleen", "LI", "Liver", "Brain")
mid_org  <- c("Other", "Skin", "Muscle", "Adipose")

# --- Linear distribution model as a matrix (amounts), same flows as file 05
build_A <- function(V_pl) {
  nm <- c("Central", "LymphNode", as.vector(rbind(paste0(organs, "V"), paste0(organs, "IS"))))
  vol <- setNames(numeric(length(nm)), nm)
  vol["Central"] <- V_pl; vol["LymphNode"] <- V_LN
  for (o in organs) { vol[paste0(o, "V")] <- V_V[o]; vol[paste0(o, "IS")] <- V_IntS[o] }
  A <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  flow <- function(rate, from, to) {
    k <- rate / vol[from]
    A[from, from] <<- A[from, from] - k
    A[to, from]   <<- A[to, from] + k
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

equilibrium <- function(A) {
  s <- svd(A); x <- abs(s$v[, which.min(s$d)]); names(x) <- rownames(A); x / sum(x)
}

lump <- function(A, x, V_pl, groups) {
  Cpe <- unname(x["Central"] / V_pl)
  G <- names(groups); k <- length(G)
  L <- matrix(0, k, k, dimnames = list(G, G))
  V <- vapply(groups, function(g) sum(x[g]) / Cpe, numeric(1))
  for (a in G) for (b in G) {
    w <- x[groups[[b]]] / sum(x[groups[[b]]])
    L[a, b] <- sum(A[groups[[a]], groups[[b]], drop = FALSE] %*% w)
  }
  list(L = L, V = V)
}

qss_free <- function(Ct, R) { b <- Ct - R - Kss; 0.5 * (b + sqrt(b^2 + 4 * Kss * Ct)) }

# --- Full model (as file 17) with a given plasma volume -------------------
i_cc <- grep("^C_central <-", lines); i_ce <- grep("^d/dt\\(Central\\)", lines)
nl <- lines
nl[i_cc] <- paste("Ctot <- Central / V_central", "bq <- Ctot - Rtot - Kss",
                  "C_central <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot))",
                  "RC <- Ctot - C_central", sep = "\n")
nl[i_ce] <- paste(lines[i_ce], "- CL_lin * C_central + ka * Depot - kint * RC * V_central")
nl <- c("Rtot(0) <- R0", "f(Depot) <- F_sc", nl, "d/dt(Depot) <- -ka * Depot",
        "d/dt(Rtot) <- ksyn - kdeg * Rtot - (kint - kdeg) * RC")
mod_full <- rxode2(paste(nl, collapse = "\n"))

# --- Lumped models (central index 1) --------------------------------------
mod_l2 <- rxode2({
  Ctot <- Ac / Vc; bq <- Ctot - Rtot - Kss
  C <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot)); RC <- Ctot - C
  Acf <- C * Vc
  d/dt(Depot) <- -ka * Depot
  d/dt(Ac) <- L11 * Acf + L12 * Ap + ka * Depot - kint * RC * Vc - CL_lin * C
  d/dt(Ap) <- L21 * Acf + L22 * Ap
  d/dt(Rtot) <- ksyn_a - kdeg * Rtot - (kint - kdeg) * RC
  Rtot(0) <- R0_a
  f(Depot) <- F_sc
})
mod_l3 <- rxode2({
  Ctot <- Ac / Vc; bq <- Ctot - Rtot - Kss
  C <- 0.5 * (bq + sqrt(bq * bq + 4 * Kss * Ctot)); RC <- Ctot - C
  Acf <- C * Vc
  d/dt(Depot) <- -ka * Depot
  d/dt(Ac) <- L11 * Acf + L12 * Ai + L13 * Ab + ka * Depot - kint * RC * Vc - CL_lin * C
  d/dt(Ai) <- L21 * Acf + L22 * Ai + L23 * Ab
  d/dt(Ab) <- L31 * Acf + L32 * Ai + L33 * Ab
  d/dt(Rtot) <- ksyn_a - kdeg * Rtot - (kint - kdeg) * RC
  Rtot(0) <- R0_a
  f(Depot) <- F_sc
})

tt <- seq(0, 24 * 300, by = 1)
base_p <- c(CL_lin = CL_lin, F_sc = F_sc, ka = ka, Kss = Kss, kint = kint, kdeg = kdeg)

for (V_pl in c(3.126, 3.126 - sum(V_V[organs]))) {
  cat(sprintf("===== central plasma %.4f L =====\n", V_pl))
  A <- build_A(V_pl); x <- equilibrium(A)
  Cpe <- unname(x["Central"] / V_pl)
  cat(sprintf("Check A (matrix vs file 05): Bone IS/plasma %.4f [0.1878], Vss %.3f L\n",
              x["BoneIS"] / V_IntS["Bone"] / Cpe, 1 / Cpe))
  central <- c("Central", "LymphNode", paste0(organs, "V"), paste0(fast_org, "IS"))
  g2 <- list(C = central, P = c(paste0(mid_org, "IS"), "BoneIS"))
  g3 <- list(C = central, I = paste0(mid_org, "IS"), B = "BoneIS")

  p_full <- c(organ_params, base_p, R0 = R0, ksyn = ksyn); p_full["V_central"] <- V_pl
  full <- lapply(c(6, 60, 210), function(dm) {
    r <- rxSolve(mod_full, p_full, et(amt = dm * 1e-3 / MW * 1e9, cmt = "Depot", time = 0) |> et(tt),
                 atol = 1e-14, rtol = 1e-10)
    qss_free(r$Central / V_pl, r$Rtot)
  })

  for (lab in c("2-cmt", "3-cmt")) {
    lp <- if (lab == "2-cmt") lump(A, x, V_pl, g2) else lump(A, x, V_pl, g3)
    L <- lp$L; V <- lp$V; Vc <- unname(V["C"])
    R0_a <- R0 * V_pl / Vc
    Qs <- sapply(names(V)[-1], function(h) L["C", h] * V[h] * 1000)
    cat(sprintf("%s: %s L | Q %s mL/h | Vss %.3f | R0_app %.2f nM\n", lab,
                paste(sprintf("V%s=%.3f", names(V), V), collapse = ", "),
                paste(sprintf("%s:%.1f", names(V)[-1], Qs), collapse = ", "), sum(V), R0_a))
    pl <- c(base_p, Vc = unname(Vc), R0_a = unname(R0_a), ksyn_a = unname(ksyn * V_pl / Vc))
    for (i in seq_len(nrow(L))) for (j in seq_len(ncol(L))) pl[sprintf("L%d%d", i, j)] <- L[i, j]
    m <- if (lab == "2-cmt") mod_l2 else mod_l3
    for (k in 1:3) {
      dm <- c(6, 60, 210)[k]
      r <- rxSolve(m, pl, et(amt = dm * 1e-3 / MW * 1e9, cmt = "Depot", time = 0) |> et(tt),
                   atol = 1e-14, rtol = 1e-10)
      b <- qss_free(r$Ac / Vc, r$Rtot); a <- full[[k]]
      w <- tt >= 24 & a * ug_per_mL > 0.01
      q <- b[w] / a[w]
      cat(sprintf("   %3d mg: lumped/full FREE median %.3f, range %.3f-%.3f | Cmax %.3f vs %.3f\n",
                  dm, median(q), min(q), max(q), max(b) * ug_per_mL, max(a) * ug_per_mL))
    }
  }
  cat("\n")
}
cat("Compare every number with the EXPECTED RESULTS block at the top of this file.\n")
