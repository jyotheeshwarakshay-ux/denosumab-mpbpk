# denosumab-mpbpk — instructions for Claude (Claude Code and Claude Science)

This repository is the code base for **Paper 2** (direction locked 26 Nov 2026; target PhD: Uppsala University, backup: Leiden).

## Paper 2 in one line
Reduce a whole-body denosumab PBPK–TMDD model to the clinical 2-compartment QSS-TMDD model, explain why Kss was fragile in Paper 1, then use PopED optimal design + simulation–estimation to find a sampling design that makes Kss estimable. Target journal: CPT:PSP (fallback JPKPD). Preprint by 30 June 2027.

## Stages
- S1 DONE (commit a5ba6fa): `R/13_whole_body_elimination.R` — linear plasma elimination, CL 3.06 mL/h (Sutjandra 2011); terminal half-life 56.6 days; Vss 5.865 L.
- S2a DONE (commit 6825cf6): `R/14_sc_dosing_rankl_plasma.R` — SC (F 0.64, ka 0.00883 /h, Sutjandra 2011) + plasma RANKL QSS (Kss 1.56 nM, R0 15.23 nM, kint 0.022 /h, ksyn 0.01 nM/h, Choi 2025); half-life 24.1 days vs label 25.4.
- S2b DONE (commit 8336af5): `R/15_rankl_in_bone.R` — same RANKL amount in bone interstitium; half-life 41.8 days; bone RANKL suppression delivery-limited (~92% of bone drug bound; bone total/plasma 0.02 vs 0.188 without target).
- S3 DONE (verified in R, all 21 values match Python): `R/16_compare_clinical_model.R` — PBPK vs Choi 2025 2-cmt QSS-TMDD at 6/60/210 mg SC (Sutjandra full text unavailable). Python: mismatch has two causes — central volume (3.126 L vs 1.58 L apparent) at low dose, distribution volume (Vss 5.87 L) at high dose. Do not tune.
- S4a DONE (commit 38c3f51): R/17_total_plasma_volume.R — central plasma 1.4121 L so total plasma = 3.126 L (Jones 2019 code adds 1.714 L organ vascular on top of 3.126 L; Chang et al. 2022 Table 1 totals 3.135 L). Vss 4.151 L; 60 mg SC Cmax 7.06 ug/mL (label 6.75), half-life 29.4 d (label 25.4), AUC0-16wk 374 (label 316); PBPK/Choi median ratio 0.20 / 1.08 / 1.20 at 6 / 60 / 210 mg. Sensitivity analysis, not a claim that Jones 2019 is wrong.
- S4b DONE: R/18_lumped_models.R — lumped Vc 4.999 L (plasma 3.126) / 3.285 L (plasma 1.4121), Vp 0.866 L, Q 7.3 mL/h; 3-cmt lumped vs full PBPK free conc within 2.3% at 6/60/210 mg; 2-cmt median within 2.3% but down to 0.55x late at 210 mg (bone merged with intermediate tissues); R0_app = R0 x V_plasma / Vc. NEXT: carry-overs (RANKL elimination at a common day; Step 2b bound-fraction figures), then Step 5 PopED.
- Carry-overs DONE: R/19_carry_overs.R — RANKL elimination at day 400 (60 mg): 54.1% (plasma 3.126), 28.1% (plasma 1.4121), 20.1% (RANKL in bone); lumped models within 0.2 points. Step 2b: 90.6–93.8% of bone drug bound, bone total/plasma 0.017–0.021. 2-cmt deviation at 210 mg is a terminal-tail effect (min at the 0.010 ug/mL cutoff, day 235–239); use '3–12% of time points off >10%' in text; 3-cmt 0%. NEXT: Step 5 PopED optimal design for Kss.
- S5a DONE (commit c7a08cc): R/20_poped_baseline.R — PopED FIM (FO) for Paper 1 designs, CRAN PopED 0.7.0 is the reference: Kss RSE v2 78.9% vs v1 35.1% (v1 keeps late 3600/4380 h samples where free drug nears Kss). IIV of Kss/kint not estimable in either design. NEXT S5b: FOCE check; optimise sampling for Kss under assay LLOQ 0.136 nM; plan in paper2/step5-poped-notes.md (claude.ai project).
- S5b DONE: R/21_poped_candidate_designs.R — free drug < Kss from ~3400 h after each dose; design D (same 24 samples, 3 moved to 3000/3600/4200 h) Kss RSE 33.6% vs 78.9% for v2 (CRAN PopED 0.7.0); FOCE not feasible in workspace; NEXT Step 6 SSE (A vs D, BLQ handling).
- S5: PopED optimal sampling design for Kss.
- S6: simulation–estimation in nlmixr2 (Kss bias, precision, fit success).

## Rules (non-negotiable)
- One change per new numbered script (`R/NN_description.R`); never edit earlier validated scripts.
- Every script starts with a summary block: what changed, assumptions, expected results.
- Every number that goes into the manuscript is checked against the script output that produced it.
- Never tune a parameter to force a match with published data. Report mismatches as fold-differences with possible causes.
- Every parameter value carries its source (paper + table/line) in a comment.
- Scripts only (R: rxode2, nlmixr2, PopED). No GUI tools.
- Explain in plain English before code. One step at a time; wait for confirmation of the output before the next step.
- Only one tool (Claude Code or Claude Science) edits this folder at a time. Commit to GitHub after each validated step.

## Key references
- Sutjandra L et al. 2011, Clin Pharmacokinet 50:793–807. doi:10.2165/11594240-000000000-00000
- Fronton L, Pilari S, Huisinga W 2014, J PKPD 41:87–107. doi:10.1007/s10928-014-9349-1
- Pilari S, Huisinga W 2010, J PKPD 37:365–405. doi:10.1007/s10928-010-9165-1
- Brekkan A et al. 2018, J PKPD 45:637–647. doi:10.1007/s10928-018-9594-9
- Nyberg J et al. 2012 (PopED), Comput Methods Programs Biomed 108:789–805. doi:10.1016/j.cmpb.2012.05.005
- Lavezzi SM et al. 2018, J PKPD 45:787–802. doi:10.1007/s10928-018-9608-7 (closest prior work; cite and differentiate)

Manuscript draft: https://claude.ai/artifact/Q1bYMs5mvcwHuF6PQepaTY
