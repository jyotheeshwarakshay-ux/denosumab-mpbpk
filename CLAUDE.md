# denosumab-mpbpk — instructions for Claude (Claude Code and Claude Science)

This repository is the code base for **Paper 2** (direction locked 26 Nov 2026; target PhD: Uppsala University, backup: Leiden).

## Paper 2 in one line
Reduce a whole-body denosumab PBPK–TMDD model to the clinical 2-compartment QSS-TMDD model, explain why Kss was fragile in Paper 1, then use PopED optimal design + simulation–estimation to find a sampling design that makes Kss estimable. Target journal: CPT:PSP (fallback JPKPD). Preprint by 30 June 2027.

## Stages
- S1 (done, pending R confirmation): `R/13_whole_body_elimination.R` — linear elimination from plasma, CL = 3.06 mL/h (Sutjandra 2011). Expected: mass balance error < 1e-8, AUC·CL/dose = 1, terminal half-life ≈ 56.6 days.
- S2: subcutaneous absorption via lymph (F = 64%, ka = 0.00883 /h, Sutjandra 2011) + RANKL binding. RANKL location and turnover must come from cited primary sources.
- S3: compare with Sutjandra 2011 typical profiles at ≥3 doses.
- S4: reduction to 2-compartment QSS-TMDD; lumped parameters vs Sutjandra.
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
