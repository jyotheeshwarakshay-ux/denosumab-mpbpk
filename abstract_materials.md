# Abstract materials — denosumab-mpbpk

Reference document for writing a PAGE abstract. This is a collection of
verified numbers, facts, and citations pulled from the actual committed
files and their outputs — **not the abstract itself**, and not written from
memory. Every number below was re-run from its source file on the current
commit before being written down here; where something could not be traced
to a file, it is marked `[UNVERIFIED]` rather than included silently.

---

## The claim (locked scope)

PAGE abstract, this repo only. The finding is a **decomposition**, not a
single validated prediction:

> A from-source antibody PBPK reproduction predicts denosumab bone exposure
> to ~78% of the published benchmark by convective distribution alone; the
> residual is attributable to FcRn kinetics, not tissue volume — established
> by ruling out the volume hypothesis.

---

## Key numbers (with source file for each)

| Number | Value | Source |
|---|---|---|
| Validated distribution Bone ABC | **5.701%** (reported rounded as 5.70%) | `R/06_bone_abc_validation.R`, re-run: `Bone ABC (this distribution-only skeleton): 5.701 %` |
| Published benchmark | **7.27%** (CV 11.8%) | Shah & Betts (2013), *mAbs* 5(2):297–305, cited in `R/06_bone_abc_validation.R` and README Sources |
| Fold difference | **1.28-fold** lower than benchmark | `R/06_bone_abc_validation.R` output: `Ratio (this model / benchmark): 0.784` → `1.28-fold LOWER` |
| "~78% of benchmark" | 5.701 / 7.27 = **0.7841** → ~78% | Computed from the two numbers above; not printed verbatim by any script as a percentage, so this specific framing is a derived figure, not a direct file quote — flagged here rather than presented as a script output |
| Mass conservation precision, 15-tissue distribution model | Relative range **9.168e-12** across 32 compartments, every timepoint of a 20,000 h run | `R/06_bone_abc_validation.R` re-run (which sources `R/05_full_15_tissues.R`): `Absolute range: 9.168e-12` / `Relative range: 9.168e-12` / `Mass conserved (rel. range < 1e-8)? TRUE` |
| Minimal E6a rescue at correct (whole-body) FcRn scarcity | 5.70% → **5.8587%** (~0.16 percentage points) | `R/11_bone_fcrn_shared_pool.R` re-run: `Bone ABC, shared FcRn pool (this file, 11): 5.8587 %`. Delta computed as 5.8587 − 5.701 = 0.1577 ≈ 0.16pp |
| The debunked near-match | **7.2764%** at bone-only-scoped FcRn pool | `R/10_bone_fcrn_real_dose.R` re-run: `Bone ABC at REAL dose (this file, 60 mg): 7.2764 %` |
| ...corrected to | **5.8587%** at whole-body-scoped FcRn pool | `R/11_bone_fcrn_shared_pool.R` re-run (same figure as rescue row above — same run) |
| FcRn pool overstatement factor | bone-only pool was **6.7659×** too concentrated (298.327 µM vs 44.0927 µM) | `R/11_bone_fcrn_shared_pool.R` re-run: `FcRn_Conc, bone-only (old): 298.327 uM` / `FcRn_Conc, whole-body (new): 44.0927 uM` / `Ratio (old / new): 6.7659` |
| Cellular-volume bounding estimates | Lower bound **20.3147%**, upper bound **82.7599%**, both overshoot 7.27% | `R/12_cellular_volume_diagnostic.R` re-run: `Lower bound (cellular = IS conc): 20.3147 %` / `Upper bound (cellular = plasma conc): 82.7599 %` |
| Real denosumab dose → plasma concentration | 60 mg → **0.13057 µM** (`C0_plasma_uM`) | `R/10_bone_fcrn_real_dose.R` re-run: `C0_plasma_uM = ... = 0.13057 uM` |

---

## Model provenance

- **Base model:** Shah & Betts (2012) 15-tissue whole-body platform PBPK
  model for monoclonal antibodies, reproduced from the Berkeley Madonna
  source code and parameter tables published as the supplement to Jones et
  al. (2019), *CPT: Pharmacometrics & Systems Pharmacology* — files
  `sources/PSP4-8-738-s001.pdf` (Tables S1–S3, model equations) and
  `sources/PSP4-8-738-s002.txt` (reference code), per README Sources and
  `R/01_parameters.R`.
- **Implementation:** built in `rxode2` (R). 15 tissues, with liver
  portal-vein routing (Spleen/Pancreas/SI/LI drain through Liver before
  reaching plasma, `R/05_full_15_tissues.R`) and one shared lymph node
  pooling lymphatic outflow across all 15 tissues (`R/04_generic_14_tissues.R`,
  `R/05_full_15_tissues.R`).
- **One documented simplification** in the validated distribution model:
  Lung is treated as a generic tissue rather than a distinct arterial
  pass-through compartment. Measured, not assumed — at steady state, Lung's
  vascular concentration differs from central plasma by **0.01%** (ratio
  0.9999), per README Limitations §2. (This is the only simplification
  remaining in the validated model as of the current README; the earlier
  plain-tissue `PLQ` vs `PLQ-LF` simplification was corrected and removed —
  see `35d563f` in git history / README Limitations discussion.)

---

## Methods summary (for the Methods section later)

- **Staged, independently-verified build**, each stage checked before the
  next was added:
  1. `R/01_parameters.R` — human physiological parameters transcribed from
     Jones et al. 2019, verified line-by-line against the paper's own Table
     S2.
  2. `R/02_single_tissue_prototype.R` — 4-compartment Bone-only skeleton;
     confirms the convective transport mechanism is mass-conserving.
  3. `R/03_two_tissue_prototype.R` — Bone + Muscle sharing one lymph node;
     confirms the pooled-node balance equation.
  4. `R/04_generic_14_tissues.R` — 14 of 15 organs, built programmatically
     (looped model-text generation).
  5. `R/05_full_15_tissues.R` — adds Liver with portal-vein routing.
  6. `R/06_bone_abc_validation.R` — computes per-tissue ABC and validates
     Bone against the 7.27% benchmark.
- **Validation approach:** steady-state (20,000 h simulation) antibody
  biodistribution coefficient (ABC) per tissue, compared against the
  published bone ABC benchmark (Shah & Betts 2013).
- **FcRn investigation** (exploratory, layered on top of the validated
  model): E6a (sorting endosome) single-site FcRn rescue for the drug only
  (`R/07`–`R/09`); real-dose conversion (`R/10`); whole-body FcRn pool
  correction (`R/11`); bulk cellular-volume bounding diagnostic (`R/12`).

---

## Citations (verified titles, copied exact from README Sources)

- Shah DK, Betts AM. *Towards a platform PBPK model to characterize the
  plasma and tissue disposition of monoclonal antibodies in preclinical
  species and human.* J Pharmacokinet Pharmacodyn. 2012;39(1):67–86.
- Shah DK, Betts AM. *Antibody biodistribution coefficients: Inferring
  tissue concentrations of monoclonal antibodies based on the plasma
  concentrations in several preclinical species and human.* mAbs.
  2013;5(2):297–305.
- Jones HM, Zhang Z, Jasper P, Luo H, Avery LB, King LE, Neubert H, Barton
  HA, Betts AM, Webster R. *A Physiologically-Based Pharmacokinetic Model
  for the Prediction of Monoclonal Antibody Pharmacokinetics From In Vitro
  Data.* CPT Pharmacometrics Syst Pharmacol. 2019;8(10):738–747.
- **Denosumab FDA label** — MW ~147 kDa, dose 60 mg every 6 months (verified
  against the official FDA prescribing information / DailyMed; re-confirm
  the exact label URL before abstract submission).

---

## Caveats to state honestly in the abstract

- **The FcRn investigation (files 07–12) is exploratory and un-calibrated**,
  not a second validated result. It sits subordinate to the validated
  5.70% distribution finding in the README, and should be presented the
  same way in the abstract.
- **"~78% of benchmark" describes reproduction of the published number**,
  not a proven mechanistic split of where the missing ~22% goes. The
  distribution model reproduces most of the benchmark on convective
  transport alone; that is a statement about model performance, not a
  claim about the physical partition of drug within bone tissue.
- **The residual gap is "attributable to" FcRn kinetics because the volume
  hypothesis was actively ruled out** (`R/12`: bounds of 20.3% and 82.8%
  both overshoot 7.27%), and a minimal FcRn mechanism was shown not to
  close it at correct FcRn scarcity (`R/11`: only +0.16pp). This is
  elimination-by-testing, not a demonstration that FcRn kinetics **does**
  close the gap — the full FcRn machinery (all endosomal compartments,
  multi-site binding, endogenous IgG competition, degradation-linked FcRn
  turnover) was never built. Phrase as "FcRn kinetics, not volume,"
  **never** as "FcRn shown to close the gap."
- The apparent 7.28% near-match to the benchmark (`R/10`) was a
  **compensating-error artifact** (6.77× overstated bone-only FcRn pool),
  not a result — worth stating explicitly if the abstract discusses the
  investigation's history at all, to preempt the obvious "but you got
  7.28%, why isn't that the headline" question.
