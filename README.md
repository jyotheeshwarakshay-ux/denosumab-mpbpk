# denosumab-mpbpk

**A from-source reconstruction (with one documented simplification) of the Shah & Betts (2012) platform whole-body PBPK model for monoclonal antibodies, built up tissue-by-tissue and validated against a published bone biodistribution benchmark, ahead of application to denosumab.**

This repository transcribes and rebuilds the 15-tissue platform antibody PBPK
model in `rxode2`, verifying each structural piece (mass conservation, shared
lymph-node pooling, liver portal-vein routing) against the original source
before trusting it. It is a companion to a connected pharmacometrics programme
centered on denosumab (see [Related work](#related-work)), though — see
[The denosumab thread](#the-denosumab-thread) — it is not yet wired into that
chain.

---

## Key finding

A pure convective distribution skeleton — vascular and interstitial
compartments only, **no FcRn, no cellular/endothelial compartment** — already
reproduces Bone's antibody biodistribution coefficient (ABC) to within
**1.28-fold of the published benchmark**: this model gives **5.70%**, against
Shah & Betts' (2013) reported **7.27%**, comfortably inside the 2-fold accuracy
band the source literature itself uses for cross-species validation.

The result is verified, not just computed once and trusted:

- **Mass conservation** holds to floating-point precision (relative range
  ~1e-11) across all 32 compartments, at every timepoint of a 20,000-hour
  simulation — not just at the end.
- **The shared lymph node** (which pools lymphatic outflow from all 15
  tissues) balances inflow against outflow exactly (relative difference
  `0`–`2e-16`) at steady state.
- **Liver's portal-vein routing** is independently checked: total blood
  entering liver's vascular compartment (arterial + the four portal organs)
  exactly equals total blood leaving it, at steady state.
- **Steady-state interstitial:plasma ratios reproduce a closed-form
  prediction exactly.** Because `PLQ >> LF` everywhere, the equilibrium ratio
  for a "plain" tissue reduces analytically to `(1-sigma_V)/(1-sigma_IS)`.
  Checking that formula against each tissue's `sigma_V` reproduces every
  simulated ratio to 3–4 decimal places — strong independent evidence the
  32-equation system generated programmatically is actually correct, not
  merely stable.

The residual 1.28-fold gap to the 7.27% benchmark is attributed to the
compartments this skeleton deliberately omits (endothelial/endosomal space,
FcRn recycling) — see [Limitations](#limitations).

---

## What this is

The core reference is Shah DK & Betts AM (2012), a 15-tissue whole-body PBPK
model for monoclonal antibodies (plasma → vascular → interstitial → shared
lymph node → plasma, per tissue), reproduced here from the Berkeley Madonna
source code and parameter tables published as the supplement to Jones et al.
(2019) — see [Sources](#sources).

Built in six staged, independently-verified steps:

| Stage | Script | What it establishes |
|-------|--------|---------------------|
| 1 | `R/01_parameters.R` | Transcribes the human physiological parameter set (tissue volumes, plasma/lymph flows, reflection coefficients) from the Jones et al. 2019 source code, verified line-by-line against the paper's own Table S2. |
| 2 | `R/02_single_tissue_prototype.R` | 4-compartment Bone-only skeleton (Central → BoneV → BoneIS → LymphNode → Central); establishes the convective transport mechanism and confirms it is closed (mass-conserving) before scaling up. |
| 3 | `R/03_two_tissue_prototype.R` | Extends to Bone + Muscle sharing ONE lymph node; confirms the pooled-node balance equation and surfaces Bone's slow (~456 h) interstitial time constant. |
| 4 | `R/04_generic_14_tissues.R` | All 15 organs except Liver, built **programmatically** (looped model-text generation, not hand-copied blocks) so the framework is reusable rather than bespoke per tissue. |
| 5 | `R/05_full_15_tissues.R` | Adds Liver with correct portal-vein routing (Spleen/Pancreas/SI/LI drain through Liver before reaching plasma) — the one genuinely different structural feature in the whole model. |
| 6 | `R/06_bone_abc_validation.R` | Computes the antibody biodistribution coefficient per tissue and validates Bone against the 7.27% published benchmark (the [key finding](#key-finding)). |

Distribution mechanics only, throughout — no target binding, no FcRn-mediated
recycling, no clearance. That is a deliberate scope boundary, not an
oversight; see [Limitations](#limitations).

---

## The denosumab thread

This repository is a **mechanistic tissue-distribution PBPK model** — a
structurally different characterization of denosumab disposition from the
**population TMDD-QSS PK model** used elsewhere in this programme
(`denosumab-tmdd-qss`, and downstream in `rankl-opg-qsp`). The two are not
interchangeable: the TMDD-QSS model captures target-mediated elimination
(RANKL binding reducing clearance) with no tissue resolution, while this
model captures whole-body tissue distribution with no target binding at all
yet. Neither currently subsumes the other.

The connection this repository is building toward:

**Shah & Betts platform PBPK (this repo) → bone-tissue antibody concentration → [not yet connected] → `rankl-opg-qsp`'s RANKL-blockade input**

A prior fact-finding check of `rankl-opg-qsp/stage3_pk_coupling.R` confirmed
that project currently drives its RANKL-blocking term (`I_O`) off **plasma**
concentration (`I_O = k_scale * C_deno`, where `C_deno` is a central-PK-
compartment concentration), with no bone-tissue-specific exposure term at
all — a placeholder already flagged as "illustrative, not calibrated" in that
repo's own README. This repository's bone-compartment output is the natural
candidate to substitute in for that plasma proxy, once FcRn is added here and
denosumab-specific binding is layered on top. That coupling has not been
built yet.

---

## Diagnostic highlights

Several steps here surfaced a real discrepancy or bug rather than confirming
what was expected on the first try:

- **A structural bug in the very first lymph-node design.** The initial
  single-tissue prototype returned interstitial lymph flow directly to
  plasma. Re-reading the source line-by-line (`sigma_IS` filter, shared
  lymph-node compartment) showed this was wrong on two counts: a missing
  `(1-sigma_IS)` filter on the interstitial outflow, and a missing
  intermediate lymph-node compartment that pools *all* tissues' outflow
  before it reaches plasma. Both were rebuilt and re-verified for mass
  conservation before any further tissue was added.

- **A "failing" check that wasn't actually a bug.** The two-tissue lymph-node
  balance check failed at the chosen 1e-6 tolerance at t=2000h. Rather than
  loosen the tolerance to force a pass, the discrepancy was traced
  analytically: Bone's interstitial time constant is ~456 h, so 2000h is only
  ~4.4 time constants — not yet equilibrated. Re-running to 20,000h
  (~44 time constants) closed the balance to `2.4e-16`, confirming a slow
  transient, not a leak. The simulation window for every subsequent file was
  set to 20,000h for this reason.

- **A discrepancy between the paper and its own supplementary code.** The
  source's Berkeley Madonna file hardcodes `PLQ[Other] = 9190/1000` mL/h, but
  Table S2 of the same paper lists `Other` plasma flow as `5521` mL/h. Rather
  than silently pick one, the 9190 figure was checked and found to be exactly
  the residual that balances total peripheral flow to cardiac output
  (`181913 - sum(other 13 organs) = 9190`) — genuinely used by the model,
  genuinely inconsistent with its own published table. Documented in
  `R/01_parameters.R` rather than resolved either way.

- **Portal-routing conservation was proven before being run.** Each gut
  organ's vascular compartment drains through two separate paths from a
  single `PLQ` arterial inflow: `(PLQ-LF)` to Liver's portal inflow, and
  `LF*(1-sigma_V)` to its own interstitial space — a dual-drain structure
  that looked, on first read, like it might leak mass. Term-by-term tracing
  showed every outflow is matched by an identical inflow elsewhere in the
  system (Liver's portal-inflow term exactly matches what each gut organ
  loses to the portal path; the "missing" `LF*sigma_V` fraction is matched
  separately by that organ's own interstitial gain, not by the portal
  path). The numerical result (`0` relative difference at steady state)
  confirmed the by-hand proof.

---

## Limitations

Stated plainly, because knowing where a model stops being trustworthy is part
of using it.

1. **No FcRn, no cellular/endothelial compartment.** This is a pure
   distribution skeleton. The 5.70% Bone ABC is *expected* to sit below the
   7.27% benchmark for exactly this reason — the full source model's ABC
   includes antibody held in the endosomal/endothelial machinery that this
   model doesn't represent yet. The 1.28-fold gap is consistent with, not
   contrary to, that omission.

2. **Lung is treated as a generic tissue rather than a distinct arterial
   pass-through.** The source routes every tissue's arterial supply through
   a separate Lung vascular compartment (Central → Lung → tissues), whereas
   this model feeds tissues from central plasma directly. This
   simplification was measured, not assumed: at steady state, Lung's
   vascular concentration differs from central plasma by 0.01% (ratio
   0.9999), because Lung carries the entire cardiac output and its lymph
   filtering is negligible against that flow. The simplification therefore
   changes tissue biodistribution coefficients by ~0.01% — far below the
   model's validation precision — and was retained deliberately rather than
   rebuilt.

3. **The `PLQ[Other]` discrepancy is unresolved by design** (see
   [Diagnostic highlights](#diagnostic-highlights)) — this model uses the
   source code's `9190` mL/h, not the paper's own Table S2 value of `5521`
   mL/h, and that conflict has not been adjudicated.

4. **Single bolus, single steady-state timepoint.** No repeat dosing, no
   time-course comparison against any clinical or preclinical denosumab
   concentration-time data. The 7.27% validation is a single equilibrium
   number, not a full PK profile match.

5. **Human parameters only.** No mouse/rat/monkey cross-species validation
   has been performed here, despite that being a central feature of the
   original Shah & Betts (2012) paper.

6. **Not yet coupled to `rankl-opg-qsp`.** See
   [The denosumab thread](#the-denosumab-thread) — the bone-tissue
   concentration this model can produce has no live connection to that
   project's RANKL-blockade term yet.

---

## Reproducing

Requires R with `rxode2`.

```r
# from the repository root
source("R/01_parameters.R")
source("R/02_single_tissue_prototype.R")
source("R/03_two_tissue_prototype.R")
source("R/04_generic_14_tissues.R")
source("R/05_full_15_tissues.R")
source("R/06_bone_abc_validation.R")
```

Each later script sources the one before it, so running `06` alone
reproduces the full chain and prints every check along the way.

---

## Sources

- **Platform PBPK model:** Shah DK, Betts AM. *Towards a platform PBPK model
  to characterize the plasma and tissue disposition of monoclonal antibodies
  in preclinical species and human.* J Pharmacokinet Pharmacodyn.
  2012;39(1):67–86. Structure and physiological-parameter design of the
  15-tissue model reproduced here.
- **Bone ABC benchmark:** Shah DK, Betts AM. *Antibody biodistribution
  coefficients: Inferring tissue concentrations of monoclonal antibodies
  based on the plasma concentrations in several preclinical species and
  human.* mAbs. 2013;5(2):297–305. Source of the 7.27% (CV 11.8%) bone ABC
  benchmark validated against in `R/06_bone_abc_validation.R`.
- **Source code and parameter tables:** Jones HM, Zhang Z, Jasper P, Luo H,
  Avery LB, King LE, Neubert H, Barton HA, Betts AM, Webster R. *A
  Physiologically-Based Pharmacokinetic Model for the Prediction of
  Monoclonal Antibody Pharmacokinetics From In Vitro Data.* CPT
  Pharmacometrics Syst Pharmacol. 2019;8(10):738–747. Supplementary files
  `PSP4-8-738-s001.pdf` (Tables S1–S3, PBPK model equations) and
  `PSP4-8-738-s002.txt` (Berkeley Madonna implementation) in `sources/`. All
  physiological parameters in `R/01_parameters.R` are transcribed directly
  from these files, verified line-by-line against the paper's own Table S2.

## Related work

- **rankl-opg-qsp** — virtual-population QSP model of denosumab's effect on
  bone remodeling; currently drives its RANKL-blockade term off plasma
  concentration, a placeholder this repository's bone-compartment output
  could eventually replace.
- **denosumab-tmdd-qss** — the denosumab TMDD / QSS parameter-recovery and
  identifiability study whose population PK model feeds `rankl-opg-qsp`
  Stage 3; structurally distinct from the mechanistic PBPK approach here.
- **tmdd-approximations** — TMDD model with QSS and Michaelis–Menten
  approximations.
- **Theophylline-popPK** — population-PK reproduction.
