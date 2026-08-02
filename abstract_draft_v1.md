DRAFT v1 - PAGE abstract, denosumab-mpbpk. Sentences drafted by author, to
be reviewed against PAGE word limit and formatting before submission.
Numbers verified in abstract_materials.md.

## Objectives

To reproduce the Shah & Betts (2012) platform monoclonal antibody PBPK
model from its original source code and parameter tables, validate it
against a published bone tissue exposure benchmark, and investigate
whether FcRn-mediated recycling or bone cellular volume can explain any
residual discrepancy, using denosumab as the candidate antibody.

## Methods

A 15-tissue whole-body PBPK model (Shah & Betts, 2012) was reconstructed
in rxode2 from its Berkeley Madonna source code and parameter tables
(Jones et al., 2019), built and verified in staged increments
(single-tissue prototype, shared lymph node, liver portal-vein routing),
validated against the published bone antibody biodistribution coefficient
(ABC) benchmark at steady state, and then extended with an exploratory,
single-site FcRn-mediated recycling mechanism in bone's sorting endosome
and a bounding calculation testing whether bone's unmodeled cellular
volume could account for any residual gap to the benchmark.

## Results

The validated distribution-only model predicted a bone antibody
biodistribution coefficient of 5.70%, ~78% of the published 7.27%
benchmark, with mass balance conserved to within 1e-11 relative error; a
minimal, single-site FcRn rescue mechanism in bone's sorting endosome
raised this by only ~0.16 percentage points (to 5.86%) once the FcRn pool
was correctly scaled to whole-body abundance, while a bounding calculation
showed that even the lower-bound assumption for bone's unmodeled cellular
volume overshot the benchmark (20.3-82.8%), ruling out missing tissue
volume as the explanation for the residual gap.

## Conclusions

Convective tissue distribution alone accounts for the majority of bone
antibody exposure predicted by the platform PBPK model, and the residual
gap to the published benchmark is attributable to FcRn-mediated kinetics
rather than tissue volume, indicating that closing it will require the
complete FcRn recycling machinery (all endosomal sub-compartments,
multi-site binding, and endogenous IgG competition) rather than further
compartmental expansion - a natural extension for future denosumab-specific
model development.
