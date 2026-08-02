DRAFT v1 - PAGE abstract, denosumab-mpbpk. Author-drafted sentences, 
tightened collaboratively. To be reviewed against PAGE word limit 
and formatting before submission. Numbers verified in 
abstract_materials.md.

## Objectives
We developed a physiologically-based pharmacokinetic model of 
denosumab tissue distribution to predict bone tissue exposure, 
because plasma exposure does not reflect the local bone exposure that 
drives denosumab's pharmacological effect.

## Methods
We reproduced the Shah & Betts 15-tissue platform antibody PBPK model 
in rxode2, applied it to denosumab, and calculated the bone antibody 
biodistribution coefficient (ABC) at steady state to compare against 
the published benchmark (Shah & Betts 2013; 7.27%).

## Results
The model predicted a bone ABC of 5.70%, reproducing ~78% of the 
benchmark by convective distribution alone; the residual gap was 
traced to FcRn recycling kinetics rather than tissue volume, the 
volume hypothesis having been excluded by a bounding analysis.

## Conclusions
This work demonstrates that platform antibody PBPK models can be 
reproduced in open-source rxode2 to evaluate tissue-level exposure. 
For target-tissue drugs like denosumab, plasma kinetics alone cannot 
predict local distribution, and the dominant mechanism - here, FcRn 
recycling rather than tissue volume - must be established by explicit 
testing rather than assumed.
