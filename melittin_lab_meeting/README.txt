Melittin lab-meeting figure pack (face-poke visual analog)
==============================================================

This is NOT poke-aligned. Baseline and injection are independently extracted.
Do not treat injection neuron #k as the same cell as baseline neuron #k.

Injection window: last N seconds of BOTH arms, N = min(baseline, injection) duration.
Per injection neuron: amp = mean dF/F - mean of all baseline neurons mean dF/F.
Up / down: two-sided t-test of 1 s dF/F bins vs that BL session mean, p<0.05.
Population p: ranksum of 1 s population mean dF/F, BL vs INJ.
mirror2 is flipped so P is left and A is right, same as face-poke maps.
