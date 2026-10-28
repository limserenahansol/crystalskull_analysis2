Subject: Re: face poking crystal skull imaging result

Hi Mark,

Answers below, and one decision I need from you and Greg before I restart imaging.

**1. Excited vs inhibited cells per mouse (face poke; 18 recordings, 4 mice, 43,156 cells)**

| Animal (day)   | Recordings | Cells  | Excited      | Inhibited     | Population sig. |
|----------------|-----------:|-------:|-------------:|--------------:|----------------:|
| m1 (13 Aug)    | 4          | 9,705  | 896 (9.2%)   | 2,406 (24.8%) | 4/4 |
| m201 (13 Aug)  | 6          | 15,467 | 1,151 (7.4%) | 4,496 (29.1%) | 5/6 |
| m3 (13 Aug)    | 2          | 3,685  | 475 (12.9%)  | 1,113 (30.2%) | 2/2 |
| m200 (18 Aug)  | 2          | 5,681  | 553 (9.7%)   | 1,105 (19.5%) | 0/2 |
| m201 (18 Aug)  | 4          | 8,618  | 644 (7.5%)   | 2,788 (32.4%) | 3/4 |
| **All**        | **18**     | **43,156** | **3,719 (8.6%)** | **11,908 (27.6%)** | **14/18** |

(m201 was imaged on both days; each recording is one mirror, so 9 movies x 2 mirrors.) The direction is
the same in every animal: inhibition is about 3x more common than excitation (per-recording
8.8 +/- 2.4% excited vs 28.2 +/- 4.6% inhibited).

**2. Allen areas** — I cannot give named areas yet, and I do not want to guess them. The crystal-skull
FOVs are not registered to the Allen CCF; "anterior/posterior" in those maps is just the image axis after
aligning mirror2 to mirror1. What is solid is the gradient itself: excitation biased anterior, inhibition
biased posterior, in 18/18 recordings (median rho = 0.13, range 0.03-0.23). Landmark-based registration to
the Allen dorsal-cortex map is straightforward to add, and I would rather do it once, on whichever
paradigm we commit to, than build it twice.

**3. Pain vs touch vs surprise** — this dataset cannot separate them, and no reanalysis will fix that.
There is a single stimulus delivered a single way, so there is no contrast to dissociate the categories;
a poke is simultaneously touch, nociception and an unexpected event. Separating them requires a stimulus
panel interleaved within a session: innocuous brush or air puff (touch), noxious pinch or heat (pain),
an unexpected non-somatosensory event such as a tone (surprise), plus stimulus-omission catch trials. That
is a different experiment from the one we have.

**4. Drug effect** — this is the crux, and it is why I need the decision.

- Face poke: no drug arm at all. Those 18 recordings were run to test whether the paradigm is usable,
  not to test a drug.
- Melittin: we do have a drug arm, from the older data.

To make the two comparable I re-ran the old melittin data through the *same* pipeline and the same
figure set as the face poke: 12 recordings, 6 mice (8 recordings / 4 mice melittin only; 4 recordings /
2 mice melittin + morphine). Combined deck:
`combined_paradigm/CrystalSkull_facepoke_vs_melittin_vs_morphine.pptx`.

| Cohort              | Rec / mice | Cells  | Up vs baseline | Down vs baseline | Net change in population mean dF/F |
|---------------------|-----------:|-------:|---------------:|-----------------:|-----------------------------------:|
| Melittin only       | 8 / 4      | 81,834 | 27.2 +/- 3.9%  | 64.4 +/- 5.8%    | -8.5 +/- 10.1%  (down in 7/8)      |
| Melittin + morphine | 4 / 2      | 62,533 | 31.4 +/- 0.4%  | 60.6 +/- 2.3%    | -0.3 +/- 5.3%   (down in 1/4)      |

**Short answer: morphine moves the melittin response toward more activity, not less.** Melittin alone
suppresses the cortex; with morphine on board that suppression is gone. Every measure that moves, moves
the same way:

| Measure | Melittin | + morphine | p (per recording) | p (per animal) |
|---|---:|---:|---:|---:|
| fraction UP | 27.2% | **31.4%** | 0.048 | 0.13 |
| fraction DOWN | 64.4% | 60.6% | 0.21 | 0.53 |
| UP-cell amplitude | 142% | 136% of baseline | 0.68 | 0.53 |
| net population dF/F | -8.5% | **-0.3%** | 0.15 | 0.27 |

So the direction is consistent and it is the physiologically sensible one - morphine reversing a
melittin-induced cortical suppression. But I would present it as a trend, not a result, and I want to be
straight with you about why:

1. Only fraction-UP crosses p < 0.05, and only if the two mirrors of one session count as independent
   samples. They are two views of the same recording in the same animal. At the animal level
   (n = 4 vs 2 mice) nothing is significant.
2. No vehicle arm, different animals in the two cohorts, different recording dates - the comparison is
   confounded with cohort.
3. Baseline and injection are separate movies with independent EXTRACT IDs, so no cell is followed across
   the drug. Everything here is population-level; we cannot say anything about "the neurons" themselves.

One thing that could strengthen this immediately with no new mice: if 1014m01b / 1014m02w are in fact the
same animals as CS1014-1b / CS1014-2w on a later day, this becomes a paired within-animal test on data we
already have. Can someone confirm the animal identities?

**The decision I need from you and Greg.** For the SBI drug work these are two different experimental
plans and we should commit to one:

- *Melittin* — already has a working drug arm and a sustained, ongoing pain read-out, which is the
  regime an analgesic should act on. Slower, terminal-ish, one injection per session.
- *Face poke* — fast, repeatable, many trials per session, trial-locked so it gives clean single-cell
  statistics, but it is a phasic mixed stimulus and we have no drug data with it at all.

This is where we left the last conversation, so I have paused crystal-skull imaging until you and Greg
confirm which one we build the SBI plan around. I would rather hold for a few days than collect another
18 recordings on the paradigm we end up not using. Once you decide I will add the Allen registration and,
if we go with face poke, design the stimulus panel in point 3 so the pain/touch/surprise question is
answerable by design rather than after the fact.

Best,
Hansol
