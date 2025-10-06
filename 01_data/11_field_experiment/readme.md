## Premise

These files pertain to a field experiment at Juleboda in Sweden, begun in Autumn 
2025.

The idea is to sow pools of known numbers of F10 seeds of the intercross lines generated
in the common-garden phenotyping experiment at two sites by the beach, and two 
sites just inland in Autumn, then return later and sample what has survived.
There are two plots per site.

Full details are given on eLabJournal:
https://vbc.elabjournal.com/members/experiments/browser/#view=experiment&nodeID=263840&page=0&userID=0&status=0&column=created&order=DESC&search=

- `seed_tube_checklist.R` processes the seeds I received from LabDeers
- `seed_tube_checklist.tsv` gives a sample sheet of the tubes processed, and
  which pools they were processed into:
  - **label**: Position label from the common-garden experiment
  - **genotype**: Expected genotype as printed on the tubes.
  - **id**: Updated lines labels
  - **i**: The order of each tube in the seed catalogue from LabDeers
  - **pool**: Integer from 1 to 8 giving the group into which seeds were pooled
  - **checked**: For some tubes I was able to go through empty tubes after
    sowing and double check I had sown what I though. Note that I hadn't kept
    all the tubes, so there are gaps.