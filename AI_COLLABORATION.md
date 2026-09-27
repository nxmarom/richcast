# AI collaboration

## Who did what

**Nimrod Marom** wrote the original code: the richness counting and the
MaxEnt hindcasting of richness surfaces from which richcast grew.

**Everything after that was developed with Claude (Anthropic)**, working in
the local repository and on GitHub. Two Claude sessions took part: one
checked that the existing package still built and ran, and one carried out
the revision released as 0.1.0. The revision covered:

* the random forest + MaxEnt ensemble, its evaluation (AUC, continuous Boyce)
  and the p10, TSS and minimum-presence thresholds;
* IUCN folder input, preset regions, the 10-degree species rule, merged taxa,
  and point and focus-area queries;
* the tests, the vignettes, the README, the reference manual, and the Middle
  East analysis.

## How the work was divided

Nimrod set the specification and made the scientific decisions. He chose the
model design and defaults, which species to include and merge, the use of
presence bands rather than a single threshold, the focus area, and the
figures. He also reviewed the results as they came in.

Claude wrote the code, tests and documentation, ran the analyses, and raised
problems as they appeared. For example, merging species ranges before
modelling erased fallow deer from the Levant, which is why merged taxa are
now combined after modelling.

## What to check

* **The Pleistocene species list** in `inst/extdata/zooarch_taxa.csv` was
  compiled by Claude from published sources and then reviewed and edited by
  Nimrod. Check the cited sources before relying on any single entry.
* **Correctness** rests on the test suite and `R CMD check`, and on the
  Middle East run being read by a specialist. Neither guarantees that the
  ecological assumptions hold. The main assumptions are niche conservatism and
  present-day IUCN ranges as the training data.

Commits made with Claude carry a `Co-Authored-By` trailer.
