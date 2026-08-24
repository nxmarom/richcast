# Handoff: a resolvability criterion for richcast

Scratch note for a fresh session. Not part of the package build
(`.Rbuildignore`d); delete once the work lands.

## The question

Some species cannot support a hindcast trajectory because their range is too
small relative to the climate grid, whatever the threshold rule. richcast
should say so rather than returning numbers that look like measurements.

## What prompted it

*Gazella gazella* in the Levant, Beyer2020 at 0.5 degrees:

| species | range (deg) | study extent (cells) | suitable | % |
|---|---|---|---|---|
| Gazella_gazella | 2.3 x 5.8 | **81** | 16 | 19.8 |
| Dama_mesopotamica | 18.3 x 8.0 | 1190 | 71 | 6.0 |
| Capra_ibex | 18.1 x 5.7 | 830 | 97 | 11.7 |
| Sus_scrofa | 151.4 x 70.4 | 30155 | 14521 | 48.2 |

Across 30 slices *Gazella* never once exceeded its present-day range, sitting
4-5 SD below it including at 2 ka BP. That reads as a striking biological
result and is not one: with 16 present-day cells and hindcast values of 0-11,
every delta is a single-digit integer and the z-score divides by an SD
computed over those same single digits. One cell is 6% of its range.

Ruled out first, so they need not be revisited:

* **Not a climate discontinuity.** Beyer2020's present slice is continuous with
  the palaeo series over the Levant: present-to-2ka steps are the same
  magnitude as 2ka-to-4ka steps across all eight bioclims.
* **Not an in-sample advantage for the baseline.** Most taxa have many slices
  *above* present -- *Capra ibex* 21 of 30, *Dama mesopotamica* 20 of 30. If
  fitting on the present slice structurally favoured it, none would.

## Why the earlier attempt was withdrawn, and why that was right

`NEWS.md` records a planned screen, then its withdrawal. The original evidence
was TSS-derived: 11 of 28 rodents gave a zero range on at least one draw,
median CV 112%. Switching the default to `p10` dissolved all of it -- zero
zero-draws, median CV 3.0%, 32 of 32 resolvable, *Arvicola amphibius* going
from 4 cells to 84561.

That withdrawal stands. But the rodent evidence covered Eurasian species with
study extents in the tens of thousands of cells; it says nothing about an
81-cell extent. The screen was dropped on evidence that did not cover the case
that motivates it now.

## The criterion should be extent-relative, not absolute

An absolute cell count is what was proposed first and what the rodent data made
look unnecessary. It is also wrong in principle: 100 cells means something
different in an 81-cell extent than in a 30000-cell one.

Candidates worth testing, in rough order of appeal:

1. **Study extent in cells.** Below some floor, no trajectory is supportable.
   *Gazella* at 81 fails; *Capra ibex* at 830 is marginal.
2. **Suitable cells as a fraction of extent**, with a floor on the absolute
   count. Catches both "too few cells" and "too small a fraction to be stable".
3. **Replicate CV**, which richcast can now measure directly via
   `run_hindcast_series(replicates = N)`. The most honest criterion -- it
   measures instability rather than proxying for it -- but costs N fits per
   species, so it cannot be the default screen.

Whatever is chosen must be validated on **both** datasets before implementing:
the 32 Tian Shan rodents (large extents, should nearly all pass) and the 8
Levantine ungulates (mixed, *Gazella* should fail). A criterion that fails
species in the rodent set is too aggressive.

## A confound to settle first

*Gazella*'s buffer hit the 1-degree minimum (`buffer_range = c(1, 8)`), so its
background is a thin ring around the range rather than a broad contrast. That
may be depressing its estimate independently of range size. Test before
attributing everything to resolution: refit with a wider minimum buffer and see
whether the 16-cell estimate moves.

## Design constraints, carried forward

* Excluded species go in a separate table with their realised range size and
  the reason. Silent shrinkage of the assemblage is the failure mode to avoid.
* Exclusions must not alter `res$species`; they belong in richness and in the
  summary, so the underlying numbers stay inspectable.
* The threshold must be measured, not picked. That is what went wrong the
  first time.

## Where things stand

richcast is at the commit that made `p10` and `nsample = 1000` the defaults,
added the product manifest, and regenerated the Tian Shan vignette. All tests
pass and `R CMD check` is clean. Nothing about the screen is implemented.

Useful data already on disk, under
`/tmp/claude-1000/.../scratchpad/` (session-scoped, may be gone):
`nsample_test.rds` (rodents, TSS), `nsample_test_p10.rds` (rodents, p10),
`adarim_threshold.rds` (ungulates, four threshold rules),
`adarim_buffer.rds` (ungulates, three background definitions),
`db.rds` (built rodent taxon database).

Regenerate rather than trust these if the scratch directory has been cleared.
