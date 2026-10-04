# richcast (development version)

* **Fossil presences**: `fit_sdm()` and `run_hindcast_series()` take an
  optional `fossils` table of dated occurrences, read by the new
  `fossil_presences()` (a `slices` column of `"ka:prob;..."`, or long form
  with `ka_bp` and `weight`). Each species gets `n_fossil` draws (as many as
  its pseudo-presences by default), spread evenly over the dated units, each
  taking the climate of a slice drawn by the unit's probabilities. Matching
  background is drawn from the same slices. The study extent grows to take
  in the sites. The `p10` threshold and the hold-out metrics use the range's
  points only. Without `fossils`, fitting is unchanged.
* Fossil data are not bundled. The deer records used to develop the feature
  come from the ROCEEH Out of Africa Database (ROAD; Kandel et al. 2023,
  CC BY-SA 4.0) via the roadDB package; both are now in
  `citation("richcast")` and `?fossil_presences`.
* `vignette("middle-east")` has a new section comparing ungulate richness in
  the Levantine corridor with dated archaeological site density, using a
  circular-shift test. The aggregated site series ships as
  `inst/extdata/levant_site_series.csv` (CC BY-SA 4.0; built from ROAD, NERD
  and p3k14c, see `inst/extdata/README.md`). Its sources are in
  `citation("richcast")` and the vignette's new reference list.
* **TSS threshold switch**: `richness_at()`, `richness_in()` and
  `suitability_grid()` take `threshold = "p10"` or `"tss"` to choose the
  cutoff that sets presence. The default, `NULL`, keeps each model's fitted
  threshold (p10 unless fitted otherwise), so existing results are unchanged.
  `fit_sdm()` now stores both cutoffs on every model in `$cutoffs`, and
  accepts `threshold = "tss"` (also through `run_hindcast_series(...)`) to fit
  and build richness surfaces with the TSS cutoff. The TSS cutoff is the one
  `presence_thresholds()` already reported. Models saved by an earlier
  version can be given their cutoffs without refitting by the new
  `refresh_thresholds(x, db)`.
* **Summed suitability is deprecated.** The `expected` surface layer, the
  `mean_expected` column and `expected_richness` in `richness_at()` and
  `richness_in()` sum the models' suitabilities, which are not calibrated
  probabilities of presence (MaxEnt on the cloglog scale; a random forest
  trained on balanced presence and background draws). Their sum is therefore
  not a richness estimate. They remain for now but will be removed;
  `richness_surface(layer = "expected")` and `layer = "both"` warn, and the
  print methods no longer show them. Compare the p10 and TSS counts instead.

# richcast 0.1.0

## Revised and simplified model

* **Input** is IUCN range polygons from a folder: `iucn_folder()` reads every
  shapefile under it (e.g. `bovid_IUCN/`, `cervid_IUCN/`, `equid_IUCN/`), and
  `build_taxon_db()` accepts a folder path directly. Only polygons coded as
  extant (`PRESENCE` 1-3) and native or reintroduced (`ORIGIN` 1-2) are kept
  by default.
* **`fit_sdm()` is now an ensemble** of a random forest (`ranger`, balanced
  down-sampling) and MaxEnt (`maxnet`), averaged with equal weight. It trains
  on 100 pseudo-presences from inside the range polygon and 1000 background
  points, over a study extent of the range's bounding box plus 30% of its
  diagonal. The `p10` threshold is applied to the ensemble.
* **Metrics**: AUC and the continuous Boyce index for each member and the
  ensemble, on a 25% hold-out, in `$metrics` and in the series' `$models`.
* **Regions**: `region()` takes a preset (`"europe"`, `"asia"`,
  `"middle_east"`, `"africa"`, `"north_america"`, `"south_america"`) or a
  box, replacing `focus_box()`, `focus_global()` and `focus_polygon()`.
* **Species selection**: `run_hindcast_series()` trains only species whose
  present range lies within 10 degrees of the region (`species_near()`). The
  distance is fixed and always applied; `species =` can narrow the set but not
  add to it.
* **Species lists**: by default `run_hindcast_series()` models the species
  reported from Pleistocene zooarchaeological and palaeontological sites in
  the preset region (`zooarch_taxa()`, with evidence and sources for every
  entry); `species =` takes your own list instead. The Middle East list (14
  bovids, cervids, equids and wild boar) is compiled; the other presets need a list
  passed in until theirs are.
* **Merged taxa**: `run_hindcast_series(merge = )` treats several species
  as one identification (e.g. `Dama_sp` for *Dama dama* and
  *D. mesopotamica*). Each member is modelled on its own range and the
  outputs combined: the taxon is present wherever any member is. The Middle
  East list merges fallow deer, gazelles, wild goat and ibex, and hartebeest
  and oryx.
* **Predictors** are the eight canonical bioclim variables (`bioclim_vars`:
  bio01, bio04, bio05, bio06, bio12, bio15, bio16, bio17), the default for
  both `fit_sdm()` and `prepare_climate()`. Other variables are refused.
* **Richness** is stacked from the models' own rasters on the climate grid,
  with an `expected` layer (the sum of suitabilities) beside the thresholded
  count.
* **`richness_in()`** pools a small focus area (a `region()` box of a few
  cells): a species is present if it clears its threshold in any cell.
* **`presence_thresholds()`** gives each model's p10 threshold alongside a
  minimum-presence and a TSS-maximising threshold, computed from its
  present-day surface, for reading presence as a band; **`suitability_grid()`**
  returns cell-level suitability for every model over an area.
* **`richness_at()`** predicts richness and lists the expected species at any
  coordinate and time.

## Documentation

* New precomputed vignette, `vignette("middle-east")`: Middle Eastern
  ungulates on Beyer2020, 120 ka to present, with a focus area on
  Mediterranean Israel. `vignettes/precompute.R` re-renders it.
  Richness is shown as a band between the strict (p10) and lenient (TSS)
  counts, each divided by its own median and kernel-smoothed separately;
  each taxon's panel shows its own strict-to-lenient band.
* README rewritten around the pipeline; PDF reference manual
  (`richcast-manual.pdf`).

## Removed

* GBIF occurrences and point ranges, replicate fits (`fit_replicates()`,
  `project_replicates()`), `ensemble_series()`, `gaussian_window()`,
  subregions, the per-species focus-margin columns, the `min_cells`
  resolvability screen, `richness_stack()` and `richness_stats()`, and the
  `tss`/`mtp` thresholds. `modEvA` and `rgbif` are no longer used.
* The precomputed Tian Shan vignette, which documented the previous method.

# richcast 0.0.0.9000 (earlier development)

Extracted and generalised from the Tian Shan rodent-richness pipeline.

## Ingest

* Three range backends: `iucn_shapefile()`, `gbif_occurrences()`,
  `sf_polygons()`.
* `build_taxon_db()` joins ranges to traits and produces the object the rest
  of the package consumes; `filter_taxa()` reports what each filter removed.
* `normalise_species()` canonicalises binomials on both sides of the trait
  join, so formatting mismatches stop masquerading as missing species.
* `rodent_traits` bundled (Ecke et al. 2022, CC BY 4.0). No IUCN or climate
  data is redistributed.

## Climate

* `prepare_climate()` streams and aggregates slices one variable at a time,
  writing straight to disk and skipping what already exists, so an interrupted
  run resumes.
* `climate_dir()` and `pastclim_climate()` are interchangeable climate sources.
* `check_climate_grids()` verifies every slice shares one grid, with zero
  tolerance.
* `gaussian_window()` describes time-averaging over a window; slices missing
  from the ends of a reconstruction are dropped and the remaining weights
  renormalised rather than silently down-weighting the result.

## Models

* `fit_sdm()` and `project_sdm()` are separate: a model is fitted once against
  present-day climate and then projected onto any number of slices.
* `suitability()` unwraps the stored surface from either object.
* `run_hindcast_series()` runs the whole pipeline over a series of slices;
  `richness_stack()` and `richness_stats()` remain usable on their own.
* `focus_box()` / `focus_global()` / `focus_polygon()` define the summary
  region. `focus_global()` degrades to medium resolution when
  rnaturalearthhires is absent instead of failing.

## Bugs fixed on the way in

Four defects in the source pipeline are fixed structurally rather than
patched, so the API makes them hard to reintroduce.

* **Multi-feature species were modelled as separate species.** The pipeline
  iterated over range *features*. Where a species held several IUCN features
  it was modelled repeatedly, each run overwriting the last, and it trained
  against a single fragment while the rest of its true range went into the
  background sample -- asking the model to discriminate the species from
  itself. `build_taxon_db()` now unions by species. Verified on a global
  Rodentia download: 3090 features collapse to 2345 species, with the 248
  trait-matched species reproduced exactly.

* **Time slices were sorted as text, mis-lagging every delta.** `arrange()` on
  a character year puts 850 and 950 *after* 1850, so the earliest slice
  differenced against the latest and the second-earliest lost its predecessor.
  In the source pipeline's own output this corrupted 64 of 352 rows.
  `run_hindcast_series()` sorts numerically throughout.

* **The present-day model was refitted once per time slice.** Identical inputs
  and a fixed seed made every refit redundant, and each slice overwrote the
  previous slice's per-species diagnostics on disk. Fitting and projection are
  now separate calls.

* **Saved results contained dead pointers.** `saveRDS()` on a list of terra
  `SpatVector`s writes external pointers that reload as invalid handles, so
  the pipeline's saved polygon files raise "external pointer is not valid" on
  first use. richcast stores vectors as `sf` and rasters via `terra::wrap()`,
  and a round-trip is covered by tests.

A fifth was latent rather than active: aggregation applied both when slices
were prepared and again inside the projection step, which would have coarsened
the projection grid relative to the fitting grid had the Gaussian branch been
enabled. Aggregation now happens only in `prepare_climate()`, and
`check_climate_grids()` catches any mismatch.

## Vignettes

* `vignette("richcast")` -- the API on synthetic data; every chunk runs at
  build time.
* `vignette("tianshan")` -- a full analysis: 32 rodent species across the Tian
  Shan, 850-1850 CE. Precomputed via `vignettes/precompute.R`, because it needs
  a global IUCN download and ~170 MB of prepared climate that cannot be
  shipped or fetched at build time.

## Richness surfaces

* `run_hindcast_series()` gains `keep_surfaces` (default `TRUE`) and retains
  the richness raster for every slice. A regional mean says how many species,
  not where they are.
* `richness_surface(series, time)` returns one slice as a `SpatRaster`;
  `richness_grid(series)` returns them all in long `x`/`y`/`time`/`richness`
  form for a faceted heatmap. Surfaces are stored wrapped, so the series still
  survives `saveRDS()`.

## Per-species occupancy of the focus

* `res$species` gains `focus_cells`, `focus_present_cells`,
  `focus_delta_from_previous` and `focus_delta_from_present`: each projected
  range clipped to `focus`, on the grid the richness surface is built at.

  `cells` is counted over the species' study extent -- its range plus the
  fitting buffer -- so for a widespread taxon it is a continental number.
  `res$richness` is a focus quantity. Nothing previously reported the two on
  the same region, so a per-species trajectory plotted beside a richness curve
  silently answered a different question, and species selected by the
  `prefilter_buffer` expansion but absent from the focus had trajectories
  indistinguishable from those of species actually present. In the Tian Shan
  vignette that was 4 of 32 species with zero cells in the study region at
  every slice, and it inverted the growth ranking: *Marmota baibacina* leads
  on `cells` at +30% net growth and is flat at -3.4% on `focus_cells`.

  The two columns are counts on different grids -- `cells` on the climate
  grid, `focus_cells` on `resolution` -- so compare each against itself across
  slices, and use area to compare one against the other. See the
  two-geographies section of `?run_hindcast_series`.

* `res$species` also gains `focus_suit_q90` and `focus_suit_margin`: the
  ninetieth percentile of suitability inside `focus` before thresholding, and
  that value less the species' own threshold.

  A cell count cannot distinguish a focus inside a species' niche from one
  sitting on its cutoff, and the two behave completely differently. In the
  Tian Shan vignette *Ellobius talpinus* has a negative margin at every slice
  and its regional footprint runs 2, 1221, 260, 120, 1418 cells across
  consecutive centuries while its range-wide count barely moves -- the
  threshold flickering, not a range responding. Nine of 26 species there have
  a negative margin, holding 0.4-17% of the region above their own cutoff.

  Screening on volatility instead is not equivalent: a fold-change rule on
  `focus_cells` catches only three of those nine, missing *Castor fiber*,
  which varies less than two-fold and led the regional growth ranking while
  holding 2-5% of the region above threshold.

* `project_replicates()` returns `representative` (the index of the
  median-extent replicate) and `suit` (that replicate's suitability surface).
  `run_hindcast_series()` now takes the representative index from there
  instead of recomputing it.

## Baselines

* `run_hindcast_series()` gains `baseline`, naming the slice `delta_from_present`
  is measured against and the `present` richness row is built from. It is
  projected through the same `project_sdm()` path as the hindcast slices, so
  the two are commensurable. `"present"` (the default) uses the fitting slice
  and reproduces the previous numbers exactly.

  This exists because the comparison was silently unsound whenever the
  present-day slice came from a different climate product than the
  palaeoclimate series -- fitting on WorldClim while projecting onto CHELSA,
  which is what the source pipeline did. Same model, same threshold, only the
  present raster swapped: *Marmota baibacina*'s present-day range came out 42%
  smaller on WorldClim than on CHELSA, flipping it from above-present in 2 of
  11 centuries to above-present in all 11. The offset is species-specific in
  sign as well as size -- across six species it ranged from -42% to +177% --
  so it distorts comparisons between taxa as well as within them.

  richcast cannot tell which product a directory of GeoTIFFs came from, so it
  cannot warn automatically; `?run_hindcast_series` documents the trap and the
  two ways out.

* The Tian Shan vignette now runs CHELSA throughout, present slice included,
  rather than fitting on WorldClim. Its richness minimum moves from 1250 to
  1350 CE, and the present-day row is comparable with the hindcast rows for
  the first time. The WorldClim-fitted run is retained as the contrast.

## Also fixed

* **Signed years.** `climate_dir()` and `prepare_climate()` applied `abs()` to
  the year when building slice directory names, so a request for 1050 BCE
  resolved to the 1050 CE directory and returned the wrong climate silently;
  `prepare_climate()` would likewise have written both to one directory and
  skipped the second as already present. `sprintf("time_%04d", -1050)` gives
  `time_-1050`, which is the correct on-disk name, and slice detection now
  parses the sign. This mattered most for the deep-time case the pipeline was
  originally written for.
* Slice detection round-trips each candidate through the label template, so
  sibling directories that merely end in digits -- `present_1985`,
  `chelsa_1950` -- are no longer mistaken for time slices.
* `check_climate_grids()` skips directories the naming scheme does not
  recognise. Download tooling leaves scratch directories such as `.staging`
  beside the data, and those are not slices with a problem; they are not
  slices.

* Geometry predicates now retry in planar mode when s2 rejects the input.
  Published range maps routinely contain rings that are valid in the plane but
  self-intersecting on the sphere, and `st_make_valid()` does not reliably
  repair them, because planar and spherical validity are different properties.
  A real Red List download trips this in both `st_union()` and
  `st_intersects()`.

## Threshold rule: `p10` replaces `tss` as the default

The binarisation rule turned out to govern almost everything else about this
pipeline, and the inherited choice was the wrong one.

`tss` maximises the true skill statistic, which is referenced to the
*background*: widen the background and discrimination gets easier, so the
optimum cutoff rises and less and less clears it. Measured across eight
ungulates, widening the background buffer moved the median threshold from 0.78
to 0.86 and the median modelled range from 18 cells to 2. Three separate
background definitions were tried and every one made things worse -- the
background was never the problem, the rule was.

That last sentence needs a qualification, from an independent check on a
different pipeline. A tidysdm BART ensemble over the same region, with a
tightly scoped background, put its TSS optimum at 0.287 against p10's 0.242 --
close enough that the two series correlated at 0.987 and gave an identical
headline count, with p10 simply rescaling area by 1.4x. So `tss` is not wrong
in itself; it is *background-sensitive*, and richcast makes that bite because
its study extents are each species' full range plus a buffer, which for a
Eurasian species is continental. Where the background is narrow the two rules
can agree closely. Where it is continental, TSS optima land at 0.7-0.94 and
strip ranges to single digits.

`p10` is the safer default precisely because it is indifferent to that choice.

`p10` takes the tenth percentile of predictions at training presences. Being
referenced to the presences, it cannot be tightened by a wider background.
`mtp` (minimum training presence) is also available, and a fixed numeric
cutoff still works.

Measured on the same eight ungulates, five seeds each:

| rule | median threshold | median cells | median CV | zero-range draws | resolvable |
|---|---|---|---|---|---|
| tss | 0.775 | 18 | 19.3% | 8 | 2/8 |
| fixed 0.5 | 0.500 | 1442 | 2.8% | 0 | 6/8 |
| **p10** | 0.520 | 1365 | **2.4%** | **0** | **7/8** |
| mtp | 0.075 | 2967 | 5.2% | 0 | 7/8 |

The contrast that matters survives: the coefficient of variation of the
LGM-minus-present difference falls from 19.8% under `tss` to 5.1% under `p10`,
and LGM/present ratios stay dispersed (0.04 to 3.1 across taxa) rather than
compressing toward 1 as they do under `mtp`.

Validated against the known ranges rather than on internal statistics alone.
Rasterising each species' IUCN polygon onto the prediction grid, `p10` recovers
a modelled range 1.2 times the known range at the median -- slightly larger, as
climatic suitability should be -- occupying about a quarter of the study
extent. `tss` recovers 0.0 times it: for half the species tested it predicts
nothing at all, including *Spermophilus pygmaeus*, whose 8329-cell range it
reduces to zero while `p10` returns 8757.

`nsample` also rises from 100 to 1000: under `p10` that cuts the worst-case
CV across 32 rodents from 28.6% to 5.5%, and it is now cheap, since dropping
`modEvA::optiThresh` made the same 120-fit grid run in 4.1 minutes rather than
23.5.

### Reinstated: the resolvability screen, on new evidence

An earlier note here withdrew a planned screen for species too small to model
stably. That withdrawal was correct on its evidence and wrong about its scope.

The evidence was 32 Tian Shan rodents, where switching to `p10` removed every
apparent problem: zero zero-range draws, median CV 3.0%, *Arvicola amphibius*
going from 4 cells to 84561. What that dataset could not show is the case the
screen was for. The narrowest of those 32 ranges, *Marmota baibacina*, covers
**3161 grid cells**. No species in the set was small enough to fail.

Eight Levantine ungulates at 0.5 degrees are:

| species | range cells | 1000 presence samples give | worst-slice CV |
|---|---|---|---|
| *Dama mesopotamica* | **11** | each cell about 91 times | 115.9% |
| *Gazella gazella* | **19** | each cell about 53 times | 194.4% |
| *Capra ibex* | **79** | each cell about 13 times | 46.3% |
| *Capra aegagrus* | 747 | 544 distinct cells | 7.9% |
| *Sus scrofa* | 12549 | 966 distinct cells | 6.9% |

Presences are drawn **with replacement**, so a range covering `k` cells gives
the model `k` distinct climate vectors however large `nsample` is. *Dama* --
which carries the *highest* AUC in the assemblage, 0.987, because separating 11
points from 1000 is easy -- is fitted on 11.

`fit_sdm()` now records `range_cells`, and `run_hindcast_series()` gains
`min_cells` (default 100). Species below it are still fitted, projected and
reported: they keep every row in `$species`, `$models`, `$ranges` and `$fits`,
and are held out of the richness surfaces only, listed in a new
`$resolvability` table with their range size and the reason. `min_cells = 0`
restores the previous behaviour.

The cutoff is measured, not picked. The ungulates leave a gap between 79 cells
(unstable) and 747 (stable), so 23 more species were measured to fill it --
rodents inside the same extent, same grid, five seeds each, 31 species from 11
to 12549 cells. All 8 below 100 cells vary by more than 25% across seeds at
their worst slice; 100 is the largest cutoff at which that holds without
exception, and at 200 it fails. All 32 Tian Shan rodents clear it 31-fold, so
the constraint that no rodent be excluded is met.

Three findings that changed the design along the way:

* **Replicate spread is the wrong diagnostic, and inverts.** Clipping one
  donor range to compact blocks of fixed size, an 8-cell range returned
  `2, 2, 2, 2, 2` cells across five seeds -- perfect agreement, worthless
  answer. Sampling with replacement from 8 cells returns the same 8 cells
  every time, so the interval narrows precisely where the estimate is least
  trustworthy. The previous advice to screen on `cells_sd` is withdrawn.
* **Study extent and suitable fraction move with the fitting buffer.** Forcing
  background rings of 1 to 8 degrees on *Gazella gazella* moved its extent
  14-fold and its suitable fraction from 19.8% to 3.6% -- on one unchanged
  species, whose trajectory stayed equally unstable throughout. `range_cells`
  does not move, because the buffer does not change the range.
* **`p10` is less background-insensitive than claimed above.** The same sweep
  moved *Gazella*'s cutoff from 0.573 to 0.217. The rule reads off the
  presences, but the model generating those predictions is fitted against the
  background, so the cutoff moves anyway. The section above overstates this;
  what survives is that `p10` does not *tighten* toward emptiness as `tss`
  does.

The screen fixes resolution, not transferability. Above the cutoff, 8 of 23
species still vary by more than 25% at their worst slice, and range size no
longer predicts it -- the 200-400 cell band is worse than the 100-200 band.
That residue concentrates in the deep slices (median across-seed CV 7.2% at
2 ka against 22.8% at 60 ka), consistent with projection beyond the fitted
climate. Clearing `min_cells` says the grid can resolve the species. It does
not say the trajectory is sound.

### Also on the way in

* The small-range guard in `fit_sdm()` was dead code. `nrow(pres_df) < nsample`
  can never be true for a polygon range, because `spatSample(replace = TRUE)`
  returns `nsample` rows for an 11-cell range and a 12000-cell one alike -- so
  a check that read as a safeguard had never fired. It now warns on
  `range_cells`, which is the quantity that varies.
* The Tian Shan vignette's diagnostics prose still described TSS-optimised
  thresholds and `nsample = 100` against output showing `p10` and 1000, and
  discussed four zero-range species that no longer exist under `p10`. Rewritten
  against what the chunk actually prints.
* `range_cells` means something slightly different for occurrence data than for
  polygons: it counts cells holding a record, so it reflects survey effort as
  well as range size. Sampling inside *Capra aegagrus*' range, 30 records report
  30 cells against the polygon's 747, while the narrow-ranged taxa converge on
  their polygon count within a few hundred records. The count is honest -- a
  model given 30 points sees 30 climate vectors -- but an exclusion there should
  be read as "too few distinct records", not "small range".
  `?richcast-resolvability` has the numbers.

## Products

* `prepare_climate()` now defaults `dataset_present` to `dataset_past`, and
  `agg_present` to `agg_past`, so one product is used throughout unless two are
  deliberately requested. The previous defaults paired a WorldClim present with
  CHELSA slices, which is precisely the mixed pipeline that made
  `delta_from_present` unsound.
* `prepare_climate()` writes `richcast_manifest.csv` recording the dataset,
  time and aggregation behind each slice, and `check_climate_products()` reads
  it and warns when the present-day slice came from a different product than
  the time slices. Directories assembled by hand carry no manifest and are
  reported as unverifiable rather than passing silently.

## Occurrence ranges could never be fitted

`gbif_occurrences()` is exported and documented, and `fit_sdm()` has a branch
for point geometry, but nothing in the test suite ever fitted one. Three
defects had accumulated there, and any one of them was fatal:

* **The background came back empty.** `terra::erase(study, range)` is how the
  background is carved out, and erasing zero-area point geometry returns zero
  features rather than the untouched extent. Every occurrence fit aborted on
  "No background area left after erasing the range" -- a diagnosis exactly
  backwards, since a point range fills nothing. Points now take the study
  extent as their background, which is the usual presence-background
  convention: the background describes availability, and a cell holding an
  occurrence was available too.
* **The model would have been fitted to a single occurrence.**
  `build_taxon_db()` deliberately skips the dissolve for point sources, so one
  row per record is the intended shape; `db_row()` then collapsed multi-row
  species to `hit[1]`. The two halves contradicted each other, and the
  warning's advice -- rebuild with `dissolve = TRUE` -- pointed at a branch
  that never runs for points. Occurrences are now combined across rows.
* **The progress message counted features, not points.** `nrow(row)` is 1 for
  a combined geometry however many records it holds, so it reported "Using 1
  occurrence point" for any dataset.

The polygon path is unchanged: ranges are still erased from their background,
and a range filling its own study extent still aborts as before.
`tests/testthat/test-point-ranges.R` covers both.

## To do

* A transferability diagnostic, for the instability the resolvability screen
  leaves behind: how far outside its fitted climate a projection has strayed,
  per species per slice.
* Plotting helpers for richness surfaces and per-species trajectories.
* `pkgdown` site.
