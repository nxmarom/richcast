# richcast 0.0.0.9000 (development)

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

### Withdrawn: the planned resolvability screen

An earlier note here proposed omitting species whose ranges were too small to
model stably, reporting them as `not resolvable`. That was based on
measurements taken under `tss`, where 11 of 28 rodents produced a zero range on
at least one draw and median CV was 112%.

Under `p10` the same 32 species give **zero** zero-range draws, a median CV of
3.0%, and all 32 clearing any sensible size threshold. Every species that
looked unresolvable was an artefact of the rule: *Arvicola amphibius* goes from
4 cells to 84561, *Castor fiber* from 0 to 37180.

The screen is therefore not implemented. Species genuinely too small for the
grid may still exist, but none of the evidence gathered so far demonstrates
one, and a screen justified by superseded measurements would remove real data
on false grounds.

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

* Implement the resolvability screen described above, once the cutoff is
  grounded in the large-range measurements.
* Plotting helpers for richness surfaces and per-species trajectories.
* `pkgdown` site.
