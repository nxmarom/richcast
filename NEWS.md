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

## To do

* Plotting helpers for richness surfaces and per-species trajectories.
* `pkgdown` site.
