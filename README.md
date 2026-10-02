# richcast

<!-- badges: start -->
[![R-CMD-check](https://github.com/nxmarom/richcast/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nxmarom/richcast/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Hindcast species distributions and assemblage richness through time.**

`richcast` models where a set of species could have lived under past climates
and adds them up into species richness. For a region, it takes each species'
present-day IUCN range, fits a species distribution model to it, projects that
model onto a palaeoclimate reconstruction slice by slice, and stacks the
results into richness surfaces. You can then ask what the richness was, and
which species were expected, anywhere in the region at any slice.

The reference manual is [`richcast-manual.pdf`](richcast-manual.pdf); a worked
analysis is in `vignette("middle-east")`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("nxmarom/richcast", build_vignettes = TRUE)
```

Without `build_vignettes = TRUE` the package installs and works, but the
vignettes are silently absent. `rnaturalearthhires`, needed only for a
high-resolution land mask, lives on r-universe:

```r
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

## The pipeline

```
IUCN polygons ──► build_taxon_db() ──┐
                                     ├─► run_hindcast_series() ──► richness, maps,
pastclim slices ─► prepare_climate() ┘        per taxon:              metrics
                                          fit_sdm() once,
region() + species list ─────────────►    project_sdm() per slice  ──► richness_at()
                                                                       richness_in()
```

### 1. Ranges

`iucn_folder()` reads every IUCN Red List shapefile under a folder, so
downloads for several groups (`bovid_IUCN/`, `cervid_IUCN/`, `sus_IUCN/`, ...)
unpacked side by side are read as one source. Only polygons coded as extant
(`PRESENCE` 1-3) and native or reintroduced (`ORIGIN` 1-2) are kept, since an
extinct or introduced patch is not part of the niche a model should learn.
`build_taxon_db()` then dissolves each species' polygons into a single range.

```r
db <- build_taxon_db(iucn_folder("UngulatePolygons"))
```

### 2. Climate

`prepare_climate()` streams the slices from
[`pastclim`](https://evolecolgroup.github.io/pastclim/) once, one variable at a
time, and writes them to disk; an interrupted run resumes. Models use the
eight canonical bioclimatic variables in `bioclim_vars` (bio01, bio04, bio05,
bio06, bio12, bio15, bio16, bio17), and nothing else.

```r
steps <- pastclim::get_time_bp_steps(dataset = "Beyer2020")
clim <- prepare_climate("beyer", times = bp_to_ce(steps[steps < 0]),
                        extent = c(-180, 180, -60, 90),
                        dataset_past = "Beyer2020", agg_past = 1)
```

Years are CE throughout; `bp_to_ce()` and `ce_to_bp()` convert.

### 3. Region and species

`region()` takes a preset (`"europe"`, `"asia"`, `"middle_east"`, `"africa"`,
`"north_america"`, `"south_america"`) or a box `c(xmin, xmax, ymin, ymax)`.

The species default to the region's list of taxa reported from **Pleistocene
zooarchaeological and palaeontological sites**, `zooarch_taxa()`, with the
evidence and source for each entry (compiled so far for the Middle East). Pass
`species =` to use your own list instead.

Whatever the list, only species whose present range lies **within 10 degrees**
of the region are trained (`species_near()`). The rule is fixed: naming a
far-away species cannot bring it in, so no llamas turn up in Britain.

Species that the record does not tell apart form one **merged taxon**
(`merge =`, e.g. `Dama_sp` for *Dama dama* and *D. mesopotamica*). Each member
is modelled on its own range and the outputs are combined: the taxon is present
wherever any member is. Modelling the union of their ranges instead lets the
larger range swamp the smaller.

### 4. The model

For each species, `fit_sdm()`:

* sets a **study extent** of the range's bounding box plus 30% of its
  diagonal;
* draws **100 pseudo-presences** inside the range polygon and **1000
  background points** from the rest of the extent;
* fits a **random forest** (`ranger`, balanced down-sampling) and **MaxEnt**
  (`maxnet`) and averages them with equal weight;
* sets presence at the **p10 threshold**, the tenth percentile of ensemble
  suitability at the training presences;
* scores **AUC and the continuous Boyce index** for each member and the
  ensemble on a 25% hold-out.

Each model is fitted once, on present-day climate, and `project_sdm()` applies
it to every slice within its own study extent.

### 5. Richness

```r
res <- run_hindcast_series(db, clim, times = bp_to_ce(steps[steps < 0]),
                           region = region("middle_east"))
res$richness   # per slice: mean, median and max richness; mean expected
res$models     # per species: AUC and Boyce for RF, MaxEnt and the ensemble
res$species    # per taxon and slice: occupied cells, and change from today
```

Each slice gives two surfaces on the climate grid: **`richness`**, the number
of taxa above their threshold in a cell, and **`expected`**, the sum of their
suitabilities, which needs no threshold. `richness_surface()` and
`richness_grid()` return them for mapping.

### 6. Points, focus areas and thresholds

* `richness_at(res, lon, lat, time, climate)` gives richness, expected
  richness and the ranked list of expected taxa at a coordinate.
* `richness_in(res, area, time, climate)` does the same for a small **focus
  area** of a few cells: a taxon is present if it clears its threshold in any
  of them. `suitability_grid()` returns the cell-level values underneath.
* `presence_thresholds(res, db)` gives each model's p10 threshold together
  with a TSS-maximising and a minimum-presence threshold, so presence can be
  read as a band from strict to lenient rather than a single line.
  `vignette("middle-east")` plots richness this way: the strict and lenient
  counts are each divided by their own median, smoothed separately, and drawn
  as the band between them.

## Where the data comes from

`richcast` bundles **no range data and no climate data**.

* **Ranges** are IUCN Red List polygons you download yourself. The Red List
  Terms of Use (v3, section 4) prohibit redistributing them "in whole, or in
  part... including within Derivative Works", so the package ships none, and
  the vignette prints none. Cite them with the Red List version and download
  date, in the format the IUCN spatial-data metadata gives: "IUCN <year>. The
  IUCN Red List of Threatened Species. <version>. https://www.iucnredlist.org.
  Downloaded on <date>." (see `?iucn_folder`).
* **Climate** comes through `pastclim` (Leonardi et al. 2023), which handles
  the downloads. Cite the reconstruction too, as pastclim asks: Beyer2020
  (Beyer et al. 2020, used in the Middle East vignette), CHELSA-TraCE21k
  (Karger et al. 2023) or WorldClim 2.1 (Fick & Hijmans 2017).
* **Land outlines**, which mask predictions and background points to land,
  are Natural Earth via `rnaturalearth`. Natural Earth is public domain and
  asks for no credit, but suggests "Made with Natural Earth" if you give one.
* **Traits** for rodents ship as `rodent_traits` (Ecke et al. 2022, CC BY 4.0).
* **Fossil occurrences**, for `fit_sdm(fossils = )`, are yours to supply;
  none ship. The dated deer records behind the fossil-calibrated Middle East
  runs come from the ROCEEH Out of Africa Database (ROAD; Kandel et al.
  2023), retrieved with the [roadDB](https://doi.org/10.32614/CRAN.package.roadDB)
  package. ROAD content is CC BY-SA 4.0: cite Kandel et al. (2023) and share
  derived tables under the same licence.
* **Levantine site series**, the one table that ships
  (`inst/extdata/levant_site_series.csv`, used in the Middle East vignette), is
  an aggregate of dated site-phases built from ROAD, NERD (Palmisano et al.
  2022) and p3k14c (Bird et al. 2022), calibrated with IntCal20 (Reimer et al.
  2020) and corrected after Surovell et al. (2009). It is CC BY-SA 4.0, not
  MIT; see `inst/extdata/README.md`.

## Vignettes

* `vignette("richcast")`: the API on synthetic data you can run yourself.
* `vignette("middle-east")`: the full pipeline on Middle Eastern ungulates,
  Beyer2020, 120 ka to the present, with a focus area on Mediterranean Israel
  and a comparison with archaeological site density in the Levantine
  corridor.
  It is precomputed, because its inputs cannot be redistributed;
  `vignettes/precompute.R` re-renders it.

## Citation

Please cite the package and every data source you used:

* the IUCN Red List version and download date of your range polygons;
* pastclim (Leonardi et al. 2023) and the palaeoclimate reconstruction:
  Beyer et al. (2020), Karger et al. (2023) or Fick & Hijmans (2017);
* the sources behind any `zooarch_taxa()` list you relied on;
* Ecke et al. (2022) if you used `rodent_traits`;
* Kandel et al. (2023) and roadDB if you trained on fossil occurrences from
  ROAD;
* Kandel et al. (2023), Palmisano et al. (2022), Bird et al. (2022), Reimer
  et al. (2020) and Surovell et al. (2009) if you used the Levantine site
  series.

`citation("richcast")` lists them with full references.

## How this was built

The original code is by Nimrod Marom; later development was done with Claude
(Anthropic). See [AI_COLLABORATION.md](AI_COLLABORATION.md).

## License

MIT for the code. `rodent_traits` is CC BY 4.0 (Ecke et al. 2022). Data you
supply remains under its own terms.
