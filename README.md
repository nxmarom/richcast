# richcast

<!-- badges: start -->
[![R-CMD-check](https://github.com/nxmarom/richcast/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nxmarom/richcast/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<!-- While the repository is private the badge renders only for signed-in
     users with access; it will show publicly if the repo is ever opened up. -->

**Hindcast species distributions and assemblage richness through time.**

`richcast` fits a species distribution model for every species near a region,
projects each onto palaeoclimate reconstructions across a series of time
slices, and stacks the results into richness surfaces for that region. It is taxon-agnostic: rodents,
ungulates, carnivorans and anything else with range data and a climate niche
go through the same functions.

## Installation

```r
# install.packages("remotes")
remotes::install_github("nxmarom/richcast", build_vignettes = TRUE)
```

`build_vignettes = TRUE` is worth the extra minute. Without it the package
installs correctly and nothing errors, but `vignette("richcast")` finds
nothing — the vignettes are simply absent, silently.

`pak::pak("nxmarom/richcast")` and `devtools::install_github()` also work and
have the same default.

Some dependencies are not on CRAN. `rnaturalearthhires` lives on r-universe and
is only needed for a high-resolution land mask, which degrades to medium
resolution with a warning when it is absent:

```r
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

## Where the data comes from

`richcast` bundles **no range data and no climate data**, by design.

* **Ranges** are IUCN Red List polygons you download yourself. Unpack the
  downloads side by side in one folder and point `iucn_folder()` at it;
  `sf_polygons()` takes any polygon `sf` object you already hold.
* **Climate** is retrieved through [`pastclim`](https://evolecolgroup.github.io/pastclim/),
  which handles WorldClim and CHELSA-TraCE21k downloads and caching.
* **Traits** for rodents ship with the package as `rodent_traits`, from
  Ecke et al. (2022), redistributable under CC BY 4.0.

### A note on IUCN Red List data

IUCN Red List range polygons **cannot** be redistributed. The Red List Terms
of Use (v3, section 4) prohibit redistribution "in whole, or in part... alone
or combined with other data, including within Derivative Works". That covers
any file this package might ship, so `iucn_shapefile()` points at a download
you make yourself, under your own acceptance of those terms. Cite the version
you used.

## Quick start

```r
library(richcast)

# 1. Ranges: every IUCN shapefile under a folder, dissolved per species
db <- build_taxon_db(iucn_folder("UngulatePolygons"))

# 2. Climate slices, prepared once (slow; resumable)
clim <- prepare_climate(
  path   = "climate/beyer",
  vars   = c("bio01", "bio04", "bio05", "bio06",
             "bio12", "bio15", "bio16", "bio17"),
  times  = bp_to_ce(seq(-2000, -20000, by = -2000)),
  extent = c(-30, 80, -40, 75),
  dataset_past = "Beyer2020", agg_past = 1
)

# 3. Fit, project and stack, for a preset region or your own box
res <- run_hindcast_series(
  db, clim,
  times  = bp_to_ce(seq(-2000, -20000, by = -2000)),
  region = region("middle_east")        # or region(c(xmin, xmax, ymin, ymax))
)

res$richness   # one row per slice: mean/median/max richness, mean expected
res$models     # AUC and Boyce for the RF, MaxEnt and the ensemble
res$species    # suitable cells per species per slice

# 4. Richness, and the expected species, at a coordinate
richness_at(res, lon = 35.5, lat = 33, time = c("present", -4050),
            climate = clim)
```

## The model

For each species:

* **Study extent**: the range polygon's bounding box, widened by 30% of its
  diagonal.
* **Training points**: 100 pseudo-presences sampled inside the range polygon,
  1000 background points from the rest of the extent.
* **Ensemble**: a random forest (`ranger`, balanced down-sampling) and MaxEnt
  (`maxnet`), averaged with equal weight.
* **Threshold**: `p10`, the tenth percentile of ensemble suitability at the
  training presences.
* **Evaluation**: AUC and the continuous Boyce index for each member and the
  ensemble, on a 25% hold-out.

When hindcasting a region, only species whose present-day range lies within 10
degrees of it are trained (`species_near()`). Each is fitted once and projected
onto every slice. Preset regions are `"europe"`, `"asia"`, `"middle_east"`,
`"africa"`, `"north_america"` and `"south_america"` (`region_presets` has their
boxes).

Every slice yields two surfaces: `richness`, the number of species above their
threshold, and `expected`, the sum of their suitabilities.

## Vignettes

* `vignette("richcast")` — the API, on synthetic data you can run yourself.

## Why the ranges get dissolved

Range databases store one species across many features — IUCN splits by
subspecies, seasonality and disjunct patches, so a global rodent download
holds ~3100 features for ~2350 species. If each feature is treated as a
species, the model trains on a single fragment while the *rest of that
species' true range* goes into the background sample: the model is asked to
discriminate the species from itself.

`build_taxon_db()` unions features by species by default. Pass
`dissolve = FALSE` if per-feature modelling is genuinely what you want.

## Citation

If you use `richcast`, please cite the package plus the data sources you
actually used — the IUCN Red List version, the relevant
`pastclim` climate reconstruction, and Ecke et al. (2022) if you used
`rodent_traits`. `citation("richcast")` lists them.

## License

MIT for the code. Bundled `rodent_traits` is CC BY 4.0 (Ecke et al. 2022).
Data you supply remains under its own terms.
