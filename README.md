# richcast

<!-- badges: start -->
<!-- badges: end -->

**Hindcast species distributions and assemblage richness through time.**

`richcast` projects species distribution models onto palaeoclimate
reconstructions across a series of time slices, then summarises the resulting
assemblage richness for a region of interest. It is taxon-agnostic: rodents,
ungulates, carnivorans and anything else with range data and a climate niche
go through the same functions.

## Installation

```r
# install.packages("pak")
pak::pak("nxmarom/richcast")
```

Or with devtools:

```r
devtools::install_github("nxmarom/richcast")
```

## Where the data comes from

`richcast` bundles **no range data and no climate data**, by design.

* **Ranges** are supplied by you, through one of three backends:
  `iucn_shapefile()` for a Red List spatial download, `gbif_occurrences()` for
  point records fetched live, or `sf_polygons()` for any `sf` object you
  already hold.
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

If you want a pipeline that runs with no manual downloads at all, use the GBIF
backend instead.

## Quick start

```r
library(richcast)

# 1. Assemble a taxon database: your ranges + a trait table
db <- build_taxon_db(
  ranges = iucn_shapefile("~/iucn_rodentia/data_0.shp"),
  traits = rodent_traits
)
#> i Reading 'data_0.shp'
#> i Dissolving 3090 features into 2345 species ranges.
#> i Trait join: 400/2345 species matched, 1945 with NA traits.

# 2. Narrow to the taxa you care about
db <- filter_taxa(db, !is.na(S_index))
#> i Filter kept 248/2345 species (2097 removed).

# 3. Prepare climate slices once (slow; resumable)
clim <- prepare_climate(
  path   = "climate/eurasia",
  vars   = c("bio01", "bio04", "bio05", "bio06",
             "bio12", "bio15", "bio16", "bio17"),
  times  = seq(850, 1850, by = 100),
  extent = c(-15, 180, 10, 82)
)
#> v 12 slices share one grid (0.1667 x 0.1667 deg).

# 4. Fit, project and summarise
res <- run_hindcast_series(
  db, clim,
  times      = seq(850, 1850, by = 100),
  focus      = focus_box(c(68, 87, 39, 46), label = "Tian Shan"),
  subregions = list(karadja = focus_box(c(74.38, 74.99, 42.58, 43.03)))
)

res$richness   # one row per slice: mean, median, max, variance
res$species    # per species per slice, with deltas
res$models     # AUC, threshold, sample sizes
```

Each species is fitted **once** and projected onto every slice, so the cost
scales with species, not with species x slices. Individual pieces —
`fit_sdm()`, `project_sdm()`, `richness_stack()`, `richness_stats()` — work on
their own if you want a different loop.

Time-averaging across chronological uncertainty is one argument:

```r
res <- run_hindcast_series(..., window = gaussian_window())
```

## Vignettes

* `vignette("richcast")` — the API, on synthetic data you can run yourself.
* `vignette("tianshan")` — a worked analysis: 32 rodent species across the
  Tian Shan, 850–1850 CE, with plague-reservoir trajectories against the
  fourteenth-century pandemic.

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
actually used — the IUCN Red List version, GBIF download DOI, the relevant
`pastclim` climate reconstruction, and Ecke et al. (2022) if you used
`rodent_traits`. `citation("richcast")` lists them.

## License

MIT for the code. Bundled `rodent_traits` is CC BY 4.0 (Ecke et al. 2022).
Data you supply remains under its own terms.
