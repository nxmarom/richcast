# Point (occurrence) ranges, as gbif_occurrences() produces them.
#
# This path had no coverage at all, which is how three interacting defects
# survived in it: the background was erased with zero-area geometry and came
# back empty, multi-row species collapsed to their first occurrence, and the
# progress message counted features rather than points.

# build_taxon_db() does not dissolve point sources, so one row per occurrence
# is the shape fit_sdm() has to cope with.
occurrence_db <- function(n = 40, seed = 1) {
  set.seed(seed)
  pts <- sf::st_sfc(
    lapply(seq_len(n), function(i) {
      sf::st_point(c(1 + stats::runif(1) * 3, 1 + stats::runif(1) * 3))
    }),
    crs = 4326
  )
  build_taxon_db(
    sf_polygons(sf::st_sf(species = rep("Genus_points", n), geometry = pts)),
    quiet = TRUE
  )
}

test_that("a point source keeps one row per occurrence", {
  db <- occurrence_db(40)
  expect_equal(nrow(db), 40)
  expect_true(all(sf::st_geometry_type(db) == "POINT"))
})

test_that("fit_sdm uses every occurrence, not just the first", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- occurrence_db(40)

  # Silent: several rows for one species is normal for occurrence data, so it
  # must not trigger the "using the first" warning meant for polygons.
  expect_no_warning(
    f <- fit_sdm(db, "Genus_points", clim, predictors = c("bio01", "bio12"),
                 land = fake_land(), quiet = TRUE)
  )
  # Presences are the occurrences themselves, not nsample draws.
  expect_gt(f$n_presence, 1)
  expect_lte(f$n_presence, 40)
  expect_equal(f$n_background, 1000)
})

test_that("more occurrences reach the model as more presences", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)

  few <- fit_sdm(occurrence_db(10, seed = 2), "Genus_points", clim,
                 predictors = c("bio01", "bio12"), land = fake_land(),
                 quiet = TRUE)
  many <- fit_sdm(occurrence_db(400, seed = 2), "Genus_points", clim,
                  predictors = c("bio01", "bio12"), land = fake_land(),
                  quiet = TRUE)
  expect_lte(few$n_presence, 10)
  expect_gt(many$n_presence, few$n_presence)

  # For occurrences, range_cells counts distinct occupied cells, so it is
  # bounded by the record count and rises with survey effort rather than
  # describing the range on its own.
  expect_lte(few$range_cells, 10)
  expect_lte(many$range_cells, many$n_presence)
  expect_gt(many$range_cells, few$range_cells)
})

test_that("point ranges get a background instead of aborting", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- occurrence_db(40)

  # Erasing zero-area geometry from the study extent yields no features, so
  # the polygon path's guard would abort here on a range that fills nothing.
  f <- fit_sdm(db, "Genus_points", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), quiet = TRUE)
  expect_s3_class(f, "richcast_sdm")
  expect_gt(f$n_background, 0)
})

test_that("a polygon range still loses its area to the background", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)

  # The polygon branch must be untouched by the point fix: the range is still
  # erased, and a range filling its own study extent still aborts.
  db <- two_species_db()
  f <- fit_sdm(db, "Genus_low", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), quiet = TRUE)
  expect_s3_class(f, "richcast_sdm")

  expect_error(
    fit_sdm(db, "Genus_low", clim, predictors = c("bio01", "bio12"),
            land = fake_land(), extent = c(0.5, 3.5, 0.5, 3.5), quiet = TRUE),
    "No background area"
  )
})

test_that("a point range can be projected like any other", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- occurrence_db(60)

  f <- fit_sdm(db, "Genus_points", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), quiet = TRUE)
  p <- project_sdm(f, clim, 850, quiet = TRUE)
  expect_s3_class(p, "richcast_projection")
  expect_true(is.numeric(p$cells))
})
