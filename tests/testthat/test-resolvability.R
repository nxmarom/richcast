# The resolvability screen: species too small for the grid are reported and
# held out of richness, but never quietly dropped from the record.

# One wide range and one narrow one on the same grid, so a single min_cells
# separates them.
mixed_db <- function() {
  build_taxon_db(
    sf_polygons(sf::st_sf(
      species = c("Genus_wide", "Genus_narrow"),
      geometry = sf::st_sfc(square(0.5, 0.5, 4), square(7, 7, 0.6), crs = 4326)
    )),
    quiet = TRUE
  )
}

test_that("range_cells counts the grid cells the range covers", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- mixed_db()

  wide <- fit_sdm(db, "Genus_wide", clim, predictors = c("bio01", "bio12"),
                  land = fake_land(), quiet = TRUE)
  narrow <- fit_sdm(db, "Genus_narrow", clim, predictors = c("bio01", "bio12"),
                    land = fake_land(), quiet = TRUE)

  # Against the fixture's own grid, independently of how fit_sdm got there.
  present <- terra::rast(file.path(dir, "present_1985", "bio01.tif"))
  expected <- function(sp) {
    v <- terra::vect(sf::st_sf(geometry = sf::st_geometry(db[db$species == sp, ])))
    as.integer(terra::global(!is.na(terra::mask(present, v)), "sum",
                             na.rm = TRUE)[1, 1])
  }
  expect_equal(wide$range_cells, expected("Genus_wide"))
  expect_equal(narrow$range_cells, expected("Genus_narrow"))
  expect_gt(wide$range_cells, narrow$range_cells)
})

test_that("range_cells does not move with nsample", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- mixed_db()

  # Sampling with replacement recovers distinct cells only partially, so a
  # count taken from the sample would drift with nsample. This one must not.
  a <- fit_sdm(db, "Genus_wide", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), nsample = 50, quiet = TRUE)
  b <- fit_sdm(db, "Genus_wide", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), nsample = 2000, quiet = TRUE)
  expect_equal(a$range_cells, b$range_cells)
})

test_that("unresolvable species leave richness but stay in the record", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- mixed_db()
  args <- list(db, clim, times = c(850, 950), focus = focus_box(c(0, 10, 0, 10)),
               predictors = c("bio01", "bio12"), land = fake_land(),
               resolution = 0.5, quiet = TRUE)

  full <- suppressWarnings(do.call(run_hindcast_series,
                                   c(args, list(min_cells = 0))))
  cut <- max(full$models$range_cells[full$models$species == "Genus_narrow"]) + 1

  screened <- suppressWarnings(do.call(run_hindcast_series,
                                       c(args, list(min_cells = cut))))

  # Held out of richness ...
  expect_true(all(screened$richness$species_contributing <
                    full$richness$species_contributing))
  # ... but still fitted, projected and reported.
  expect_setequal(screened$species$species, full$species$species)
  expect_setequal(screened$models$species, full$models$species)
  expect_true("Genus_narrow" %in% names(screened$fits))
  expect_equal(
    screened$species$cells[screened$species$species == "Genus_narrow"],
    full$species$cells[full$species$species == "Genus_narrow"]
  )

  res <- screened$resolvability
  expect_setequal(res$species, c("Genus_wide", "Genus_narrow"))
  expect_false(res$resolvable[res$species == "Genus_narrow"])
  expect_true(res$resolvable[res$species == "Genus_wide"])
  expect_match(res$reason[res$species == "Genus_narrow"], "below min_cells")
  expect_true(is.na(res$reason[res$species == "Genus_wide"]))
  expect_false(screened$models$resolvable[screened$models$species == "Genus_narrow"])
})

test_that("min_cells = 0 reproduces the unscreened result exactly", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()
  args <- list(db, clim, times = c(850, 950), focus = focus_box(c(0, 10, 0, 10)),
               predictors = c("bio01", "bio12"), land = fake_land(),
               resolution = 0.5, quiet = TRUE)

  a <- suppressWarnings(do.call(run_hindcast_series, c(args, list(min_cells = 0))))
  b <- suppressWarnings(do.call(run_hindcast_series, c(args, list(min_cells = 1))))
  expect_equal(a$richness$mean_richness, b$richness$mean_richness)
  expect_true(all(a$resolvability$resolvable))
})

test_that("a screen nothing clears is an error, not an empty richness surface", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()

  expect_error(
    suppressWarnings(run_hindcast_series(
      db, clim, times = 850, focus = focus_box(c(0, 10, 0, 10)),
      predictors = c("bio01", "bio12"), land = fake_land(),
      resolution = 0.5, min_cells = 1e6, quiet = TRUE
    )),
    "No species clears"
  )
})

test_that("the screen still applies when replicates are requested", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- mixed_db()
  args <- list(db, clim, times = c(850, 950), focus = focus_box(c(0, 10, 0, 10)),
               predictors = c("bio01", "bio12"), land = fake_land(),
               resolution = 0.5, replicates = 3, quiet = TRUE)

  full <- suppressWarnings(do.call(run_hindcast_series,
                                   c(args, list(min_cells = 0))))
  cut <- full$models$range_cells[full$models$species == "Genus_narrow"] + 1
  res <- suppressWarnings(do.call(run_hindcast_series,
                                  c(args, list(min_cells = cut))))

  expect_false(res$resolvability$resolvable[
    res$resolvability$species == "Genus_narrow"])
  # Replicate richness intervals must exclude the same species the point
  # estimate does, or the interval describes a different assemblage.
  expect_true(all(res$richness$species_contributing[-1] <
                    full$richness$species_contributing[-1]))
  expect_true(all(is.finite(
    res$richness$mean_richness_lo[res$richness$period != "present"])))
})

test_that("min_cells is validated", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  expect_error(
    run_hindcast_series(db, clim, times = 850, focus = focus_box(c(0, 10, 0, 10)),
                        predictors = c("bio01", "bio12"), land = fake_land(),
                        min_cells = -1, quiet = TRUE),
    "min_cells"
  )
})
