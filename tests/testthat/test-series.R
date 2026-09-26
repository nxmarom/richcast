run_fixture_series <- function(dir, times = c(950, 850), ...) {
  clim <- structured_climate(dir, times = sort(times))
  res <- run_hindcast_series(
    two_species_db(), clim, times = times, region = region(c(0, 10, 0, 10)),
    land = fake_land(), num_trees = 50, quiet = TRUE, ...
  )
  list(res = res, clim = clim)
}

test_that("run_hindcast_series orders time chronologically", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  x <- run_fixture_series(withr::local_tempdir())
  res <- x$res
  expect_s3_class(res, "richcast_series")
  expect_equal(res$times, c(850, 950))
  expect_equal(res$richness$period, c("present", "hindcast_850", "hindcast_950"))
  expect_equal(res$richness$species_modelled, c(2L, 2L, 2L))
  expect_equal(nrow(res$species), 4)
  expect_equal(res$species$delta_from_present,
               res$species$cells - res$species$present_cells)
})

test_that("models carry AUC and Boyce for each member and the ensemble", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  res <- run_fixture_series(withr::local_tempdir())$res
  expect_true(all(c("auc_maxent", "auc_rf", "auc_ensemble",
                    "boyce_maxent", "boyce_rf", "boyce_ensemble")
                  %in% names(res$models)))
  expect_setequal(res$models$species, c("Genus_low", "Genus_high"))
})

test_that("richness surfaces count thresholded species and sum suitability", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  res <- run_fixture_series(withr::local_tempdir())$res
  s <- richness_surface(res, "present", layer = "both")
  expect_equal(names(s), c("richness", "expected"))

  # The present richness surface is exactly the sum of the two binary maps.
  manual <- Reduce(`+`, lapply(res$fits, function(f) {
    b <- terra::resample(suitability(f, binary = TRUE), s, method = "near")
    terra::ifel(is.na(b), 0, b)
  }))
  v <- terra::values(s[["richness"]], mat = FALSE)
  ok <- !is.na(v)
  expect_equal(v[ok], terra::values(manual, mat = FALSE)[ok])
  expect_lte(max(v, na.rm = TRUE), 2)

  e <- terra::values(s[["expected"]], mat = FALSE)
  expect_true(all(e[ok] >= 0 & e[ok] <= 2))

  # Region cells per species agree with the surface.
  low_present <- res$species$region_present_cells[res$species$species == "Genus_low"][1]
  expect_equal(low_present, as.integer(terra::global(
    terra::resample(suitability(res$fits$Genus_low, binary = TRUE), s, method = "near"),
    "sum", na.rm = TRUE)[1, 1]))
})

test_that("surfaces survive saveRDS and tidy up into a long grid", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  res <- run_fixture_series(withr::local_tempdir())$res
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(res, path)
  back <- readRDS(path)
  expect_equal(terra::values(richness_surface(back, 850)),
               terra::values(richness_surface(res, 850)))

  g <- richness_grid(res)
  expect_named(g, c("x", "y", "time", "richness", "expected"))
  expect_type(g$time, "double")
  expect_setequal(unique(g$time), c(850, 950))
  expect_type(richness_grid(res, c("present", 850))$time, "character")

  expect_error(richness_surface(res, 1234), "No surface")
  expect_error(richness_surface(list(), 850), "run_hindcast_series")
})

test_that("keep_surfaces = FALSE makes the accessors explain themselves", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  res <- run_fixture_series(withr::local_tempdir(), keep_surfaces = FALSE)$res
  expect_error(richness_surface(res, 850), "keep_surfaces")
})

test_that("richness_at matches the surfaces and lists the expected species", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  x <- run_fixture_series(withr::local_tempdir())
  res <- x$res
  pt <- richness_at(res, lon = c(2.25, 7.75), lat = c(2.25, 7.75),
                    time = c("present", 850), climate = x$clim)
  expect_s3_class(pt, "richcast_point")
  expect_equal(nrow(pt$richness), 4)
  expect_equal(nrow(pt$species), 8)

  for (tt in c("present", "850")) {
    s <- richness_surface(res, tt, layer = "both")
    at <- terra::extract(s, cbind(c(2.25, 7.75), c(2.25, 7.75)))
    got <- pt$richness[pt$richness$time == tt, ]
    got <- got[order(got$lon), ]
    expect_equal(got$richness, as.integer(at$richness))
    expect_equal(got$expected_richness, at$expected, tolerance = 1e-6)
  }

  # Each corner holds the species whose range sits there.
  low <- pt$species[pt$species$lon == 2.25 & pt$species$present, ]
  expect_true("Genus_low" %in% low$species)
  high <- pt$species[pt$species$lon == 7.75 & pt$species$present, ]
  expect_true("Genus_high" %in% high$species)

  # A single model works too.
  one <- richness_at(res$fits$Genus_low, 2.25, 2.25, climate = x$clim)
  expect_equal(one$species$species, "Genus_low")
  expect_error(richness_at(list(1), 0, 0, climate = x$clim), "richcast_series")
})

test_that("only species within the distance are trained", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  # Genus_low sits at 0.5-3.5, Genus_high at 6-9. A region at the far corner
  # with distance 1 reaches Genus_high only.
  res <- run_hindcast_series(
    two_species_db(), clim, times = 850, region = region(c(9.5, 10, 9.5, 10)),
    distance = 1, land = fake_land(), num_trees = 50, quiet = TRUE
  )
  expect_equal(names(res$fits), "Genus_high")
})

test_that("run_hindcast_series validates its inputs", {
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  expect_error(run_hindcast_series(db, clim, times = "850",
                                   region = region(c(0, 10, 0, 10))), "times")
  expect_error(run_hindcast_series(db, clim, times = 850,
                                   region = c(0, 10, 0, 10)), "region")
  expect_error(run_hindcast_series(db, clim, times = 850,
                                   region = region(c(100, 110, 0, 10)),
                                   quiet = TRUE), "No species")
})

test_that("a failing species is skipped rather than killing the run", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  # A species whose range lies off the climate grid cannot be fitted.
  off <- build_taxon_db(
    sf::st_sf(species = "Genus_off",
              geometry = sf::st_sfc(square(30, 30, 2), crs = 4326)),
    quiet = TRUE
  )
  db <- rbind(db, off)
  expect_warning(
    res <- run_hindcast_series(
      db, clim, times = 850, region = region(c(0, 10, 0, 10)),
      species = c("Genus_low", "Genus_off"), land = fake_land(),
      num_trees = 50, quiet = TRUE
    ),
    "Failed to fit"
  )
  expect_equal(names(res$fits), "Genus_low")
  expect_error(
    suppressWarnings(run_hindcast_series(
      db, clim, times = 850, region = region(c(0, 10, 0, 10)),
      species = "Genus_off", land = fake_land(), quiet = TRUE
    )),
    "No species could be modelled"
  )
})
