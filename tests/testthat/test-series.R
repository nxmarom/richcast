run_fixture_series <- function(dir, times = c(950, 850), ...) {
  clim <- structured_climate(dir, times = sort(times))
  res <- run_hindcast_series(
    two_species_db(), clim, times = times, region = region(c(0, 10, 0, 10)),
    species = c("Genus_low", "Genus_high"), land = fake_land(),
    num_trees = 50, quiet = TRUE, ...
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

test_that("the 10-degree rule cannot be overridden by naming a species", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  llama <- build_taxon_db(
    sf::st_sf(species = "Lama_glama",
              geometry = sf::st_sfc(square(40, 40, 2), crs = 4326)),
    quiet = TRUE
  )
  db <- rbind(two_species_db(), llama)
  expect_false("distance" %in% names(formals(run_hindcast_series)))
  expect_warning(
    res <- run_hindcast_series(
      db, clim, times = 850, region = region(c(0, 10, 0, 10)),
      species = c("Genus_low", "Lama_glama"), land = fake_land(),
      num_trees = 50, quiet = TRUE
    ),
    "Lama_glama"
  )
  expect_equal(names(res$fits), "Genus_low")
  expect_error(
    suppressWarnings(run_hindcast_series(
      db, clim, times = 850, region = region(c(0, 10, 0, 10)),
      species = "Lama_glama", quiet = TRUE
    )),
    "No species to model"
  )
})

test_that("run_hindcast_series validates its inputs", {
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  expect_error(run_hindcast_series(db, clim, times = "850",
                                   region = region(c(0, 10, 0, 10))), "times")
  expect_error(run_hindcast_series(db, clim, times = 850,
                                   region = c(0, 10, 0, 10)), "region")
  expect_error(suppressWarnings(run_hindcast_series(
    db, clim, times = 850, region = region(c(100, 110, 0, 10)),
    species = "Genus_low", quiet = TRUE)), "No species")
  # A custom box has no bundled zooarchaeological list, so one must be given.
  expect_error(run_hindcast_series(db, clim, times = 850,
                                   region = region(c(0, 10, 0, 10)),
                                   quiet = TRUE), "zooarchaeological")
})

test_that("a failing species is skipped rather than killing the run", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  # A species near the region but off the climate grid cannot be fitted.
  off <- build_taxon_db(
    sf::st_sf(species = "Genus_off",
              geometry = sf::st_sfc(square(12, 12, 2), crs = 4326)),
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

test_that("a preset region defaults to its zooarchaeological species list", {
  me <- zooarch_taxa("middle_east")
  expect_true(all(c("species", "evidence", "source") %in% names(me)))
  expect_true(all(c("Dama_sp", "Gazella_gazella", "Capra_aegagrus",
                    "Equus_hemionus") %in% me$species))
  expect_false(any(duplicated(me$species)))
  expect_true(all(nzchar(me$source)))
  expect_equal(zooarch_taxa(region("middle_east")), me)
  expect_equal(nrow(zooarch_taxa("europe")), 0)
  expect_equal(richcast:::default_taxa(region("middle_east")), me$species)
  expect_error(richcast:::default_taxa(region("europe")), "Pass your own")
})

test_that("richness_in pools a block of cells and agrees with richness_at", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  x <- run_fixture_series(withr::local_tempdir())
  # A 2 x 2 block of 0.5-degree cells around (2.5, 2.5).
  block <- region(c(2, 3, 2, 3), label = "block")
  fa <- richness_in(x$res, block, time = c("present", 850), climate = x$clim)
  expect_s3_class(fa, "richcast_area")
  expect_equal(fa$richness$cells, c(4L, 4L))

  centres <- expand.grid(lon = c(2.25, 2.75), lat = c(2.25, 2.75))
  pt <- richness_at(x$res, centres$lon, centres$lat, time = c("present", 850),
                    climate = x$clim)
  per_sp <- pt$species |>
    dplyr::group_by(time, species) |>
    dplyr::summarise(n = sum(present), present = any(present),
                     suit = if (all(is.na(suitability))) NA_real_ else
                       max(suitability, na.rm = TRUE),
                     .groups = "drop")
  got <- dplyr::inner_join(fa$species, per_sp, by = c("time", "species"))
  expect_equal(nrow(got), 4)
  expect_equal(got$suitability, got$suit, tolerance = 1e-8)
  expect_equal(got$present.x, got$present.y)
  expect_equal(got$cells_present, got$n)
  expect_equal(fa$species$present,
               dplyr::coalesce(fa$species$suitability > fa$species$threshold, FALSE))
  expect_error(richness_in(x$res, c(2, 3, 2, 3), climate = x$clim), "region")
})

test_that("merged taxa are modelled per member and counted once", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  res <- run_hindcast_series(
    two_species_db(), clim, times = 850, region = region(c(0, 10, 0, 10)),
    species = "Genus_sp", merge = list(Genus_sp = c("Genus_low", "Genus_high")),
    land = fake_land(), num_trees = 50, quiet = TRUE
  )
  # Each member keeps its own model...
  expect_setequal(names(res$fits), c("Genus_low", "Genus_high"))
  expect_equal(res$taxa, list(Genus_sp = c("Genus_low", "Genus_high")))
  expect_equal(unique(res$models$taxon), "Genus_sp")
  # ...but the taxon counts once, present wherever either member is.
  s <- richness_surface(res, "present")
  expect_lte(max(terra::values(s), na.rm = TRUE), 1)
  either <- Reduce(`|`, lapply(res$fits, function(f) {
    b <- terra::resample(suitability(f, binary = TRUE), s, method = "near")
    terra::ifel(is.na(b), 0, b) > 0
  }))
  v <- terra::values(s, mat = FALSE)
  ok <- !is.na(v)
  expect_equal(v[ok], as.numeric(terra::values(either, mat = FALSE)[ok]))
  expect_equal(unique(res$species$species), "Genus_sp")

  pt <- richness_at(res, c(2.25, 7.75), c(2.25, 7.75), climate = clim)
  expect_equal(unique(pt$species$species), "Genus_sp")
  expect_equal(pt$richness$richness, c(1L, 1L))
  expect_equal(pt$species$present, pt$species$suitability > pt$species$threshold)
})
