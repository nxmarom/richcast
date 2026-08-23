# End-to-end tests on synthetic data. Shared fixtures -- fake_land(),
# structured_climate(), two_species_db() -- live in helper-fixtures.R so that
# test-surfaces.R can use them too.

test_that("a model is fitted once and projected onto many slices", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950, 1050))
  db <- two_species_db()

  m <- fit_sdm(db, "Genus_low", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), quiet = TRUE)
  expect_s3_class(m, "richcast_sdm")
  expect_equal(m$species, "Genus_low")
  expect_true(m$present_cells > 0)

  # The same fitted object serves every slice; no refitting.
  p1 <- project_sdm(m, clim, 850, quiet = TRUE)
  p2 <- project_sdm(m, clim, 1050, quiet = TRUE)
  expect_s3_class(p1, "richcast_projection")
  expect_equal(p1$threshold, m$threshold)
  expect_equal(p2$threshold, m$threshold)
})

test_that("a fitted model survives saveRDS", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()

  m <- fit_sdm(db, "Genus_low", clim, predictors = c("bio01", "bio12"),
               land = fake_land(), quiet = TRUE)

  # terra objects held bare would reload as invalid external pointers, which
  # is why rasters are wrapped and vectors kept as sf.
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(m, f)
  m2 <- readRDS(f)

  expect_equal(m2$species, m$species)
  expect_s4_class(suitability(m2), "SpatRaster")
  expect_gt(as.numeric(sum(sf::st_area(sf::st_as_sf(m2$present_range)))), 0)
  expect_s3_class(project_sdm(m2, clim, 950, quiet = TRUE), "richcast_projection")
})

test_that("run_hindcast_series orders time chronologically and lags correctly", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  # 850 and 950 sort AFTER 1850 as character, which is exactly the bug that
  # mis-lagged deltas in the source pipeline. Feed them out of order too.
  times <- c(1850, 950, 850, 1050)
  clim <- structured_climate(dir, times = times)
  db <- two_species_db()

  res <- suppressWarnings(run_hindcast_series(
    db, clim, times = times,
    focus = focus_box(c(0, 10, 0, 10), label = "test"),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, quiet = TRUE
  ))
  expect_s3_class(res, "richcast_series")

  # Slices come back in chronological order, not lexical.
  hind <- dplyr::filter(res$richness, .data$period != "present")
  expect_equal(hind$time, c(850, 950, 1050, 1850))

  # The earliest slice per species has no predecessor; later ones do.
  first_rows <- res$species |>
    dplyr::group_by(.data$species) |>
    dplyr::slice_min(.data$time, n = 1) |>
    dplyr::ungroup()
  expect_true(all(is.na(first_rows$delta_from_previous)))
  expect_true(all(first_rows$time == 850))

  later <- dplyr::filter(res$species, .data$time > 850)
  expect_true(all(!is.na(later$delta_from_previous)))

  # And the lag really is the previous chronological slice.
  one <- res$species |>
    dplyr::filter(.data$species == .data$species[1]) |>
    dplyr::arrange(.data$time)
  expect_equal(
    one$delta_from_previous[-1],
    diff(one$cells)
  )
})

test_that("run_hindcast_series validates its inputs", {
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  f <- focus_box(c(0, 10, 0, 10))

  expect_error(run_hindcast_series(db, clim, times = character(0), focus = f),
               "non-empty numeric")
  expect_error(
    run_hindcast_series(db, clim, times = 850, focus = focus_box(c(150, 160, 60, 70))),
    "No species ranges intersect"
  )
})

test_that("a failing species is skipped rather than killing the run", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()

  # Slice 950 exists but 1750 does not, so every projection onto it fails.
  res <- suppressWarnings(run_hindcast_series(
    db, clim, times = c(850, 950),
    focus = focus_box(c(0, 10, 0, 10)),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, on_error = "warn", quiet = TRUE
  ))
  expect_gt(nrow(res$models), 0)
  expect_true(all(c("auc", "threshold", "present_cells") %in% names(res$models)))
})
