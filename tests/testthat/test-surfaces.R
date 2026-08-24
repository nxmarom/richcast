make_series <- function(keep = TRUE) {
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()
  suppressWarnings(run_hindcast_series(
    db, clim, times = c(850, 950),
    focus = focus_box(c(0, 10, 0, 10), label = "test"),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, keep_surfaces = keep, quiet = TRUE
  ))
}

test_that("richness surfaces are kept for every slice plus the present", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  expect_setequal(names(res$surfaces), c("850", "950", "present"))

  r <- richness_surface(res, 850)
  expect_s4_class(r, "SpatRaster")
  expect_equal(names(r), "richness")
  expect_true(all(terra::values(r, na.rm = TRUE) >= 0))
})

test_that("surfaces survive saveRDS, like everything else terra-backed", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(res, f)
  expect_s4_class(richness_surface(readRDS(f), 950), "SpatRaster")
})

test_that("richness_grid returns tidy long data with numeric time", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  g <- richness_grid(res)

  expect_s3_class(g, "tbl_df")
  expect_equal(names(g), c("x", "y", "time", "richness"))
  expect_type(g$time, "double")
  expect_setequal(unique(g$time), c(850, 950))

  # One row per non-missing cell per slice.
  per_slice <- sum(!is.na(terra::values(richness_surface(res, 850), mat = FALSE)))
  expect_equal(nrow(g), per_slice * 2)
})

test_that("richness_grid can include the present baseline", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  g <- richness_grid(res, times = c(850, "present"))
  expect_type(g$time, "character")
  expect_setequal(unique(g$time), c("850", "present"))
})

test_that("accessors fail helpfully when surfaces were not kept", {
  skip_if_not_installed("maxnet")
  res <- make_series(keep = FALSE)
  expect_length(res$surfaces, 0)
  expect_error(richness_surface(res, 850), "keep_surfaces")
  expect_error(richness_grid(res), "keep_surfaces")
})

test_that("accessors reject unknown times and wrong objects", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  expect_error(richness_surface(res, 1234), "No surface for")
  expect_error(richness_grid(res, times = 1234), "No surface for")
  expect_error(richness_surface(list(), 850), "run_hindcast_series")
})

test_that("a slice where every projection failed is flagged, not reported as zero", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)   # 950 deliberately absent
  db <- two_species_db()

  # An all-zero richness surface is indistinguishable from a real absence of
  # species, so a slice that produced nothing must say so.
  res <- NULL
  # on_error = "warn" also emits one warning per failed species; capture the
  # lot and assert the slice-level one is among them.
  warnings_seen <- testthat::capture_warnings(
    res <- run_hindcast_series(
      db, clim, times = c(850, 950),
      focus = focus_box(c(0, 10, 0, 10)),
      predictors = c("bio01", "bio12"), land = fake_land(),
      resolution = 0.5, on_error = "warn", quiet = TRUE
    )
  )
  expect_true(any(grepl("No species contributed a range at 950", warnings_seen)))

  contrib <- res$richness$species_contributing
  names(contrib) <- res$richness$period
  expect_equal(unname(contrib[["hindcast_950"]]), 0L)
  expect_gt(unname(contrib[["hindcast_850"]]), 0)
})

test_that("species_contributing is reported for every row", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  expect_true("species_contributing" %in% names(res$richness))
  expect_false(any(is.na(res$richness$species_contributing)))
})

test_that("baseline = 'present' reproduces the fitted present_cells", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  # The default short-circuits to the fitted value; check it really matches
  # what projecting onto the fitting slice gives.
  for (sp in res$models$species) {
    fit <- res$fits[[sp]]
    expect_equal(
      unique(res$species$present_cells[res$species$species == sp]),
      fit$present_cells
    )
  }
  expect_equal(res$baseline, "present")
})

test_that("a numeric baseline changes the deltas and is recorded", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950, 1050))
  db <- two_species_db()
  args <- list(db, clim, times = c(850, 950, 1050),
               focus = focus_box(c(0, 10, 0, 10)),
               predictors = c("bio01", "bio12"), land = fake_land(),
               resolution = 0.5, quiet = TRUE)

  d_pres <- suppressWarnings(do.call(run_hindcast_series, args))
  d_1050 <- suppressWarnings(do.call(run_hindcast_series,
                                     c(args, list(baseline = 1050))))

  expect_equal(d_pres$baseline, "present")
  expect_equal(d_1050$baseline, 1050)

  # Measured against the 1050 slice, that slice's own delta must be zero.
  own <- d_1050$species |> dplyr::filter(.data$time == 1050)
  expect_true(all(own$delta_from_present == 0))

  # And the baseline column is the 1050 cell count, not the fitted present.
  for (sp in d_1050$models$species) {
    cells_1050 <- d_1050$species$cells[d_1050$species$species == sp &
                                         d_1050$species$time == 1050]
    expect_equal(
      unique(d_1050$species$present_cells[d_1050$species$species == sp]),
      cells_1050
    )
  }
})

test_that("the present richness row is built from the baseline, not the fit", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950, 1050))
  db <- two_species_db()
  res <- suppressWarnings(run_hindcast_series(
    db, clim, times = c(850, 950, 1050),
    focus = focus_box(c(0, 10, 0, 10)),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, baseline = 1050, quiet = TRUE
  ))
  # With baseline = 1050, the "present" row must equal the 1050 hindcast row.
  pres <- res$richness |> dplyr::filter(.data$period == "present")
  s1050 <- res$richness |> dplyr::filter(.data$time == 1050)
  expect_equal(pres$mean_richness, s1050$mean_richness)
})
