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
