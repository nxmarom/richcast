# A tiny on-disk climate source, so these tests need no downloads.
make_fake_climate <- function(dir, times, vars = c("bio01", "bio12"),
                              present = "present_1985", value_offset = 0) {
  for (lbl in c(present, sprintf("time_%04d", times))) {
    d <- file.path(dir, lbl)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    for (v in vars) {
      r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10,
                       ymin = 0, ymax = 10, crs = "EPSG:4326")
      terra::values(r) <- seq_len(100) + value_offset
      names(r) <- v
      terra::writeRaster(r, file.path(d, paste0(v, ".tif")), overwrite = TRUE)
    }
  }
  climate_dir(dir, present = present)
}

test_that("climate_dir discovers slices from directory names", {
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = c(850, 950, 1850))
  expect_s3_class(clim, "richcast_climate")
  expect_equal(clim$times, c(850, 950, 1850))
})

test_that("climate_at returns requested variables and errors informatively", {
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = c(850, 950))

  r <- richcast:::climate_at(clim, 850, c("bio01", "bio12"))
  expect_s4_class(r, "SpatRaster")
  expect_setequal(names(r), c("bio01", "bio12"))

  expect_s4_class(richcast:::climate_at(clim, "present", "bio01"), "SpatRaster")
  expect_error(richcast:::climate_at(clim, 1234, "bio01"), "No climate rasters")
  expect_error(richcast:::climate_at(clim, 850, "bio99"), "missing")
})

test_that("gaussian_window normalises weights and centres on zero", {
  w <- gaussian_window()
  expect_equal(sum(w$weights), 1)
  expect_equal(w$offsets, c(-200, -100, 0, 100, 200))
  expect_equal(which.max(w$weights), 3L)

  w2 <- gaussian_window(step = 50, weights = c(1, 2, 1))
  expect_equal(w2$offsets, c(-50, 0, 50))
  expect_equal(w2$weights, c(0.25, 0.5, 0.25))

  expect_error(gaussian_window(weights = c(1, 1)), "odd-length")
  expect_error(gaussian_window(weights = c(1, -1, 1)), "non-negative")
})

test_that("a window averages slices and renormalises around missing ones", {
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = c(850, 950, 1050))

  # All three slices hold identical values, so any correctly normalised
  # weighting must return those same values -- this catches the classic
  # weights-do-not-sum-to-one bug.
  avg <- richcast:::climate_for_projection(
    clim, 950, "bio01", NULL,
    window = gaussian_window(step = 100, weights = c(1, 1, 1)),
    quiet = TRUE
  )
  expect_equal(
    terra::values(avg, mat = FALSE),
    as.numeric(seq_len(100)),
    tolerance = 1e-6
  )

  # A window running off the end of the reconstruction should renormalise the
  # surviving weights, not silently shrink the values.
  edge <- richcast:::climate_for_projection(
    clim, 1050, "bio01", NULL,
    window = gaussian_window(step = 100, weights = c(1, 1, 1)),
    quiet = TRUE
  )
  expect_equal(
    terra::values(edge, mat = FALSE),
    as.numeric(seq_len(100)),
    tolerance = 1e-6
  )
})

test_that("check_climate_grids catches a mismatched slice", {
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = c(850, 950))
  expect_true(check_climate_grids(clim, quiet = TRUE))

  # Rewrite one slice at a different resolution.
  bad <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 10,
                     ymin = 0, ymax = 10, crs = "EPSG:4326")
  terra::values(bad) <- seq_len(25)
  terra::writeRaster(bad, file.path(dir, "time_0950", "bio01.tif"),
                     overwrite = TRUE)

  expect_error(check_climate_grids(clim, quiet = TRUE), "not consistent")
})

test_that("prepare_climate rejects a file where a directory is required", {
  # Pointing `path` at the source NetCDF is an easy confusion, and left
  # unchecked it surfaces much later as "[writeRaster] cannot write file".
  f <- withr::local_tempfile(fileext = ".nc")
  file.create(f)
  expect_error(
    prepare_climate(path = f, vars = "bio01", times = 850,
                    extent = c(0, 10, 0, 10)),
    "existing file"
  )
  expect_error(
    prepare_climate(path = f, vars = "bio01", times = 850,
                    extent = c(0, 10, 0, 10)),
    "set_data_path"
  )
})
