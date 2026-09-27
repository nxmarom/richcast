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

# --- Product provenance -----------------------------------------------------

write_manifest_fixture <- function(dir, present_ds, past_ds) {
  utils::write.csv(
    data.frame(slice = c("present", "time_0850", "time_0950"),
               time_ce = c(1950, 850, 950),
               dataset = c(present_ds, past_ds, past_ds),
               aggregation = 20,
               variables = "bio01"),
    file.path(dir, "richcast_manifest.csv"), row.names = FALSE)
}

test_that("a mixed-product pipeline is detected and explained", {
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = c(850, 950), present = "present")
  write_manifest_fixture(dir, "WorldClim_2.1_5m", "CHELSA_trace21k_1.0_0.5m_vsi")

  expect_warning(check_climate_products(clim), "different products")
  expect_warning(check_climate_products(clim), "does not cancel")
  expect_false(suppressWarnings(check_climate_products(clim)))
})

test_that("a single-product pipeline passes", {
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = c(850, 950), present = "present")
  write_manifest_fixture(dir, "CHELSA_trace21k_1.0_0.5m_vsi",
                         "CHELSA_trace21k_1.0_0.5m_vsi")
  expect_message(check_climate_products(clim), "come from")
  expect_true(check_climate_products(clim, quiet = TRUE))
})

test_that("slices assembled by hand are reported as unverifiable, not as passing", {
  # No manifest means richcast genuinely cannot tell; saying so is better than
  # implying a check happened.
  dir <- withr::local_tempdir()
  clim <- make_fake_climate(dir, times = 850, present = "present")
  expect_message(check_climate_products(clim), "cannot verify")
  expect_true(check_climate_products(clim, quiet = TRUE))
})
