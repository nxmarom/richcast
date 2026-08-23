# Deep-time hindcasting means BCE years, and BCE years mean signed slice
# labels. Taking abs() of the year would send a request for 1050 BCE to the
# 1050 CE directory and return the wrong climate silently.

test_that("slice labels keep their sign", {
  expect_equal(sprintf("time_%04d", -1050), "time_-1050")
  expect_equal(sprintf("time_%04d", -50), "time_-050")
  expect_equal(sprintf("time_%04d", 850), "time_0850")
})

test_that("parse_slice_times recovers CE and BCE years", {
  dirs <- c("time_0850", "time_1850", "time_-050", "time_-1050")
  expect_setequal(
    richcast:::parse_slice_times(dirs, "time_%04d"),
    c(850, 1850, -50, -1050)
  )
})

test_that("parse_slice_times ignores directories that merely end in digits", {
  # present_1985 and chelsa_1950 sit alongside the slices and must not be
  # mistaken for them; the round-trip through sprintf is what rules them out.
  dirs <- c("time_0850", "present_1985", "chelsa_1950", "notes", "figure")
  expect_equal(richcast:::parse_slice_times(dirs, "time_%04d"), 850)
})

test_that("climate_dir addresses BCE slices without collision", {
  dir <- withr::local_tempdir()

  write_slice <- function(label, value) {
    d <- file.path(dir, label)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                     ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- value
    terra::writeRaster(r, file.path(d, "bio01.tif"), overwrite = TRUE)
  }

  write_slice("present_1985", 0)
  write_slice("time_0850", 850)     # 850 CE
  write_slice("time_-850", -850)    # 850 BCE, distinct directory

  clim <- climate_dir(dir)
  expect_setequal(clim$times, c(-850, 850))

  ce  <- richcast:::climate_at(clim, 850, "bio01")
  bce <- richcast:::climate_at(clim, -850, "bio01")

  # The two must not resolve to the same raster.
  expect_equal(unique(terra::values(ce, mat = FALSE)), 850)
  expect_equal(unique(terra::values(bce, mat = FALSE)), -850)
})

test_that("a non-default present directory can be used as the fitting slice", {
  # Fitting on a palaeoclimate model's own present-day slice, rather than a
  # different product, avoids a cross-dataset step between fit and projection.
  dir <- withr::local_tempdir()
  for (lbl in c("chelsa_1950", "time_0850")) {
    d <- file.path(dir, lbl)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                     ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- if (lbl == "chelsa_1950") 1 else 2
    terra::writeRaster(r, file.path(d, "bio01.tif"), overwrite = TRUE)
  }

  clim <- climate_dir(dir, present = "chelsa_1950")
  expect_equal(clim$times, 850)
  expect_equal(
    unique(terra::values(richcast:::climate_at(clim, "present", "bio01"), mat = FALSE)),
    1
  )
})

test_that("bp_to_ce and ce_to_bp round-trip on the 1950 convention", {
  # Palaeoclimate reconstructions quote BP from 1950, so their slices land on
  # 1950, 950, -50, -1050 rather than round thousands. Passing a raw BP year
  # where CE is expected is an off-by-1950 error that returns real data from
  # the wrong millennium.
  expect_equal(bp_to_ce(0), 1950)
  expect_equal(bp_to_ce(-20000), -18050)
  expect_equal(ce_to_bp(1950), 0)
  expect_equal(ce_to_bp(-18050), -20000)
  expect_equal(ce_to_bp(bp_to_ce(-12345)), -12345)

  # Beyer2020's millennial steps, as richcast wants them.
  expect_equal(bp_to_ce(seq(0, -3000, by = -1000)), c(1950, 950, -50, -1050))
})

test_that("deep-time slice labels survive the round trip", {
  # A five-digit negative year still names a directory that parses back.
  ce <- bp_to_ce(-20000)
  lbl <- sprintf("time_%04d", ce)
  expect_equal(lbl, "time_-18050")
  expect_equal(richcast:::parse_slice_times(lbl, "time_%04d"), ce)
})
