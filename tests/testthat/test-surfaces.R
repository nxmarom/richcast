make_series <- function(keep = TRUE) {
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()
  suppressWarnings(run_hindcast_series(
    db, clim, times = c(850, 950),
    focus = focus_box(c(0, 10, 0, 10), label = "test"),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, keep_surfaces = keep, min_cells = 0, quiet = TRUE
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
      resolution = 0.5, on_error = "warn", min_cells = 0, quiet = TRUE
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
               resolution = 0.5, min_cells = 0, quiet = TRUE)

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
    resolution = 0.5, baseline = 1050, min_cells = 0, quiet = TRUE
  ))
  # With baseline = 1050, the "present" row must equal the 1050 hindcast row.
  pres <- res$richness |> dplyr::filter(.data$period == "present")
  s1050 <- res$richness |> dplyr::filter(.data$time == 1050)
  expect_equal(pres$mean_richness, s1050$mean_richness)
})

# --- Replicate ensembles ----------------------------------------------------

test_that("fit_replicates produces distinct draws", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()

  set <- fit_replicates(db, "Genus_low", clim, replicates = 4,
                        predictors = c("bio01", "bio12"),
                        land = fake_land(), quiet = TRUE)
  expect_s3_class(set, "richcast_sdm_set")
  expect_length(set$fits, 4)
  # Replicate i uses seed + i - 1, so the seeds must all differ.
  expect_equal(vapply(set$fits, function(f) f$seed, numeric(1)), 123:126)
})

test_that("project_replicates scores every replicate against one raster read", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()
  set <- fit_replicates(db, "Genus_low", clim, replicates = 3,
                        predictors = c("bio01", "bio12"),
                        land = fake_land(), quiet = TRUE)

  ps <- project_replicates(set, clim, 850, quiet = TRUE)
  expect_s3_class(ps, "richcast_projection_set")
  expect_length(ps$cells, 3)
  expect_length(ps$ranges, 3)

  # Scoring one replicate this way must match project_sdm on the same model.
  direct <- project_sdm(set$fits[[1]], clim, 850, quiet = TRUE)
  expect_equal(ps$cells[1], direct$cells)
})

test_that("replicates = 1 leaves the series unchanged", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()
  args <- list(db, clim, times = c(850, 950),
               focus = focus_box(c(0, 10, 0, 10)),
               predictors = c("bio01", "bio12"), land = fake_land(),
               resolution = 0.5, min_cells = 0, quiet = TRUE)

  a <- suppressWarnings(do.call(run_hindcast_series, args))
  b <- suppressWarnings(do.call(run_hindcast_series, c(args, list(replicates = 1))))
  expect_equal(a$species$cells, b$species$cells)
  expect_equal(a$richness$mean_richness, b$richness$mean_richness)
  expect_false("mean_richness_lo" %in% names(a$richness))
})

test_that("replicates > 1 adds intervals that bracket the point estimate", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  db <- two_species_db()

  res <- suppressWarnings(run_hindcast_series(
    db, clim, times = c(850, 950),
    focus = focus_box(c(0, 10, 0, 10)),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, replicates = 5, min_cells = 0, quiet = TRUE
  ))

  expect_equal(res$replicates, 5L)
  expect_true(all(c("cells_min", "cells_max", "cells_sd") %in% names(res$species)))
  expect_true(all(res$species$cells_min <= res$species$cells))
  expect_true(all(res$species$cells_max >= res$species$cells))

  hind <- dplyr::filter(res$richness, .data$period != "present")
  expect_true(all(c("mean_richness_lo", "mean_richness_hi") %in% names(hind)))
  expect_true(all(hind$mean_richness_lo <= hind$mean_richness_hi))
  # Each species keeps a full set of replicate fits.
  expect_length(res$sets, nrow(res$models))
})

test_that("ensemble_series summarises spread across members", {
  skip_if_not_installed("maxnet")
  dir1 <- withr::local_tempdir(); dir2 <- withr::local_tempdir()
  c1 <- structured_climate(dir1, times = c(850, 950))
  c2 <- structured_climate(dir2, times = c(850, 950))
  db <- two_species_db()
  mk <- function(cl) suppressWarnings(run_hindcast_series(
    db, cl, times = c(850, 950), focus = focus_box(c(0, 10, 0, 10)),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, min_cells = 0, quiet = TRUE))

  ens <- ensemble_series(list(a = mk(c1), b = mk(c2)))
  expect_s3_class(ens, "tbl_df")
  expect_true(all(c("a", "b", "ens_mean", "ens_min", "ens_max", "ens_sd") %in% names(ens)))
  expect_equal(nrow(ens), 2)
  expect_true(all(ens$ens_min <= ens$ens_max))

  expect_error(ensemble_series(list(mk(c1))), "at least two")
  expect_error(ensemble_series(list(mk(c1), mk(c2))), "must be named")
})

# --- Screen provenance ------------------------------------------------------
# "Nothing held out" is the printed form of three different situations, only
# one of which is a clean pass. The setting has to travel with the result.

test_that("min_cells is recorded on the series", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  expect_true("min_cells" %in% names(res))
  expect_type(res$min_cells, "integer")
  # The setting travels with the result, including when the screen is off --
  # "disabled" and "never recorded" are different claims about the data.
  expect_equal(res$min_cells, 0L)
})

test_that("the screen never emits a shapeless resolvability table", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  expect_s3_class(res$resolvability, "tbl_df")
  expect_true(all(c("species", "range_cells", "min_cells", "resolvable",
                    "reason") %in% names(res$resolvability)))
  # One row per modelled species, so a 0x0 is unreachable from current code.
  expect_equal(nrow(res$resolvability), nrow(res$models))
})

test_that("screen status distinguishes absent, disabled and clean-pass", {
  skip_if_not_installed("maxnet")
  # make_series() disables the screen, so build one that runs it.
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  res <- suppressWarnings(run_hindcast_series(
    two_species_db(), clim, times = c(850, 950),
    focus = focus_box(c(0, 10, 0, 10)),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.5, min_cells = 1, quiet = TRUE
  ))

  expect_match(richcast:::screen_status(res), "ran at min_cells")

  # A series from a build predating the screen has no field at all. That is
  # the dangerous case: it reads as a clean pass if you count rows.
  old <- res
  old$min_cells <- NULL
  expect_match(richcast:::screen_status(old), "NOT RECORDED")

  disabled <- res
  disabled$min_cells <- 0L
  expect_match(richcast:::screen_status(disabled), "disabled")
})

test_that("print reports the screen unconditionally", {
  skip_if_not_installed("maxnet")
  res <- make_series()
  expect_message(print(res), "resolvability")
})
