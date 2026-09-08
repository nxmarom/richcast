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
    resolution = 0.5, min_cells = 0, quiet = TRUE
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

test_that("focus_cells reports the focus, cells reports the study extent", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  times <- c(850, 950, 1050)
  clim <- structured_climate(dir, times = times)
  db <- two_species_db()

  # A focus around Genus_low only. Genus_high sits at 6-9 on both axes, so it
  # is selected by the 8-degree prefilter and modelled, but has no business
  # appearing in this region's assemblage.
  res <- suppressWarnings(run_hindcast_series(
    db, clim, times = times,
    focus = focus_box(c(0, 4, 0, 4), label = "low corner"),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.25, min_cells = 0, quiet = TRUE
  ))

  expect_true(all(c("focus_cells", "focus_present_cells",
                    "focus_delta_from_previous", "focus_delta_from_present")
                  %in% names(res$species)))
  expect_type(res$species$focus_cells, "integer")
  expect_true(all(res$species$focus_cells >= 0))

  # Both species are modelled; only one is in the region.
  expect_setequal(unique(res$species$species), c("Genus_low", "Genus_high"))
  low <- dplyr::filter(res$species, .data$species == "Genus_low")
  expect_true(all(low$focus_cells > 0))

  # The focus lag follows the same chronological rule as the study-extent one.
  expect_true(all(is.na(low$focus_delta_from_previous[low$time == 850])))
  ordered <- dplyr::arrange(low, .data$time)
  expect_equal(ordered$focus_delta_from_previous[-1], diff(ordered$focus_cells))
  expect_equal(ordered$focus_delta_from_present,
               ordered$focus_cells - ordered$focus_present_cells)

  # focus_cells is commensurable with the richness surface; cells is not,
  # being counted over each species' own study extent on the climate grid.
  per_slice <- res$species |>
    dplyr::group_by(.data$time) |>
    dplyr::summarise(total = sum(.data$focus_cells), .groups = "drop") |>
    dplyr::arrange(.data$time)
  from_richness <- res$richness |>
    dplyr::filter(.data$period != "present") |>
    dplyr::arrange(.data$time)
  expect_equal(per_slice$total,
               from_richness$mean_richness * from_richness$cells)
})

test_that("focus_suit_margin separates a focus inside the niche from one on its edge", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  times <- c(850, 950, 1050)
  clim <- structured_climate(dir, times = times)
  db <- two_species_db()

  # The climate is a plane increasing with x and y, and the two species sit at
  # opposite ends of it. A focus over Genus_low's own corner is inside its
  # niche and outside Genus_high's.
  res <- suppressWarnings(run_hindcast_series(
    db, clim, times = times,
    focus = focus_box(c(0, 3, 0, 3), label = "low corner"),
    predictors = c("bio01", "bio12"), land = fake_land(),
    resolution = 0.25, min_cells = 0, quiet = TRUE
  ))

  expect_true(all(c("focus_suit_q90", "focus_suit_margin") %in% names(res$species)))
  q <- res$species$focus_suit_q90
  expect_true(all(is.na(q) | (q >= 0 & q <= 1)))

  # The margin is the q90 less that species' own threshold, per species. The
  # models table carries names on its vapply-built columns, so strip them
  # before comparing against the unnamed column.
  chk <- dplyr::left_join(res$species, res$models[, c("species", "threshold")],
                          by = "species")
  expect_equal(chk$focus_suit_margin,
               unname(chk$focus_suit_q90 - chk$threshold))

  by_sp <- res$species |>
    dplyr::group_by(.data$species) |>
    dplyr::summarise(margin = min(.data$focus_suit_margin), .groups = "drop")

  # Positive for the species whose corner of the gradient this is.
  expect_gt(by_sp$margin[by_sp$species == "Genus_low"], 0)

  # Genus_high is fitted at the far end and its study extent -- range plus a
  # buffer of 15%, floored at one degree -- never reaches this focus, so it is
  # projected nowhere near it. That is NA, not a low margin: the distinction
  # matters, because "not modelled here" and "modelled here and unsuitable"
  # are different claims and only the second is evidence of anything.
  expect_true(is.na(by_sp$margin[by_sp$species == "Genus_high"]))
})

test_that("focus_suit_q90 summarises the focus and reports no overlap as NA", {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10,
                   ymin = 0, ymax = 10, crs = "EPSG:4326")

  # Constant surface: every quantile is that constant, whatever the focus.
  terra::values(r) <- 0.8
  expect_equal(focus_suit_q90(r, focus_box(c(1, 4, 1, 4))), 0.8)

  # A tenth of the focus at 0.9 and the rest at 0.1 puts the 90th percentile
  # at the top group -- the point of using q90 rather than the maximum, which
  # a single cell would set, or the mean, which the unsuitable bulk would drag.
  terra::values(r) <- ifelse(seq_len(terra::ncell(r)) %% 10 == 0, 0.9, 0.1)
  expect_gt(focus_suit_q90(r, focus_box(c(0, 10, 0, 10))), 0.1)

  # A focus the projection does not reach at all: terra aborts on the crop, and
  # the answer is "not measured here" rather than a number.
  expect_true(is.na(focus_suit_q90(r, focus_box(c(50, 55, 50, 55)))))
})

test_that("a species that never enters the focus is visible as zero", {
  # This is the bug the column exists for: without it such a species has a
  # perfectly ordinary-looking trajectory and no way to tell it contributes
  # nothing to the richness it is plotted beside.
  focus <- focus_box(c(0, 4, 0, 4))
  n <- focus_cells_each(
    list(here = sf::st_sfc(square(1, 1, 1), crs = 4326),
         elsewhere = sf::st_sfc(square(20, 20, 1), crs = 4326)),
    focus, resolution = 0.5
  )
  expect_gt(n[["here"]], 0L)
  expect_equal(n[["elsewhere"]], 0L)
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
    resolution = 0.5, on_error = "warn", min_cells = 0, quiet = TRUE
  ))
  expect_gt(nrow(res$models), 0)
  expect_true(all(c("auc", "threshold", "present_cells") %in% names(res$models)))
})

test_that("presence-referenced thresholds are available and ordered sensibly", {
  skip_if_not_installed("maxnet")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()

  fits <- lapply(c("tss", "p10", "mtp"), function(rule)
    fit_sdm(db, "Genus_low", clim, predictors = c("bio01", "bio12"),
            land = fake_land(), threshold = rule, quiet = TRUE))
  names(fits) <- c("tss", "p10", "mtp")

  expect_true(all(vapply(fits, function(f) f$threshold > 0 && f$threshold < 1,
                         logical(1))))
  # mtp admits everything p10 does, so it can never be the stricter cutoff.
  expect_lte(fits$mtp$threshold, fits$p10$threshold)
  # A more permissive cutoff cannot yield a smaller range.
  expect_gte(fits$mtp$present_cells, fits$p10$present_cells)

  expect_equal(fits$p10$threshold_rule, "p10")
  expect_error(
    fit_sdm(db, "Genus_low", clim, predictors = c("bio01", "bio12"),
            land = fake_land(), threshold = "nonsense", quiet = TRUE),
    "must be"
  )
})
