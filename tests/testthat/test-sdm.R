test_that("AUC is the Mann-Whitney probability, ties counting half", {
  expect_equal(richcast:::auc_score(c(0.9, 0.8), c(0.1, 0.2)), 1)
  expect_equal(richcast:::auc_score(c(0.1, 0.2), c(0.9, 0.8)), 0)
  expect_equal(richcast:::auc_score(c(0.5, 0.5), c(0.5, 0.5)), 0.5)
  expect_equal(richcast:::auc_score(c(0.9, 0.3), c(0.5, 0.1)), 0.75)
  expect_true(is.na(richcast:::auc_score(numeric(0), 0.5)))
})

test_that("Boyce is high when presences follow suitability and low when not", {
  set.seed(1)
  fit <- runif(2000)
  good <- sample(fit, 300, prob = fit^3)
  flat <- sample(fit, 300)
  bad <- sample(fit, 300, prob = (1 - fit)^3)
  expect_gt(richcast:::boyce_index(good, fit), 0.9)
  expect_lt(abs(richcast:::boyce_index(flat, fit)), 0.6)
  expect_lt(richcast:::boyce_index(bad, fit), -0.9)
  expect_true(is.na(richcast:::boyce_index(0.5, rep(0.5, 10))))
})

test_that("p10 is the tenth percentile of presence predictions", {
  obs <- c(rep(1, 10), rep(0, 5))
  pred <- c(seq(0.1, 1, by = 0.1), rep(0, 5))
  expect_equal(richcast:::resolve_threshold("p10", obs, pred),
               unname(stats::quantile(seq(0.1, 1, by = 0.1), 0.1)))
  expect_equal(richcast:::resolve_threshold(0.3, obs, pred), 0.3)
  expect_error(richcast:::resolve_threshold(1.5, obs, pred), "between 0 and 1")
  expect_error(richcast:::resolve_threshold("mtp", obs, pred), "tss")
  expect_error(richcast:::check_threshold(c("p10", "tss")), "tss")
})

test_that("the study extent is the bbox plus 30% of its diagonal", {
  g <- sf::st_sfc(square(10, 10, 3), square(13, 14, 1), crs = 4326)
  e <- as.vector(richcast:::study_extent(g, 0.3))
  buf <- 0.3 * sqrt(4^2 + 5^2)
  expect_equal(unname(e), c(10 - buf, 14 + buf, 10 - buf, 15 + buf))

  polar <- sf::st_sfc(square(0, 80, 10), crs = 4326)
  expect_lte(as.vector(richcast:::study_extent(polar, 0.3))[4], 90)
})

test_that("fit_sdm builds an RF + MaxEnt ensemble with member and ensemble metrics", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()

  f <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50,
               quiet = TRUE)
  expect_s3_class(f, "richcast_sdm")
  expect_named(f$members, c("maxent", "rf"))
  expect_equal(f$metrics$model, c("maxent", "rf", "ensemble"))
  expect_true(all(f$metrics$auc > 0.5 & f$metrics$auc <= 1))
  expect_true(all(f$metrics$boyce >= -1 & f$metrics$boyce <= 1, na.rm = TRUE))
  expect_equal(f$threshold_rule, "p10")

  # The ensemble is the mean of its members.
  d <- as.data.frame(setNames(lapply(bioclim_vars, function(v) c(5, 20)),
                               bioclim_vars))
  p <- richcast:::predict_members(f$members, d)
  expect_equal(p[, "ensemble"], (p[, "maxent"] + p[, "rf"]) / 2)

  # The 3 x 3 degree range holds 36 cells at 0.5 degrees, fewer than the 100
  # pseudo-presences asked for, so every cell is used once.
  expect_equal(f$range_cells, 36L)
  expect_equal(f$n_presence, 36L)

  s <- suitability(f)
  expect_s4_class(s, "SpatRaster")
  expect_equal(f$present_cells,
               as.integer(terra::global(s > f$threshold, "sum", na.rm = TRUE)[1, 1]))
  expect_equal(sort(unique(stats::na.omit(terra::values(suitability(f, binary = TRUE),
                                                          mat = FALSE)))),
               c(0, 1))
})

test_that("pseudo-presences and background default to 100 and 1000", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  # A fine grid, so the range has more cells than the sample asks for.
  r <- terra::rast(nrows = 200, ncols = 200, xmin = 0, xmax = 10, ymin = 0,
                   ymax = 10, crs = "EPSG:4326")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  for (lbl in c("present_1985", "time_0850")) {
    d <- file.path(dir, lbl)
    dir.create(d)
    b1 <- r; terra::values(b1) <- xy[, 1] * 3 + xy[, 2]; names(b1) <- "bio01"
    b2 <- r; terra::values(b2) <- xy[, 2] * 2 - xy[, 1]; names(b2) <- "bio12"
    terra::writeRaster(b1, file.path(d, "bio01.tif"))
    terra::writeRaster(b2, file.path(d, "bio12.tif"))
  }
  clim <- climate_dir(dir)
  f <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               predictors = c("bio01", "bio12"), num_trees = 50, quiet = TRUE)
  expect_equal(f$n_presence, 100L)
  expect_equal(f$n_background, 1000L)
})

test_that("test_frac = 0 skips evaluation", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  f <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, test_frac = 0, quiet = TRUE)
  expect_true(all(is.na(f$metrics$auc)))
  expect_error(fit_sdm(two_species_db(), "Genus_low", clim, test_frac = 1),
               "test_frac")
})

test_that("a model is fitted once and projected onto many slices", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 950))
  f <- fit_sdm(two_species_db(), "Genus_high", clim, land = fake_land(),
               num_trees = 50, quiet = TRUE)
  p1 <- project_sdm(f, clim, 850, quiet = TRUE)
  p2 <- project_sdm(f, clim, 950, quiet = TRUE)
  expect_s3_class(p1, "richcast_projection")
  expect_equal(p1$threshold, f$threshold)
  expect_false(isTRUE(all.equal(terra::values(suitability(p1)),
                                terra::values(suitability(p2)))))
  expect_error(project_sdm(list(), clim, 850), "fit_sdm")
})

test_that("a fitted model survives saveRDS", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  f <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, quiet = TRUE)
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(f, path)
  g <- readRDS(path)
  expect_equal(terra::values(suitability(g)), terra::values(suitability(f)))
  expect_equal(project_sdm(g, clim, 850, quiet = TRUE)$cells,
               project_sdm(f, clim, 850, quiet = TRUE)$cells)
})

test_that("models use the eight canonical bioclim variables and nothing else", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  expect_equal(bioclim_vars, c("bio01", "bio04", "bio05", "bio06",
                               "bio12", "bio15", "bio16", "bio17"))
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  # A ninth variable on disk is ignored by default...
  for (lbl in c("present_1985", "time_0850")) {
    b <- terra::rast(file.path(dir, lbl, "bio01.tif"))
    names(b) <- "bio02"
    terra::writeRaster(b, file.path(dir, lbl, "bio02.tif"))
  }
  f <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, quiet = TRUE)
  expect_equal(f$predictors, bioclim_vars)
  # ...and refused if asked for.
  expect_error(fit_sdm(two_species_db(), "Genus_low", clim,
                       predictors = c("bio01", "bio02")), "canonical")
  expect_error(prepare_climate(withr::local_tempdir(), vars = "bio19",
                               times = 850, extent = c(0, 1, 0, 1)),
               "canonical")
})

test_that("MaxEnt can take background from the whole extent; the forest keeps its pseudo-absences", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()

  fd <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50,
                quiet = TRUE)
  f0 <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50,
                maxnet_background = "outside", quiet = TRUE)
  fe <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50,
                maxnet_background = "extent", quiet = TRUE)

  # The default is the extent design, recorded on the model, and naming it
  # changes nothing; "outside" is recorded too.
  expect_equal(fd$maxnet_background, "extent")
  expect_equal(terra::values(suitability(fd)), terra::values(suitability(fe)))
  expect_equal(fd$metrics, fe$metrics)
  expect_equal(fd$cutoffs, fe$cutoffs)
  expect_equal(f0$maxnet_background, "outside")
  expect_equal(f0$n_maxnet_background, f0$n_background)

  # The extent holds fewer land cells than either sample asks for, so the
  # pseudo-absences are every cell outside the range and MaxEnt's background
  # is every cell, the range's included.
  expect_equal(fe$maxnet_background, "extent")
  expect_equal(fe$n_background, f0$n_background)
  expect_equal(fe$n_maxnet_background, fe$n_background + fe$range_cells)

  # Only MaxEnt changes: the forest is fitted on the same points as before.
  d <- as.data.frame(terra::rast(file.path(dir, "present_1985",
                                            paste0(bioclim_vars, ".tif"))))
  p0 <- richcast:::predict_members(f0$members, d)
  pe <- richcast:::predict_members(fe$members, d)
  expect_equal(pe[, "rf"], p0[, "rf"])
  # maxnet adds the presences to its background, and here they are every
  # range cell, so on this coarse grid the two designs give MaxEnt the same
  # points; the fine-grid test below shows MaxEnt changing.
  expect_equal(pe[, "maxent"], p0[, "maxent"])
  expect_equal(fe$metrics$auc[fe$metrics$model == "rf"],
               f0$metrics$auc[f0$metrics$model == "rf"])

  expect_error(fit_sdm(db, "Genus_low", clim, maxnet_background = "inside"),
               "should be one of")
})

test_that("on a fine grid, extent background changes MaxEnt only, up to its cap", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  # The range holds far more cells than the 100 pseudo-presences, so the
  # extent background reaches cells inside it that the default never sees.
  r <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10, ymin = 0,
                   ymax = 10, crs = "EPSG:4326")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  for (lbl in c("present_1985", "time_0850")) {
    d <- file.path(dir, lbl)
    dir.create(d)
    b1 <- r; terra::values(b1) <- xy[, 1] * 3 + xy[, 2]; names(b1) <- "bio01"
    b2 <- r; terra::values(b2) <- xy[, 2] * 2 - xy[, 1]; names(b2) <- "bio12"
    terra::writeRaster(b1, file.path(d, "bio01.tif"))
    terra::writeRaster(b2, file.path(d, "bio12.tif"))
  }
  clim <- climate_dir(dir)
  fit <- function(...) {
    fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
            predictors = c("bio01", "bio12"), num_trees = 50, quiet = TRUE, ...)
  }
  f0 <- fit(maxnet_background = "outside")
  fe <- fit(maxnet_background = "extent", n_maxnet_background = 2000)
  expect_equal(fe$n_maxnet_background, 2000L)
  d <- as.data.frame(terra::rast(file.path(dir, "present_1985",
                                            c("bio01.tif", "bio12.tif"))))
  p0 <- richcast:::predict_members(f0$members, d)
  pe <- richcast:::predict_members(fe$members, d)
  expect_equal(pe[, "rf"], p0[, "rf"])
  expect_false(isTRUE(all.equal(pe[, "maxent"], p0[, "maxent"])))
})

test_that("regmult sets MaxEnt's regularization and leaves the forest alone", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = 850)
  db <- two_species_db()
  f1 <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50, quiet = TRUE)
  f2 <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50, regmult = 1,
                quiet = TRUE)
  f5 <- fit_sdm(db, "Genus_low", clim, land = fake_land(), num_trees = 50, regmult = 5,
                quiet = TRUE)
  expect_equal(f1$regmult, 1)
  expect_equal(f5$regmult, 5)
  expect_equal(terra::values(suitability(f2)), terra::values(suitability(f1)))

  d <- as.data.frame(terra::rast(file.path(dir, "present_1985",
                                            paste0(bioclim_vars, ".tif"))))
  p1 <- richcast:::predict_members(f1$members, d)
  p5 <- richcast:::predict_members(f5$members, d)
  expect_equal(p5[, "rf"], p1[, "rf"])
  expect_false(isTRUE(all.equal(p5[, "maxent"], p1[, "maxent"])))
  # Stronger regularization keeps fewer features.
  expect_lte(length(f5$members$maxent$betas), length(f1$members$maxent$betas))

  expect_error(fit_sdm(db, "Genus_low", clim, regmult = 0), "regmult")
  expect_error(fit_sdm(db, "Genus_low", clim, regmult = c(1, 2)), "regmult")
})
