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
  expect_error(richcast:::resolve_threshold("tss", obs, pred), "p10")
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
  d <- data.frame(bio01 = c(5, 20), bio12 = c(10, 30))
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
               num_trees = 50, quiet = TRUE)
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
