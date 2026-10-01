test_that("fossil_presences reads wide and long tables alike", {
  wide <- data.frame(species = "Genus low", unit_id = c("A", "B"),
                     lon = c(8.25, 8.75), lat = 1.25,
                     slices = c("1.1:3;1.2:1", "1.2:1"))
  f <- fossil_presences(wide)
  expect_s3_class(f, "richcast_fossils")
  expect_equal(f$species, rep("Genus_low", 3))
  expect_equal(f$ka_bp, c(1.1, 1.2, 1.2))
  expect_equal(f$prob, c(0.75, 0.25, 1))

  long <- data.frame(species = "Genus_low", unit_id = c("A", "A", "B"),
                     lon = c(8.25, 8.25, 8.75), lat = 1.25,
                     ka_bp = c(1.1, 1.2, 1.2), weight = c(3, 1, 5))
  expect_equal(fossil_presences(long)$prob, c(0.75, 0.25, 1))
  expect_identical(fossil_presences(f), f)
  expect_error(fossil_presences(data.frame(species = "x", lon = 1, lat = 1)),
               "slices")
  expect_equal(richcast:::ka_to_time(c(0, 1.1, 120)),
               c("present", "850", "-118050"))
})

test_that("fit_sdm without fossils is unchanged", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 750))
  a <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, quiet = TRUE)
  b <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, fossils = NULL, quiet = TRUE)
  expect_identical(a$threshold, b$threshold)
  expect_identical(a$metrics, b$metrics)
  expect_identical(a$study_extent, b$study_extent)
  expect_null(a$n_fossil)
  # Fossils of another species leave the fit alone too.
  other <- data.frame(species = "Genus_high", lon = 1, lat = 1, slices = "1.1:1")
  c <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, fossils = other, quiet = TRUE)
  expect_identical(a$metrics, c$metrics)
})

test_that("fossil presences widen the extent and the niche, not the threshold rule", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 750))
  # Two units far outside the range, in the slices 850 and 750 CE.
  fos <- data.frame(species = "Genus_low", unit_id = c("A", "B"),
                    lon = c(8.25, 8.75), lat = c(1.25, 1.75),
                    slices = c("1.1:0.5;1.2:0.5", "1.2:1"))
  base <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
                  num_trees = 50, quiet = TRUE)
  f <- fit_sdm(two_species_db(), "Genus_low", clim, land = fake_land(),
               num_trees = 50, fossils = fos, quiet = TRUE)

  expect_equal(f$n_fossil, 100L)
  # 1000 asked for; the 20 x 20 test grid has fewer distinct cells per slice.
  expect_gt(f$n_fossil_background, 100L)
  expect_lte(f$n_fossil_background, 1000L)
  expect_equal(f$fossil_units, 2L)
  # Draws are shared equally between units, and stay within each unit's slices.
  expect_equal(as.vector(table(f$fossil_draws$unit_id)), c(50L, 50L))
  expect_true(all(f$fossil_draws$time[f$fossil_draws$unit_id == "B"] == "750"))
  expect_setequal(unique(f$fossil_draws$time), c("850", "750"))
  # The study extent now reaches the sites.
  expect_gt(f$study_extent[2], 8.75)
  expect_lt(base$study_extent[2], 8.25)
  # p10 still comes from the range's pseudo-presences only.
  pres <- terra::extract(suitability(f), cbind(8.25, 1.25))[1, 1]
  expect_false(is.na(pres))
  expect_equal(f$threshold_rule, "p10")

  pt <- richness_at(f, 8.25, 1.25, time = 850, climate = clim)
  expect_gt(pt$species$suitability, f$threshold)

  s <- cli::cli_fmt(print(f))
  expect_true(any(grepl("fossil presences", s)))
})

test_that("run_hindcast_series passes fossils to the species they name", {
  skip_if_not_installed("maxnet")
  skip_if_not_installed("ranger")
  dir <- withr::local_tempdir()
  clim <- structured_climate(dir, times = c(850, 750))
  fos <- data.frame(species = c("Genus_low", "Genus_nope"), lon = 8.25, lat = 1.25,
                    slices = "1.1:1")
  expect_warning(
    res <- run_hindcast_series(two_species_db(), clim, times = 850,
                               region = region(c(0, 10, 0, 10)),
                               species = c("Genus_low", "Genus_high"),
                               merge = list(), fossils = fos,
                               land = fake_land(), num_trees = 50, quiet = TRUE),
    "Genus_nope"
  )
  expect_equal(res$models$n_fossil[res$models$species == "Genus_low"], 100L)
  expect_true(is.na(res$models$n_fossil[res$models$species == "Genus_high"]))
})
