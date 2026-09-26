test_that("every preset region builds and is a valid box", {
  for (nm in names(region_presets)) {
    r <- region(nm)
    expect_s3_class(r, "richcast_region")
    expect_equal(r$label, nm)
    b <- r$box
    expect_true(b[1] < b[2] && b[3] < b[4])
  }
  expect_setequal(
    names(region_presets),
    c("europe", "asia", "middle_east", "africa", "north_america", "south_america")
  )
})

test_that("preset names are matched loosely", {
  expect_equal(region("Middle East")$box, region_presets$middle_east)
  expect_equal(region("north-america")$box, region_presets$north_america)
  expect_error(region("atlantis"), "Unknown region")
})

test_that("a user-defined box is validated", {
  r <- region(c(68, 87, 39, 46), label = "Tian Shan")
  expect_equal(r$label, "Tian Shan")
  expect_equal(as.vector(richcast:::region_ext(r)), c(xmin = 68, xmax = 87, ymin = 39, ymax = 46))
  expect_equal(region(c(0, 1, 0, 1))$label, "custom")
  expect_error(region(c(87, 68, 39, 46)), "xmin < xmax")
  expect_error(region(c(68, 39, 87, 46)), "bbox order")
})

test_that("check_extent diagnoses bbox order and bad latitudes", {
  expect_error(richcast:::check_extent(c(68, 39, 87, 46)), "Did you mean")
  expect_error(richcast:::check_extent(c(0, 10, -100, 50)), "outside")
  expect_invisible(richcast:::check_extent(c(0, 10, 0, 10)))
})

test_that("species_near keeps ranges within the distance and drops the rest", {
  db <- build_taxon_db(
    sf::st_sf(
      species = c("Genus_inside", "Genus_near", "Genus_far"),
      geometry = sf::st_sfc(square(6, 1), square(18, 1), square(40, 1),
                            crs = 4326)
    ),
    quiet = TRUE
  )
  reg <- region(c(5, 10, 0, 5))
  expect_setequal(species_near(db, reg, quiet = TRUE),
                  c("Genus_inside", "Genus_near"))
  expect_equal(species_near(db, reg, distance = 0, quiet = TRUE), "Genus_inside")
  expect_error(species_near(db, c(5, 10, 0, 5)), "region")
})
