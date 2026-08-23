test_that("overlapping ranges sum to the right richness", {
  focus <- focus_box(c(0, 4, 0, 4), label = "test")
  ranges <- list(
    a = sf::st_sfc(square(0, 0, 3), crs = 4326),   # covers 0-3
    b = sf::st_sfc(square(2, 0, 2), crs = 4326),   # covers 2-4; overlaps a on 2-3
    c = NULL                                        # modelled to extinction
  )
  r <- richness_stack(ranges, focus, resolution = 0.5, quiet = TRUE)

  # Peak richness is 2 where a and b overlap, never 3.
  expect_equal(max(terra::values(r, mat = FALSE, na.rm = TRUE)), 2)
  expect_true(min(terra::values(r, mat = FALSE, na.rm = TRUE)) >= 0)
})

test_that("richness_stack reports species that contributed nothing", {
  focus <- focus_box(c(0, 4, 0, 4))
  ranges <- list(a = sf::st_sfc(square(0, 0, 2), crs = 4326), gone = NULL)
  expect_message(richness_stack(ranges, focus, resolution = 1), "contributed no range")
})

test_that("richness_stats returns one tidy row and keeps time numeric", {
  focus <- focus_box(c(0, 4, 0, 4), label = "test")
  ranges <- list(a = sf::st_sfc(square(0, 0, 2), crs = 4326))
  r <- richness_stack(ranges, focus, resolution = 0.5, quiet = TRUE)

  s <- richness_stats(r, time = 1350, focus = focus)
  expect_s3_class(s, "tbl_df")
  expect_equal(nrow(s), 1)
  expect_type(s$time, "double")
  expect_equal(s$focus, "test")
  expect_true(s$max_richness == 1)
})

test_that("subregions add their own columns", {
  focus <- focus_box(c(0, 4, 0, 4), label = "big")
  ranges <- list(a = sf::st_sfc(square(0, 0, 2), crs = 4326))
  r <- richness_stack(ranges, focus, resolution = 0.5, quiet = TRUE)

  s <- richness_stats(
    r, time = 1350, focus = focus,
    subregions = list(site = focus_box(c(0.1, 0.9, 0.1, 0.9)))
  )
  expect_true(all(c("site_max", "site_mean", "site_cells") %in% names(s)))
  expect_equal(s$site_max, 1)
})

test_that("richness_stack rejects a bare bbox in place of a focus", {
  expect_error(
    richness_stack(list(), c(0, 4, 0, 4)),
    "focus_box"
  )
})
