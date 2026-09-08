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

test_that("focus_cells_each counts each range separately, clipped to the focus", {
  focus <- focus_box(c(0, 4, 0, 4), label = "test")
  ranges <- list(
    inside  = sf::st_sfc(square(0, 0, 2), crs = 4326),   # wholly within
    partial = sf::st_sfc(square(3, 3, 2), crs = 4326),   # half hangs out
    outside = sf::st_sfc(square(50, 50, 2), crs = 4326), # nowhere near
    gone    = NULL
  )
  n <- focus_cells_each(ranges, focus, resolution = 0.5)

  expect_named(n, c("inside", "partial", "outside", "gone"))
  expect_type(n, "integer")
  # A 2-degree square on a 0.5-degree grid is 4x4 cells; touches = TRUE adds
  # the row and column the boundary grazes, so assert the ordering rather than
  # exact counts, which is what the column is actually read for.
  expect_gt(n[["inside"]], 0L)
  expect_lt(n[["partial"]], n[["inside"]])
  # The species that never enters the focus contributes nothing -- the case
  # that made a trajectory panel appear for a species absent from the region.
  expect_equal(n[["outside"]], 0L)
  expect_equal(n[["gone"]], 0L)
})

test_that("focus_cells_each and richness_stack agree on the same ranges", {
  focus <- focus_box(c(0, 4, 0, 4))
  ranges <- list(
    a = sf::st_sfc(square(0, 0, 3), crs = 4326),
    b = sf::st_sfc(square(2, 0, 2), crs = 4326),
    c = NULL
  )
  r <- richness_stack(ranges, focus, resolution = 0.5, quiet = TRUE)
  n <- focus_cells_each(ranges, focus, resolution = 0.5)

  # Richness is the sum over species of exactly these counts, so the totals
  # must match. If they drift, the two columns are being measured differently
  # and reading focus_cells against $richness stops being valid.
  expect_equal(sum(n), sum(terra::values(r, mat = FALSE, na.rm = TRUE)))
})

test_that("focus_cells_each survives an empty range list", {
  expect_length(focus_cells_each(list(), focus_box(c(0, 4, 0, 4))), 0)
})

test_that("richness_stack rejects a bare bbox in place of a focus", {
  expect_error(
    richness_stack(list(), c(0, 4, 0, 4)),
    "focus_box"
  )
})

test_that("a subregion outside the focus yields NA rather than an error", {
  # A site drifting outside the study region is a legitimate finding, not a
  # crash -- but it must be visible as zero cells, not a plausible number.
  focus <- focus_box(c(0, 4, 0, 4), label = "big")
  ranges <- list(a = sf::st_sfc(square(0, 0, 2), crs = 4326))
  r <- richness_stack(ranges, focus, resolution = 0.5, quiet = TRUE)

  expect_warning(
    s <- richness_stats(r, time = 1350, focus = focus,
                        subregions = list(elsewhere = focus_box(c(50, 51, 50, 51)))),
    "does not overlap"
  )
  expect_equal(s$elsewhere_cells, 0L)
  expect_true(is.na(s$elsewhere_max))
})

test_that("several subregions each get their own columns", {
  focus <- focus_box(c(0, 4, 0, 4))
  ranges <- list(a = sf::st_sfc(square(0, 0, 2), crs = 4326))
  r <- richness_stack(ranges, focus, resolution = 0.5, quiet = TRUE)

  s <- richness_stats(r, time = 1350, focus = focus, subregions = list(
    west = focus_box(c(0.1, 1.9, 0.1, 1.9)),
    east = focus_box(c(2.1, 3.9, 0.1, 1.9))
  ))
  expect_true(all(c("west_max", "west_mean", "west_cells",
                    "east_max", "east_mean", "east_cells") %in% names(s)))
  # The range covers the west box only.
  expect_equal(s$west_max, 1)
  expect_equal(s$east_max, 0)
})
