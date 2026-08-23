test_that("focus_box validates and stores coordinates", {
  f <- focus_box(c(68, 87, 39, 46), label = "Tian Shan")
  expect_s3_class(f, "richcast_focus")
  expect_equal(f$label, "Tian Shan")
  bb <- sf::st_bbox(f$geometry)
  expect_equal(as.numeric(bb[c("xmin", "ymin", "xmax", "ymax")]), c(68, 39, 87, 46))
})

test_that("focus_box rejects reversed or malformed coords", {
  # A c(xmin, ymin, xmax, ymax) ordering mistake is easy to make and produces
  # a silently empty focus, so it must be an error rather than a warning.
  expect_error(focus_box(c(68, 39, 87, 46)), "xmin < xmax")
  expect_error(focus_box(c(1, 2, 3)), "length 4")
  expect_error(focus_box("a"), "numeric")
})

test_that("focus_polygon accepts sf, sfc and SpatVector", {
  poly <- sf::st_sfc(square(0, 0, 4), crs = 4326)
  expect_s3_class(focus_polygon(poly), "richcast_focus")
  expect_s3_class(focus_polygon(sf::st_sf(geometry = poly)), "richcast_focus")
  expect_s3_class(
    focus_polygon(terra::vect(sf::st_sf(geometry = poly))),
    "richcast_focus"
  )
})

test_that("focus helpers convert to terra without losing the extent", {
  f <- focus_box(c(-10, 10, -5, 5))
  e <- richcast:::focus_ext(f)
  expect_equal(unname(as.numeric(as.vector(e))), c(-10, 10, -5, 5))
  expect_s4_class(richcast:::focus_vect(f), "SpatVector")
})
