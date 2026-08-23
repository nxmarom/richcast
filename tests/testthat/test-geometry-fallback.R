# Published range maps routinely carry rings that are valid in the plane but
# self-intersecting on the sphere. sf's s2 backend rejects those outright, so
# every geometry predicate richcast runs has to survive them.

test_that("with_planar_fallback returns the spherical result when it works", {
  expect_equal(richcast:::with_planar_fallback(function() 42), 42)
  # s2 setting is left as it was found.
  before <- sf::sf_use_s2()
  richcast:::with_planar_fallback(function() 1)
  expect_equal(sf::sf_use_s2(), before)
})

test_that("with_planar_fallback retries in planar mode and restores s2", {
  before <- sf::sf_use_s2()
  attempts <- 0

  out <- NULL
  expect_message(
    out <- richcast:::with_planar_fallback(function() {
      attempts <<- attempts + 1
      if (sf::sf_use_s2()) stop("s2 refused this geometry")
      "planar result"
    }),
    "planar mode"
  )

  expect_equal(out, "planar result")
  expect_equal(attempts, 2)
  expect_equal(sf::sf_use_s2(), before)
})

test_that("s2 setting is restored even when the planar retry also fails", {
  before <- sf::sf_use_s2()
  expect_error(
    suppressMessages(
      richcast:::with_planar_fallback(function() stop("broken either way"))
    ),
    "broken either way"
  )
  expect_equal(sf::sf_use_s2(), before)
})

test_that("quiet suppresses the fallback warning", {
  expect_silent(
    richcast:::with_planar_fallback(
      function() if (sf::sf_use_s2()) stop("nope") else TRUE,
      quiet = TRUE
    )
  )
})
