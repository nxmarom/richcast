# ==============================================================================
# Assemblage richness
#
# Richness is built from the models' own rasters rather than from polygons:
# each species' ensemble suitability is laid onto one grid covering the
# region, thresholded, and summed. The grid is the climate grid itself, cropped
# to the region, so nothing is resampled across resolutions. The same pass
# sums the raw suitabilities into the deprecated `expected` layer, which is
# not a richness estimate because suitabilities are not calibrated
# probabilities.
# ==============================================================================

#' Maximum distance, in degrees, between a species' range and the region
#' @noRd
near_distance <- 10

#' Species whose ranges lie near a region
#'
#' The species richcast will train for a region: those whose present-day range
#' polygon comes within 10 degrees of the region's box. The distance is fixed,
#' and [run_hindcast_series()] always applies it, so a species far from the
#' region can never enter its richness -- a model extrapolated to a continent
#' the species has never occupied will happily place a llama in Britain.
#'
#' Distance is measured on the box expanded by 10 degrees in longitude and
#' latitude, which is quick and conservative near the poles.
#'
#' @param db A `richcast_db`.
#' @param region A [region()].
#' @param quiet Suppress the report.
#' @return A character vector of species names.
#' @examples
#' db <- build_taxon_db(
#'   sf::st_sf(
#'     species = c("Genus_near", "Genus_far"),
#'     geometry = sf::st_sfc(
#'       sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
#'       sf::st_polygon(list(cbind(c(60, 61, 61, 60, 60), c(0, 0, 1, 1, 0)))),
#'       crs = 4326
#'     )
#'   ),
#'   quiet = TRUE
#' )
#' species_near(db, region(c(5, 10, 0, 5)))
#' @export
species_near <- function(db, region, quiet = FALSE) {
  check_region(region)
  d <- near_distance
  b <- region$box
  wide <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = max(b[1] - d, -180), xmax = min(b[2] + d, 180),
      ymin = max(b[3] - d, -90),  ymax = min(b[4] + d, 90)),
    crs = sf::st_crs(4326)
  ))
  hits <- with_planar_fallback(
    function() {
      suppressMessages(lengths(sf::st_intersects(sf::st_geometry(db), wide)) > 0)
    },
    what = "overlap test", quiet = quiet
  )
  if (!quiet) {
    cli::cli_alert_info(
      "{sum(hits)}/{nrow(db)} species lie within {d} deg of {.emph {region$label}}."
    )
  }
  db$species[hits]
}

#' Empty accumulator for one slice's richness
#' @noRd
new_stack <- function(template) {
  zero <- terra::rast(template)
  terra::values(zero) <- 0
  list(richness = zero, expected = zero, covered = zero, n = 0L)
}

#' Add one taxon's suitability to a slice's richness
#'
#' A taxon is one species or several merged ones. For a merged taxon each
#' member's surface is laid on the grid, and in every cell the member furthest
#' above its own threshold stands for the taxon: the taxon is present where
#' any member is, and contributes that member's suitability to `expected`.
#'
#' @param suits List of member suitability `SpatRaster`s.
#' @param thresholds Numeric vector of the members' thresholds.
#' @return The updated accumulator, with `region_cells` -- this taxon's
#'   above-threshold cells inside the region -- as an attribute.
#' @noRd
stack_add <- function(acc, suits, thresholds) {
  if (inherits(suits, "SpatRaster")) suits <- list(suits)
  on_grid <- lapply(suits, function(s) tryCatch(
    terra::resample(s, acc$richness, method = "near"),
    error = function(e) NULL
  ))
  keep <- !vapply(on_grid, is.null, logical(1))
  if (!any(keep)) {
    # No member's study extent reaches the region at all.
    attr(acc, "region_cells") <- 0L
    return(acc)
  }
  s <- terra::rast(on_grid[keep])
  thr <- unname(thresholds[keep])
  if (terra::nlyr(s) == 1) {
    suit <- s
    present <- s > thr
  } else {
    margin <- s - thr
    best <- terra::which.max(margin)
    suit <- terra::selectRange(s, best)
    present <- terra::app(margin, max, na.rm = TRUE) > 0
  }
  acc$richness <- acc$richness + terra::ifel(is.na(present), 0, present)
  acc$expected <- acc$expected + terra::ifel(is.na(suit), 0, suit)
  acc$covered  <- acc$covered + !is.na(suit)
  acc$n <- acc$n + 1L
  n <- terra::global(present, "sum", na.rm = TRUE)[1, 1]
  attr(acc, "region_cells") <- if (is.na(n)) 0L else as.integer(n)
  acc
}

#' Finish an accumulator into a two-layer surface
#'
#' Cells no species' model covered -- sea, or land outside every study
#' extent -- are NA rather than zero: nothing was predicted there, which is not
#' the same as predicting no species.
#' @noRd
stack_finish <- function(acc) {
  r <- c(acc$richness, acc$expected)
  names(r) <- c("richness", "expected")
  terra::mask(r, acc$covered, maskvalues = 0)
}

#' Summarise one richness surface as a table row
#' @noRd
richness_row <- function(surface, period, time, n_species) {
  rv <- terra::values(surface[["richness"]], mat = FALSE, na.rm = TRUE)
  ev <- terra::values(surface[["expected"]], mat = FALSE, na.rm = TRUE)
  have <- length(rv) > 0
  tibble::tibble(
    period          = period,
    time            = time,
    cells           = length(rv),
    mean_richness   = if (have) mean(rv) else NA_real_,
    median_richness = if (have) stats::median(rv) else NA_real_,
    max_richness    = if (have) max(rv) else NA_real_,
    mean_expected   = if (have) mean(ev) else NA_real_,
    species_modelled = n_species
  )
}
