# ==============================================================================
# Modelling regions
#
# A region is the area richness is mapped and summarised over, and the area
# that decides which species are worth training at all. Kept as sf throughout:
# terra SpatVectors hold external pointers and do not survive saveRDS(), so
# any object a user might reasonably save stores sf and converts on demand.
# ==============================================================================

#' Preset modelling regions
#'
#' Bounding boxes, as `c(xmin, xmax, ymin, ymax)` in decimal degrees, for the
#' regions [region()] accepts by name. They are deliberately generous and
#' overlap where continents meet (the Middle East sits inside both Asia's and
#' Africa's neighbourhoods); pass your own box to [region()] for anything
#' tighter.
#'
#' @format A named list of numeric vectors of length 4.
#' @examples
#' names(region_presets)
#' region_presets$europe
#' @export
region_presets <- list(
  europe        = c(-25,   45,  34,  72),
  asia          = c( 25,  180, -11,  82),
  middle_east   = c( 25,   63,  12,  42),
  africa        = c(-20,   55, -36,  38),
  north_america = c(-170, -50,   7,  84),
  south_america = c(-92,  -30, -57,  13)
)

#' Define a modelling region
#'
#' Either one of the preset regions, by name, or a user-defined box.
#'
#' @param x A preset name -- one of `"europe"`, `"asia"`, `"middle_east"`,
#'   `"africa"`, `"north_america"`, `"south_america"` -- or a numeric box
#'   `c(xmin, xmax, ymin, ymax)` in decimal degrees. Case, spaces and hyphens
#'   in names are ignored, so `"North America"` works.
#' @param label Optional name, used in output tables. Defaults to the preset
#'   name, or `"custom"` for a box.
#' @return A `richcast_region`.
#' @seealso [region_presets]
#' @examples
#' region("europe")
#' region("Middle East")
#' region(c(68, 87, 39, 46), label = "Tian Shan")
#' @export
region <- function(x, label = NULL) {
  if (is.character(x)) {
    if (length(x) != 1) rc_abort("{.arg x} must be a single region name.")
    key <- gsub("[ -]+", "_", tolower(trimws(x)))
    if (!key %in% names(region_presets)) {
      rc_abort(c(
        "Unknown region {.val {x}}.",
        "i" = "Presets: {.val {names(region_presets)}}.",
        "i" = "Or pass a box as {.code c(xmin, xmax, ymin, ymax)}."
      ))
    }
    coords <- region_presets[[key]]
    label <- label %||% key
  } else {
    check_extent(x, arg = "x")
    coords <- x
    label <- label %||% "custom"
  }
  bb <- sf::st_bbox(
    c(xmin = coords[1], xmax = coords[2], ymin = coords[3], ymax = coords[4]),
    crs = sf::st_crs(4326)
  )
  structure(
    list(geometry = sf::st_as_sfc(bb), box = unname(coords), label = label),
    class = "richcast_region"
  )
}

#' @export
print.richcast_region <- function(x, ...) {
  b <- x$box
  cli::cli_text("{.cls richcast_region} {.emph {x$label}}")
  cli::cli_text("  lon [{b[1]}, {b[2]}] x lat [{b[3]}, {b[4]}]")
  invisible(x)
}

#' Check an argument is a region
#' @noRd
check_region <- function(x, arg = "region") {
  if (!inherits(x, "richcast_region")) {
    rc_abort("{.arg {arg}} must come from {.fn region}.")
  }
  invisible(x)
}

#' Coerce a region to a terra SpatVector
#' @noRd
region_vect <- function(region) {
  terra::vect(sf::st_sf(geometry = region$geometry))
}

#' Coerce a region to a terra SpatExtent
#' @noRd
region_ext <- function(region) {
  b <- region$box
  terra::ext(b[1], b[2], b[3], b[4])
}

#' Natural Earth land outline, degrading gracefully
#'
#' `scale = "large"` needs rnaturalearthhires, which is not on CRAN; an absent
#' package costs resolution, not a failed run.
#' @noRd
land_outline <- function(scale = c("medium", "large")) {
  scale <- match.arg(scale)
  rlang::check_installed("rnaturalearth", "to mask predictions to land.")

  if (scale == "large" && !requireNamespace("rnaturalearthhires", quietly = TRUE)) {
    cli::cli_warn(c(
      "{.pkg rnaturalearthhires} is not installed; falling back to {.val medium} resolution.",
      "i" = 'Install it with {.code install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")}.'
    ))
    scale <- "medium"
  }

  world <- rnaturalearth::ne_countries(scale = scale, returnclass = "sf")
  world <- world[!world$name_long %in% "Antarctica", ]
  sf::st_geometry(world)
}
