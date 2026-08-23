# ==============================================================================
# Geographic foci
#
# A focus is the region richness is summarised over. Kept as sf throughout:
# terra SpatVectors hold external pointers and do not survive saveRDS(), so
# any object a user might reasonably save stores sf and converts on demand.
# ==============================================================================

new_focus <- function(geom, type, label = NULL) {
  structure(
    list(geometry = geom, type = type, label = label %||% type),
    class = "richcast_focus"
  )
}

#' Define a geographic focus
#'
#' The region that richness is cropped, masked and summarised to. Three
#' constructors cover the usual cases.
#'
#' @param coords Numeric vector `c(xmin, xmax, ymin, ymax)` in decimal degrees.
#' @param label Optional name, used in output tables and plot titles.
#' @return A `richcast_focus`.
#' @examples
#' tianshan <- focus_box(c(68, 87, 39, 46), label = "Tian Shan")
#' tianshan
#' @export
focus_box <- function(coords, label = NULL) {
  check_extent(coords, arg = "coords")
  bb <- sf::st_bbox(
    c(xmin = coords[1], xmax = coords[2], ymin = coords[3], ymax = coords[4]),
    crs = sf::st_crs(4326)
  )
  new_focus(sf::st_as_sfc(bb), "box", label)
}

#' @rdname focus_box
#' @param scale Natural Earth resolution: `"medium"` (no extra dependency) or
#'   `"large"` (needs \pkg{rnaturalearthhires}, which lives on r-universe
#'   rather than CRAN). Falls back to `"medium"` with a warning if the hires
#'   data is unavailable.
#' @export
focus_global <- function(scale = c("medium", "large"), label = NULL) {
  scale <- match.arg(scale)
  new_focus(land_outline(scale), "global", label)
}

#' @rdname focus_box
#' @param x An `sf`, `sfc` or `SpatVector` object.
#' @export
focus_polygon <- function(x, label = NULL) {
  geom <- if (inherits(x, "SpatVector")) {
    sf::st_as_sfc(sf::st_as_sf(x))
  } else if (inherits(x, "sf")) {
    sf::st_geometry(x)
  } else if (inherits(x, "sfc")) {
    x
  } else {
    rc_abort("{.arg x} must be {.cls sf}, {.cls sfc} or {.cls SpatVector}.")
  }
  if (is.na(sf::st_crs(geom))) sf::st_crs(geom) <- 4326
  geom <- sf::st_transform(geom, 4326)
  new_focus(geom, "polygon", label)
}

#' @export
print.richcast_focus <- function(x, ...) {
  bb <- sf::st_bbox(x$geometry)
  cli::cli_text("{.cls richcast_focus} <{x$type}> {.emph {x$label}}")
  cli::cli_text(
    "  bbox: [{round(bb[1], 2)}, {round(bb[3], 2)}] x [{round(bb[2], 2)}, {round(bb[4], 2)}]"
  )
  invisible(x)
}

#' Coerce a focus to a terra SpatVector
#' @noRd
focus_vect <- function(focus) {
  terra::vect(sf::st_sf(geometry = focus$geometry))
}

#' Coerce a focus to a terra SpatExtent
#' @noRd
focus_ext <- function(focus) {
  bb <- sf::st_bbox(focus$geometry)
  terra::ext(bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]])
}

#' Natural Earth land outline, degrading gracefully
#'
#' `scale = "large"` needs rnaturalearthhires, which is not on CRAN. The source
#' pipeline hard-required it; here an absent package costs resolution, not a
#' failed run.
#' @noRd
land_outline <- function(scale = c("medium", "large")) {
  scale <- match.arg(scale)
  rlang::check_installed("rnaturalearth", "to build a global focus.")

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
