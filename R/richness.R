# ==============================================================================
# Assemblage richness
# ==============================================================================

#' Stack species ranges into a richness surface
#'
#' Rasterises each species range onto a common grid and sums them, giving the
#' number of species whose modelled range covers each cell.
#'
#' @param ranges A named list of `sf`/`sfc` geometries, one per species.
#'   `NULL` entries (species with no suitable area) are skipped and counted.
#' @param focus A [focus_box()] / [focus_global()] / [focus_polygon()] object.
#'   The grid covers this region only; there is no reason to rasterise a
#'   global grid to summarise one valley.
#' @param resolution Cell size in degrees.
#' @param touches Count a cell as occupied if the range touches it at all.
#'   `TRUE` matches the source pipeline and slightly inflates small ranges.
#' @param quiet Suppress progress messages.
#' @return A `SpatRaster` of species counts, cropped and masked to `focus`.
#' @seealso [richness_stats()]
#' @export
richness_stack <- function(ranges,
                           focus,
                           resolution = 0.1,
                           touches = TRUE,
                           quiet = FALSE) {

  if (!inherits(focus, "richcast_focus")) {
    rc_abort("{.arg focus} must come from {.fn focus_box}, {.fn focus_global} or {.fn focus_polygon}.")
  }

  fe <- focus_ext(focus)
  template <- terra::rast(fe, resolution = resolution, crs = "EPSG:4326")
  terra::values(template) <- 0

  kept <- 0L
  empty <- character(0)

  for (nm in names(ranges)) {
    g <- ranges[[nm]]
    if (is.null(g) || length(g) == 0) {
      empty <- c(empty, nm)
      next
    }
    v <- terra::vect(sf::st_sf(geometry = sf::st_geometry(g)))
    layer <- terra::rasterize(v, template, field = 1, background = 0,
                              touches = touches)
    template <- template + layer
    kept <- kept + 1L
  }

  if (!quiet) {
    cli::cli_alert_info("Stacked {kept} range{?s} at {resolution} deg resolution.")
    if (length(empty) > 0) {
      cli::cli_alert_warning(
        "{length(empty)} species contributed no range: {.val {empty}}."
      )
    }
  }

  out <- terra::mask(template, focus_vect(focus))
  names(out) <- "richness"
  out
}

#' Summarise a richness surface
#'
#' @details
#' # Subregions
#'
#' `subregions` reports a small area without shrinking the analysis to it.
#' Richness is computed once, over `focus`; each subregion is then a second
#' readout of that same surface -- cropped, masked, summarised -- so it costs
#' a crop rather than another model run.
#'
#' This matters when the area of interest is small relative to the climate
#' grid. A 0.7-degree site box against a 0.5-degree reconstruction contains
#' about four cells, and a mean over four cells is not a regional signal. Make
#' the focus regional and demote the site to a subregion, and one run gives
#' both:
#'
#' ```r
#' run_hindcast_series(
#'   db, clim, times = times,
#'   focus      = focus_box(c(33, 40, 29, 37.5), label = "Levant"),
#'   subregions = list(galilee = focus_box(c(35.05, 35.75, 32.55, 33.30)))
#' )
#' ```
#'
#' Each named entry adds three columns. `<name>_cells` is worth reading
#' alongside the other two: a maximum over a handful of cells jumps around
#' between slices in a way the mean does not.
#'
#' A subregion must lie inside `focus`; one that does not yields `NA` with
#' `cells = 0` rather than an error, since a site drifting outside the study
#' region is a legitimate thing to discover.
#'
#' @param richness A `SpatRaster` from [richness_stack()].
#' @param time Time label carried into the output row. Numeric years stay
#'   numeric so that downstream ordering is chronological.
#' @param focus The focus used to build `richness`, for labelling.
#' @param subregions Optional named list of [focus_box()] objects. Each adds
#'   `<name>_max`, `<name>_mean` and `<name>_cells` columns. See details.
#' @return A one-row tibble.
#' @seealso [richness_stack()]
#' @export
richness_stats <- function(richness, time = NA, focus = NULL, subregions = NULL) {

  vals <- terra::values(richness, mat = FALSE, na.rm = TRUE)
  if (length(vals) == 0) {
    rc_abort("The richness surface has no non-missing cells.")
  }

  out <- tibble::tibble(
    time            = time,
    focus           = focus$label %||% NA_character_,
    cells           = length(vals),
    mean_richness   = mean(vals),
    median_richness = stats::median(vals),
    max_richness    = max(vals),
    variance        = stats::var(vals)
  )

  for (nm in names(subregions %||% list())) {
    sub <- subregions[[nm]]

    # terra aborts with "extents do not overlap" when a subregion falls outside
    # the focus. A site sitting outside the study region is a real thing to
    # discover, so report it as zero cells rather than killing the run.
    sub_r <- tryCatch(
      terra::mask(terra::crop(richness, focus_ext(sub)), focus_vect(sub)),
      error = function(e) NULL
    )
    if (is.null(sub_r)) {
      cli::cli_warn(c(
        "Subregion {.val {nm}} does not overlap the focus.",
        "i" = "Its columns are NA with {.code {nm}_cells = 0}."
      ))
    }
    sub_vals <- if (is.null(sub_r)) {
      numeric(0)
    } else {
      terra::values(sub_r, mat = FALSE, na.rm = TRUE)
    }
    out[[paste0(nm, "_max")]] <- if (length(sub_vals)) max(sub_vals) else NA_real_
    out[[paste0(nm, "_mean")]] <- if (length(sub_vals)) mean(sub_vals) else NA_real_
    out[[paste0(nm, "_cells")]] <- length(sub_vals)
  }

  out
}
