# ==============================================================================
# Climate sources
#
# A climate source answers one question: "give me these variables, at this
# time, over this extent". Where the rasters come from -- a directory of
# pre-downloaded GeoTIFFs, or a live pastclim query -- is the source's problem.
#
# Aggregation belongs to the source and happens exactly once, when the slices
# are prepared. The source pipeline aggregated again inside the projection
# step, so its Gaussian-averaging branch silently coarsened the hindcast grid
# by a further factor of 10 relative to the grid the model was fitted on.
# ==============================================================================

#' Recover the years a set of directory names encodes
#'
#' Extracts a trailing signed integer, then keeps only those names the label
#' template actually reproduces. The round-trip check is what stops sibling
#' directories that merely end in a number -- `present_1985`, `chelsa_1950` --
#' from being mistaken for time slices.
#'
#' @param dirs Character vector of directory base names.
#' @param slice_fmt The [sprintf()] template used to build slice names.
#' @return A numeric vector of years, possibly empty.
#' @noRd
parse_slice_times <- function(dirs, slice_fmt) {
  keep <- grepl("-?[0-9]+$", dirs)
  dirs <- dirs[keep]
  if (length(dirs) == 0) return(numeric(0))

  candidates <- suppressWarnings(
    as.numeric(regmatches(dirs, regexpr("-?[0-9]+$", dirs)))
  )
  ok <- !is.na(candidates) &
    vapply(candidates, function(x) !is.na(x) && is.finite(x), logical(1))
  candidates <- candidates[ok]
  dirs <- dirs[ok]
  if (length(candidates) == 0) return(numeric(0))

  round_trips <- sprintf(slice_fmt, candidates) == dirs
  candidates[round_trips]
}

new_climate_source <- function(type, ..., resolver, times) {
  structure(
    list(type = type, params = list(...), resolver = resolver, times = times),
    class = "richcast_climate"
  )
}

#' Climate source: a directory of pre-downloaded rasters
#'
#' Expects one subdirectory per time slice, each holding one GeoTIFF per
#' variable named after that variable (`bio01.tif`, `bio12.tif`, ...). This is
#' the layout written by [prepare_climate()].
#'
#' ```
#' climate/
#'   present_1985/  bio01.tif bio04.tif ...
#'   time_0850/     bio01.tif bio04.tif ...
#'   time_1850/     bio01.tif bio04.tif ...
#' ```
#'
#' All slices must share a grid. This is checked on first use rather than
#' assumed, because a mismatch between the fitting grid and the projection
#' grid produces plausible-looking but meaningless suitability surfaces.
#'
#' # Choosing the fitting slice
#'
#' `present` names the directory a model is fitted against, and it need not be
#' an observational product. Fitting on WorldClim and projecting onto CHELSA
#' means the model crosses datasets between fitting and projection, so any
#' systematic difference between the two products is indistinguishable from
#' climate change. Where the palaeoclimate reconstruction publishes its own
#' present-day slice, fitting on that instead keeps one product throughout:
#'
#' ```r
#' climate_dir("climate", present = "chelsa_1950")
#' ```
#'
#' The trade-off is that a reconstruction's present-day slice is itself
#' modelled, so it inherits the model's biases rather than an observational
#' network's. Neither choice is free; running both and comparing is cheap once
#' the slices are on disk.
#'
#' # BCE slices
#'
#' Years are signed. `sprintf("time_%04d", -1050)` gives `time_-1050`, keeping
#' 1050 BCE in its own directory rather than colliding with 1050 CE.
#'
#' @param path Directory holding the time-slice subdirectories.
#' @param present Subdirectory name for the slice models are fitted against.
#' @param slice_fmt [sprintf()] template mapping a year CE to a subdirectory
#'   name. The default zero-pads to four digits, matching [prepare_climate()].
#' @return A `richcast_climate`.
#' @examples
#' clim <- climate_dir("Inputs/Climate/eurasia_slices")
#' clim
#' @export
climate_dir <- function(path,
                        present = "present_1985",
                        slice_fmt = "time_%04d") {

  # No abs(): sprintf("time_%04d", -1050) gives "time_-1050", which is how BCE
  # slices are named on disk. Taking the absolute value would send a request
  # for 1050 BCE to the 1050 CE directory and return the wrong climate without
  # any error.
  label_for <- function(time) {
    if (identical(time, "present")) present else sprintf(slice_fmt, time)
  }

  available <- if (dir.exists(path)) {
    dirs <- setdiff(basename(list.dirs(path, recursive = FALSE)), present)
    sort(parse_slice_times(dirs, slice_fmt))
  } else {
    numeric(0)
  }

  new_climate_source(
    "climate_dir",
    path = path, present = present, slice_fmt = slice_fmt,
    times = available[!is.na(available)],
    resolver = function(params, time, vars, extent) {
      dir_path <- file.path(params$path, label_for(time))
      if (!dir.exists(dir_path)) {
        rc_abort(c(
          "No climate rasters for time {.val {time}}.",
          "x" = "Expected directory {.path {dir_path}}.",
          "i" = "Prepare slices with {.fn prepare_climate}."
        ))
      }
      files <- list.files(dir_path, pattern = "\\.tif$", full.names = TRUE)
      if (length(files) == 0) {
        rc_abort("Directory {.path {dir_path}} contains no GeoTIFFs.")
      }
      r <- terra::rast(files)
      names(r) <- tools::file_path_sans_ext(basename(files))

      if (!is.null(vars)) {
        missing <- setdiff(vars, names(r))
        if (length(missing) > 0) {
          # qty() pins pluralisation to `missing`; without it cli sees two
          # vectors in the string and refuses to choose.
          rc_abort(c(
            "{cli::qty(length(missing))}Variable{?s} {.val {missing}} missing from time {.val {time}}.",
            "i" = "Available: {.val {names(r)}}."
          ))
        }
        r <- r[[vars]]
      }
      if (!is.null(extent)) r <- terra::crop(r, extent)
      r
    }
  )
}

#' Climate source: live pastclim queries
#'
#' Retrieves slices through \pkg{pastclim} on demand. Convenient for
#' exploration; for a full run over many species and time slices, prepare the
#' slices once with [prepare_climate()] and use [climate_dir()] instead, since
#' every species would otherwise re-stream the same rasters.
#'
#' @param dataset_present,dataset_past pastclim dataset names.
#' @param present_time Year CE treated as "present" when fitting.
#' @param times Numeric vector of years CE this source can serve.
#' @return A `richcast_climate`.
#' @references Leonardi, M., Hallet, E. Y., Beyer, R., Krapp, M., & Manica,
#'   A. (2023). pastclim 1.2: an R package to easily access and use
#'   paleoclimatic reconstructions. *Ecography*, 2023, e06481.
#'   \doi{10.1111/ecog.06481}
#'
#'   Cite the reconstruction you use as well, as \pkg{pastclim} asks:
#'
#'   Beyer, R. M., Krapp, M., & Manica, A. (2020). High-resolution terrestrial
#'   climate, bioclimate and vegetation for the last 120,000 years.
#'   *Scientific Data*, 7, 236. \doi{10.1038/s41597-020-0552-1}
#'
#'   Karger, D. N., Nobis, M. P., Normand, S., Graham, C. H., & Zimmermann,
#'   N. E. (2023). CHELSA-TraCE21k -- high-resolution (1 km) downscaled
#'   transient temperature and precipitation data since the Last Glacial
#'   Maximum. *Climate of the Past*, 19, 439-456.
#'   \doi{10.5194/cp-19-439-2023}
#'
#'   Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1-km spatial
#'   resolution climate surfaces for global land areas. *International
#'   Journal of Climatology*, 37(12), 4302-4315. \doi{10.1002/joc.5086}
#' @examples
#' clim <- pastclim_climate(times = seq(850, 1850, by = 100))
#' clim
#' @export
pastclim_climate <- function(dataset_present = "WorldClim_2.1_5m",
                             dataset_past = "CHELSA_trace21k_1.0_0.5m_vsi",
                             present_time = 1985,
                             times = NULL) {
  new_climate_source(
    "pastclim",
    dataset_present = dataset_present, dataset_past = dataset_past,
    present_time = present_time,
    times = times %||% numeric(0),
    resolver = function(params, time, vars, extent) {
      rlang::check_installed("pastclim", "to query climate reconstructions.")
      is_present <- identical(time, "present")
      ds <- if (is_present) params$dataset_present else params$dataset_past
      tt <- if (is_present) params$present_time else time
      if (is.null(vars)) {
        rc_abort("{.arg vars} must be given for a pastclim source.")
      }
      r <- pastclim::region_slice(
        time_ce = tt, bio_variables = vars, dataset = ds,
        ext = if (is.null(extent)) NULL else terra::ext(extent)
      )
      if (!is.null(extent)) r <- terra::crop(r, extent)
      r
    }
  )
}

#' @export
print.richcast_climate <- function(x, ...) {
  cli::cli_text("{.cls richcast_climate} <{x$type}>")
  if (length(x$times)) {
    cli::cli_text("  {length(x$times)} slice{?s}: {min(x$times)} to {max(x$times)} CE")
  } else {
    cli::cli_text("  slices: {.emph none detected}")
  }
  invisible(x)
}

#' Fetch climate for one time slice
#'
#' @param climate A `richcast_climate`.
#' @param time Year CE, or the string `"present"`.
#' @param vars Character vector of variable names, or `NULL` for all.
#' @param extent Optional `SpatExtent` to crop to.
#' @return A `SpatRaster`.
#' @noRd
climate_at <- function(climate, time, vars = NULL, extent = NULL) {
  if (!inherits(climate, "richcast_climate")) {
    rc_abort(c(
      "{.arg climate} must be a climate source.",
      "i" = "Build one with {.fn climate_dir} or {.fn pastclim_climate}."
    ))
  }
  climate$resolver(climate$params, time, vars, extent)
}
