# ==============================================================================
# Preparing climate slices
#
# Streaming palaeoclimate is the slowest part of any run and the part most
# likely to exhaust memory, so it is done once, up front, one variable at a
# time, written straight to disk and resumable after an interruption.
#
# Aggregation happens HERE and nowhere else. Doing it again downstream -- as
# the source pipeline's Gaussian branch did -- coarsens the projection grid
# relative to the grid the model was fitted on, which produces a suitability
# surface that looks fine and means nothing.
# ==============================================================================

#' Download and prepare climate slices for hindcasting
#'
#' Streams climate rasters through \pkg{pastclim}, aggregates them to a common
#' grid, and writes one GeoTIFF per variable per time slice in the layout
#' [climate_dir()] expects.
#'
#' By default the present-day slice is drawn from the same product as the
#' palaeoclimate slices. Mixing them -- fitting on WorldClim while projecting
#' onto CHELSA, as is easy to do by accident -- folds the step between the two
#' products into every `delta_from_present`, and that step is species-specific
#' in sign, so it does not cancel. Overriding `dataset_present` is supported
#' (an observational present has its own arguments in its favour) but it is now
#' a deliberate act, and [check_climate_products()] will say so afterwards.
#'
#' Where the two products differ in native resolution, set `agg_present` and
#' `agg_past` so both land on the same cell size; `check_grids = TRUE` verifies
#' that they did.
#'
#' Existing files are skipped, so an interrupted run can simply be restarted.
#'
#' @param path Output directory. One subdirectory per slice is created.
#' @param vars Bioclimatic variables to prepare: the eight [bioclim_vars] by
#'   default, or a subset of them.
#' @param times Numeric vector of years CE for the palaeoclimate slices.
#' @param extent Numeric `c(xmin, xmax, ymin, ymax)` to clip to. Clipping at
#'   this stage is what keeps the whole thing tractable.
#' @param present_time Year CE for the present-day fitting slice. Defaults to
#'   1950, which is 0 BP and the present-day slice most reconstructions publish.
#' @param dataset_past pastclim dataset for the palaeoclimate slices.
#' @param dataset_present pastclim dataset for the present-day slice. Defaults
#'   to `dataset_past`, so one product is used throughout unless you
#'   deliberately ask for two. See details.
#' @param agg_past,agg_present Aggregation factors; `agg_present` defaults to
#'   `agg_past` for the same reason.
#' @param present_dir Subdirectory name for the present-day slice.
#' @param slice_fmt [sprintf()] template for palaeoclimate subdirectory names.
#' @param check_grids Verify that every written slice shares one grid.
#' @param quiet Suppress progress messages.
#' @return A [climate_dir()] source pointing at `path`, invisibly.
#' @seealso [climate_dir()]
#' @examples
#' \dontrun{
#' prepare_climate(
#'   path   = "climate/eurasia",
#'   times  = seq(850, 1850, by = 100),
#'   extent = c(-15, 180, 10, 82)
#' )
#' }
#' @export
prepare_climate <- function(path,
                            vars = bioclim_vars,
                            times,
                            extent,
                            present_time = 1950,
                            dataset_past = "CHELSA_trace21k_1.0_0.5m_vsi",
                            dataset_present = dataset_past,
                            agg_past = 20,
                            agg_present = agg_past,
                            present_dir = "present",
                            slice_fmt = "time_%04d",
                            check_grids = TRUE,
                            quiet = FALSE) {

  rlang::check_installed("pastclim", "to download climate reconstructions.")
  vars <- check_predictors(vars, arg = "vars")
  say <- function(...) if (!quiet) cli::cli_alert_info(...)

  check_extent(extent)
  ext_obj <- terra::ext(extent[1], extent[2], extent[3], extent[4])

  # `path` is where slices are WRITTEN. Pointing it at the source NetCDF is an
  # easy confusion, and dir.create() fails silently on an existing file, so the
  # run limps on until writeRaster reports "cannot write file" several steps
  # later with no clue as to why.
  if (file.exists(path) && !dir.exists(path)) {
    rc_abort(c(
      "{.arg path} must be a directory to write slices into.",
      "x" = "{.path {path}} is an existing file.",
      "i" = "Source data is located with {.code pastclim::set_data_path()}; {.arg path} is only for richcast's output.",
      "i" = "Try something like {.path {file.path(dirname(path), 'beyer_slices')}}."
    ))
  }
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) {
      rc_abort("Could not create output directory {.path {path}}.")
    }
  }

  ensure_dataset(dataset_present, quiet)
  ensure_dataset(dataset_past, quiet)

  jobs <- c(
    list(list(dir = present_dir, time = present_time,
              dataset = dataset_present, agg = agg_present)),
    # No abs(): a signed year keeps BCE and CE slices in separate directories.
    # With abs(), -850 and 850 would collide on one directory and whichever ran
    # second would be skipped as already present.
    lapply(sort(unique(times)), function(tt) {
      list(dir = sprintf(slice_fmt, tt), time = tt,
           dataset = dataset_past, agg = agg_past)
    })
  )

  for (job in jobs) {
    out_dir <- file.path(path, job$dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    done <- length(list.files(out_dir, pattern = "\\.tif$"))
    if (done >= length(vars)) {
      say("{job$dir}: complete ({done}/{length(vars)}), skipping.")
      next
    }
    say("{job$dir}: {done}/{length(vars)} present, fetching the rest.")

    for (v in vars) {
      fname <- file.path(out_dir, paste0(v, ".tif"))
      if (file.exists(fname)) next

      r <- tryCatch(
        pastclim::region_slice(time_ce = job$time, bio_variables = v,
                               dataset = job$dataset, ext = ext_obj),
        error = function(e) {
          cli::cli_warn("{job$dir}/{v}: {conditionMessage(e)}")
          NULL
        }
      )
      if (is.null(r)) next

      # Datasets already at the target grid need no aggregation. terra warns
      # "nothing to do" for fact = 1, once per variable per slice, which for a
      # 0.5-degree reconstruction like Beyer2020 is pure noise.
      if (job$agg > 1) {
        r <- terra::aggregate(r, fact = job$agg, fun = mean, na.rm = TRUE)
      }
      terra::writeRaster(r, fname, overwrite = TRUE,
                         wopt = list(gdal = "COMPRESS=LZW"))
      rm(r)
      gc()
      say("  wrote {job$dir}/{v}.tif")
    }
  }

  write_manifest(path, jobs, vars)

  src <- climate_dir(path, present = present_dir, slice_fmt = slice_fmt)
  if (check_grids) check_climate_grids(src, vars, quiet = quiet)
  check_climate_products(src, quiet = quiet)
  invisible(src)
}

#' Record which dataset each slice came from
#'
#' Without this richcast cannot tell a CHELSA slice from a WorldClim one --
#' they are both just GeoTIFFs on a matching grid -- and a mixed pipeline
#' passes every check while quietly making `delta_from_present` meaningless.
#' A manifest turns that into something detectable.
#'
#' @param path Slice directory.
#' @param jobs The per-slice job list built by [prepare_climate()].
#' @param vars Variables requested.
#' @return Invisibly, the manifest data frame.
#' @noRd
write_manifest <- function(path, jobs, vars) {
  man <- data.frame(
    slice   = vapply(jobs, function(j) j$dir, character(1)),
    time_ce = vapply(jobs, function(j) as.numeric(j$time), numeric(1)),
    dataset = vapply(jobs, function(j) j$dataset, character(1)),
    aggregation = vapply(jobs, function(j) as.numeric(j$agg), numeric(1)),
    stringsAsFactors = FALSE
  )
  man$variables <- paste(vars, collapse = ",")
  utils::write.csv(man, file.path(path, "richcast_manifest.csv"),
                   row.names = FALSE)
  invisible(man)
}

#' Warn when the present slice and the time slices come from different products
#'
#' Fitting on one climate product and projecting onto another folds the step
#' between them into every `delta_from_present`. Measured on one species, that
#' step was 42% of its present-day range and flipped it from above-present in
#' 2 of 11 centuries to above-present in all 11 -- and the offset is
#' species-specific in sign, so it does not cancel.
#'
#' Only works where the slices were prepared by [prepare_climate()], which
#' leaves a manifest. Directories assembled by hand are unreadable in this
#' respect and pass silently.
#'
#' @param climate A [climate_dir()] source.
#' @param quiet Suppress the all-clear message.
#' @return `TRUE` if consistent (or undeterminable), `FALSE` if mixed.
#' @export
check_climate_products <- function(climate, quiet = FALSE) {

  if (!identical(climate$type, "climate_dir")) {
    rc_abort("{.fn check_climate_products} only applies to a {.fn climate_dir} source.")
  }
  man_path <- file.path(climate$params$path, "richcast_manifest.csv")
  if (!file.exists(man_path)) {
    if (!quiet) {
      cli::cli_alert_info(
        "No manifest in {.path {climate$params$path}}; cannot verify that the slices share one climate product."
      )
    }
    return(invisible(TRUE))
  }

  man <- utils::read.csv(man_path, stringsAsFactors = FALSE)
  present_ds <- man$dataset[man$slice == climate$params$present]
  past_ds    <- unique(man$dataset[man$slice != climate$params$present])

  if (length(present_ds) == 0 || length(past_ds) == 0) {
    return(invisible(TRUE))
  }
  if (!all(past_ds == present_ds)) {
    cli::cli_warn(c(
      "Present-day and palaeoclimate slices come from different products.",
      "x" = "present: {.val {present_ds}}; slices: {.val {past_ds}}.",
      "i" = "Every {.code delta_from_present} then contains the step between the two products, which is species-specific in sign and does not cancel.",
      "i" = "Prepare the present slice from {.val {past_ds}} so that fitting and projection share one product."
    ))
    return(invisible(FALSE))
  }
  if (!quiet) {
    cli::cli_alert_success("All slices come from {.val {present_ds}}.")
  }
  invisible(TRUE)
}

#' Verify that every prepared slice shares one grid
#'
#' A silent resample between the fitting grid and the projection grid is one of
#' the easier ways to get confident nonsense out of an SDM, so this compares
#' geometry with zero tolerance and fails loudly.
#'
#' @param climate A [climate_dir()] source.
#' @param vars Variables that must be present in every slice.
#' @param quiet Suppress the success message.
#' @return `TRUE`, invisibly, or an error.
#' @export
check_climate_grids <- function(climate, vars = NULL, quiet = FALSE) {

  if (!identical(climate$type, "climate_dir")) {
    rc_abort("{.fn check_climate_grids} only applies to a {.fn climate_dir} source.")
  }
  path <- climate$params$path
  dirs <- list.dirs(path, recursive = FALSE)
  if (length(dirs) == 0) rc_abort("No slice directories in {.path {path}}.")

  # Only inspect directories that are actually slices. Download tooling leaves
  # scratch directories (".staging", "tmp") alongside the data, and a folder
  # the naming scheme does not recognise is not a slice with a problem -- it is
  # not a slice.
  recognised <- c(
    climate$params$present,
    sprintf(climate$params$slice_fmt,
            parse_slice_times(basename(dirs), climate$params$slice_fmt))
  )
  skipped <- setdiff(basename(dirs), recognised)
  dirs <- dirs[basename(dirs) %in% recognised]

  if (length(skipped) > 0 && !quiet) {
    cli::cli_alert_info("Ignoring {length(skipped)} non-slice director{?y/ies}: {.val {skipped}}.")
  }
  if (length(dirs) == 0) {
    rc_abort(c(
      "No recognised slice directories in {.path {path}}.",
      "i" = "Expected names matching {.val {climate$params$slice_fmt}} or {.val {climate$params$present}}."
    ))
  }

  ref <- NULL
  ref_name <- NULL
  problems <- character(0)

  for (d in dirs) {
    files <- list.files(d, pattern = "\\.tif$", full.names = TRUE)
    if (length(files) == 0) {
      problems <- c(problems, paste0(basename(d), ": no rasters"))
      next
    }
    # Variables within a slice can disagree with each other, not just across
    # slices, and terra refuses to stack those at all. Report it as an
    # inconsistency rather than letting the read abort the check.
    r <- tryCatch(terra::rast(files), error = function(e) NULL)
    if (is.null(r)) {
      problems <- c(problems, paste0(
        basename(d), ": variables within the slice have mismatched grids"
      ))
      next
    }
    names(r) <- tools::file_path_sans_ext(basename(files))

    if (!is.null(vars)) {
      missing <- setdiff(vars, names(r))
      if (length(missing) > 0) {
        problems <- c(problems,
                      paste0(basename(d), ": missing ", paste(missing, collapse = ", ")))
      }
    }
    if (is.null(ref)) {
      ref <- r
      ref_name <- basename(d)
    } else if (!terra::compareGeom(ref, r, stopOnError = FALSE, messages = FALSE)) {
      problems <- c(problems, sprintf(
        "%s: grid differs from %s (res %s vs %s)",
        basename(d), ref_name,
        paste(round(terra::res(r), 5), collapse = "x"),
        paste(round(terra::res(ref), 5), collapse = "x")
      ))
    }
  }

  if (length(problems) > 0) {
    rc_abort(c(
      "Climate slices are not consistent:",
      stats::setNames(problems, rep("x", length(problems))),
      "i" = "Delete the offending slices and re-run {.fn prepare_climate}."
    ))
  }
  if (!quiet) {
    cli::cli_alert_success(
      "{length(dirs)} slices share one grid ({paste(round(terra::res(ref), 4), collapse = ' x ')} deg)."
    )
  }
  invisible(TRUE)
}

#' @noRd
ensure_dataset <- function(dataset, quiet = FALSE) {
  have <- pastclim::get_downloaded_datasets()
  if (is.null(have[[dataset]])) {
    if (!quiet) cli::cli_alert_info("Downloading {.val {dataset}} (this is slow).")
    pastclim::download_dataset(dataset)
  }
  invisible(TRUE)
}
