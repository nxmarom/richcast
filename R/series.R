# ==============================================================================
# The whole pipeline, over a series of time slices
# ==============================================================================

#' Hindcast a whole assemblage across a series of time slices
#'
#' Fits one model per species, projects each onto every requested time slice,
#' and summarises assemblage richness at each. This is the top-level entry
#' point; the pieces ([fit_sdm()], [project_sdm()], [richness_stack()],
#' [richness_stats()]) remain usable on their own.
#'
#' Each species is fitted **once**, then projected repeatedly. Fitting inside
#' the time loop -- as the pipeline this package was extracted from did --
#' repeats identical work once per slice, since the model depends only on
#' present-day climate.
#'
#' @section Replicate intervals:
#'
#' A single fit is one draw. Redrawing the pseudo-occurrences moves the answer
#' a great deal: for one Tian Shan marmot, six seeds at the default
#' `nsample = 100` spanned a 2.4-fold range in present-day extent and
#' disagreed on the *sign* of the change at one slice. With `replicates > 1`,
#' `res$species` gains `cells_min`, `cells_max` and `cells_sd`, `cells` becomes
#' the replicate median, and `res$richness` gains `mean_richness_lo/hi/sd`.
#'
#' Richness intervals are built by stacking replicate *r* of every species into
#' richness surface *r*, then taking quantiles across replicates -- not by
#' combining per-species intervals, which would be wrong for a sum.
#'
#' A wide replicate interval is not always the right output. Where a species
#' occupies a handful of grid cells, whether its suitability clears the
#' threshold at all depends on which cells the sample happened to hit, and the
#' estimate is an artefact of the draw rather than a measurement with error
#' around it. Widening the interval implies the quantity exists and is merely
#' imprecise. A future release will omit such species by default, reporting
#' them as `not resolvable`; for now, screen on `present_cells` and on the
#' `cells_sd` column these replicates produce.
#'
#' Read these as **precision, not accuracy**. They describe how much the answer
#' moves when the sample is redrawn from the same range polygon under the same
#' reconstruction. They exclude the range polygon being wrong, the
#' reconstruction being wrong, niche conservatism, and the choice of model --
#' and in this pipeline the climate-product difference has been *larger* than
#' the replicate spread. Use [ensemble_series()] for that layer.
#'
#' @section Choosing a baseline:
#'
#' `delta_from_present` subtracts a baseline range size from each slice, so it
#' only means anything if the baseline and the slices are commensurable. They
#' are not when the present-day slice comes from a different climate product
#' than the palaeoclimate series -- fitting on WorldClim and projecting onto
#' CHELSA, say. The product step then enters every delta, and it does not
#' cancel: for one species it inflated the present-day range by 42%, flipping
#' that species from "above present in 2 of 11 centuries" to "above present in
#' all 11". The offset is species-specific in both size and sign, so it
#' distorts comparisons between taxa as well as within them.
#'
#' Two ways to keep the comparison honest:
#'
#' * Use a climate source that is one product throughout, including its
#'   present-day slice -- `climate_dir(path, present = "chelsa_1950")` rather
#'   than a WorldClim present beside CHELSA slices. Then the default
#'   `baseline = "present"` is already consistent.
#' * Where the reconstruction publishes no present-day slice, set `baseline`
#'   to its youngest time slice and describe results as change relative to
#'   that, not to the present.
#'
#' richcast cannot detect which product a directory of GeoTIFFs came from, so
#' it cannot warn you automatically. The choice is yours to make explicitly.
#'
#' @param db A `richcast_db` from [build_taxon_db()].
#' @param climate A climate source from [climate_dir()] or [pastclim_climate()].
#' @param times Numeric vector of years CE.
#' @param focus Region to summarise richness over, from [focus_box()] etc.
#' @param subregions Optional named list of sub-foci, passed to
#'   [richness_stats()].
#' @param species Optional subset of species names. Defaults to every species
#'   in `db` whose range intersects `focus` (expanded by `prefilter_buffer`).
#' @param prefilter_buffer Degrees to expand `focus` by when deciding which
#'   species are relevant. Should be at least the maximum fitting buffer.
#' @param window Optional [gaussian_window()] for time-averaging.
#' @param resolution Richness grid resolution in degrees.
#' @param replicates Number of replicate fits per species. Each redraws the
#'   presence and background samples and re-optimises the threshold, so the
#'   spread across replicates can be reported. `1` (default) fits once, exactly
#'   as before. See details.
#' @param conf Interval width for replicate summaries. `0.9` gives the 5th and
#'   95th percentiles across replicates.
#' @param baseline Slice that `delta_from_present` is measured against, and
#'   that the `present` richness row is built from. `"present"` (default) uses
#'   the slice the models were fitted on; a numeric year uses that slice
#'   instead. Either way the baseline is projected through the same path as
#'   the hindcast slices. See details.
#' @param keep_surfaces Retain the richness raster for every slice, so maps can
#'   be drawn afterwards. Stored wrapped, so the result still survives
#'   `saveRDS()`. Set `FALSE` for very large foci where only the summary
#'   statistics are wanted.
#' @param on_error `"warn"` skips a failing species and carries on;
#'   `"stop"` aborts the run.
#' @param quiet Suppress per-species progress.
#' @param ... Further arguments passed to [fit_sdm()].
#' @return A `richcast_series` with three tibbles -- `richness` (one row per
#'   time slice), `species` (one row per species per slice, with deltas), and
#'   `models` (fit diagnostics) -- plus the fitted models and projected ranges.
#' @seealso [fit_sdm()], [richness_stats()]
#' @export
run_hindcast_series <- function(db,
                                climate,
                                times,
                                focus,
                                subregions = NULL,
                                species = NULL,
                                prefilter_buffer = 8,
                                window = NULL,
                                resolution = 0.1,
                                baseline = "present",
                                replicates = 1,
                                conf = 0.9,
                                keep_surfaces = TRUE,
                                on_error = c("warn", "stop"),
                                quiet = FALSE,
                                ...) {

  on_error <- match.arg(on_error)
  if (!is.numeric(times) || length(times) == 0) {
    rc_abort("{.arg times} must be a non-empty numeric vector of years CE.")
  }
  # Chronological from the outset. The source pipeline sorted slice labels as
  # character, which puts 850 and 950 after 1850 and mis-lags every delta.
  times <- sort(unique(as.numeric(times)))

  targets <- select_species(db, focus, species, prefilter_buffer, quiet)

  models <- list()
  sets <- list()          # replicate sets, when replicates > 1
  ranges <- list()        # ranges[[time]][[species]] -- representative draw
  rep_ranges <- list()    # rep_ranges[[time]][[species]][[replicate]]
  baselines <- list()     # baseline range per species
  counts <- list()

  if (!is.numeric(replicates) || length(replicates) != 1 || replicates < 1) {
    rc_abort("{.arg replicates} must be a single positive number.")
  }
  replicates <- as.integer(replicates)
  for (tt in as.character(times)) {
    ranges[[tt]] <- list()
    rep_ranges[[tt]] <- list()
  }

  # A progress bar redraws on every update, which is useful at a console and
  # pure noise in a knitted document, so it follows `quiet` like everything else.
  if (!quiet) {
    cli::cli_progress_bar("Fitting and projecting", total = nrow(targets),
                          .envir = environment())
  }

  for (i in seq_len(nrow(targets))) {
    sp <- targets$species[i]
    if (!quiet) cli::cli_progress_update(.envir = environment())

    if (replicates > 1) {
      set <- try_step(
        fit_replicates(db, sp, climate, replicates = replicates,
                       quiet = TRUE, ...),
        what = "fit", species = sp, on_error = on_error
      )
      if (is.null(set)) next
      sets[[sp]] <- set
      fitted <- set$fits[[1]]
    } else {
      fitted <- try_step(
        fit_sdm(db, sp, climate, quiet = quiet, ...),
        what = "fit", species = sp, on_error = on_error
      )
      if (is.null(fitted)) next
    }
    models[[sp]] <- fitted

    # The baseline every delta is measured against. Projected through the same
    # path as the past slices, so the two are commensurable. When it is the
    # fitting slice this reproduces fitted$present_cells exactly, so the
    # default costs nothing.
    if (identical(baseline, "present")) {
      base_cells <- fitted$present_cells
      base_range <- fitted$present_range
    } else {
      base <- try_step(
        project_sdm(fitted, climate, baseline, window = window, quiet = quiet),
        what = paste("project baseline", baseline), species = sp,
        on_error = on_error
      )
      base_cells <- if (is.null(base)) NA_integer_ else base$cells
      base_range <- if (is.null(base)) NULL else base$range
    }
    baselines[[sp]] <- base_range

    for (tt in times) {
      if (replicates > 1) {
        ps <- try_step(
          project_replicates(sets[[sp]], climate, tt, window = window,
                             quiet = quiet),
          what = paste("project onto", tt), species = sp, on_error = on_error
        )
        if (is.null(ps)) next
        rep_ranges[[as.character(tt)]][[sp]] <- ps$ranges
        # The representative range is the replicate whose extent is the median,
        # so maps show a typical draw rather than an arbitrary one.
        med <- which.min(abs(ps$cells - stats::median(ps$cells)))[1]
        ranges[[as.character(tt)]][[sp]] <- ps$ranges[[med]]
        cell_val <- stats::median(ps$cells)
        counts[[length(counts) + 1L]] <- tibble::tibble(
          species = sp, time = tt, cells = cell_val,
          present_cells = base_cells,
          cells_min = min(ps$cells), cells_max = max(ps$cells),
          cells_sd = stats::sd(ps$cells)
        )
      } else {
        proj <- try_step(
          project_sdm(fitted, climate, tt, window = window, quiet = quiet),
          what = paste("project onto", tt), species = sp, on_error = on_error
        )
        if (is.null(proj)) next
        ranges[[as.character(tt)]][[sp]] <- proj$range
        counts[[length(counts) + 1L]] <- tibble::tibble(
          species = sp, time = tt, cells = proj$cells,
          present_cells = base_cells
        )
      }
    }
  }
  if (!quiet) cli::cli_progress_done(.envir = environment())

  if (length(models) == 0) {
    rc_abort("No species could be modelled; nothing to summarise.")
  }

  # --- Richness per slice -------------------------------------------------
  # The surfaces are kept, not just their summaries: a mean over a region is a
  # poor substitute for seeing where in that region the species actually are.
  surfaces <- list()
  richness <- lapply(times, function(tt) {
    slice_ranges <- ranges[[as.character(tt)]]

    # A slice where every projection failed produces an all-zero surface, which
    # is indistinguishable from a genuine absence of species. Say so: a missing
    # climate slice reads as "richness collapsed to nothing" otherwise.
    if (length(slice_ranges) == 0) {
      cli::cli_warn(c(
        "No species contributed a range at {tt} CE.",
        "x" = "Reported richness for this slice is zero because nothing was projected, not because the assemblage was empty.",
        "i" = "Check that climate rasters exist for {tt} CE."
      ))
    }

    r <- richness_stack(slice_ranges, focus,
                        resolution = resolution, quiet = TRUE)
    if (keep_surfaces) surfaces[[as.character(tt)]] <<- terra::wrap(r)
    out <- richness_stats(r, time = tt, focus = focus, subregions = subregions)
    out$species_contributing <- length(slice_ranges)
    out
  })
  richness <- dplyr::bind_rows(richness)

  # --- Richness intervals across replicates -------------------------------
  # Richness is a sum over species, so per-taxon intervals cannot simply be
  # added. Replicate r of every species is stacked together to give richness
  # replicate r, and the spread is taken across those. Pairing independent
  # draws is arbitrary but valid as a Monte Carlo sample of the joint; adding
  # marginals is neither.
  if (replicates > 1) {
    lo_p <- (1 - conf) / 2
    hi_p <- 1 - lo_p
    reps <- lapply(seq_len(replicates), function(r) {
      vapply(times, function(tt) {
        key <- as.character(tt)
        per_sp <- lapply(rep_ranges[[key]], function(x) x[[r]])
        per_sp <- per_sp[!vapply(per_sp, is.null, logical(1))]
        rr <- richness_stack(per_sp, focus, resolution = resolution,
                             quiet = TRUE)
        mean(terra::values(rr, mat = FALSE, na.rm = TRUE))
      }, numeric(1))
    })
    m <- do.call(rbind, reps)              # replicates x times
    richness$mean_richness_lo <- apply(m, 2, stats::quantile, probs = lo_p,
                                       na.rm = TRUE)
    richness$mean_richness_hi <- apply(m, 2, stats::quantile, probs = hi_p,
                                       na.rm = TRUE)
    richness$mean_richness_sd <- apply(m, 2, stats::sd, na.rm = TRUE)
  }

  # --- Present-day baseline ----------------------------------------------
  present_ranges <- baselines[names(models)]
  present_r <- richness_stack(present_ranges, focus,
                              resolution = resolution, quiet = TRUE)
  if (keep_surfaces) surfaces[["present"]] <- terra::wrap(present_r)
  present_row <- richness_stats(present_r, time = NA_real_, focus = focus,
                                subregions = subregions)
  present_row$species_contributing <- sum(
    !vapply(present_ranges, is.null, logical(1))
  )
  if (replicates > 1) {
    present_row$mean_richness_lo <- NA_real_
    present_row$mean_richness_hi <- NA_real_
    present_row$mean_richness_sd <- NA_real_
  }
  present_row$period <- "present"
  richness$period <- paste0("hindcast_", richness$time)
  richness <- dplyr::bind_rows(present_row, richness)

  # --- Species-level deltas ----------------------------------------------
  species_tbl <- dplyr::bind_rows(counts) |>
    dplyr::arrange(.data$species, .data$time) |>   # numeric time: correct lag
    dplyr::group_by(.data$species) |>
    dplyr::mutate(
      delta_from_previous = .data$cells - dplyr::lag(.data$cells),
      delta_from_present  = .data$cells - .data$present_cells
    ) |>
    dplyr::ungroup()

  models_tbl <- tibble::tibble(
    species       = names(models),
    auc           = vapply(models, function(m) m$auc, numeric(1)),
    threshold     = vapply(models, function(m) m$threshold, numeric(1)),
    n_presence    = vapply(models, function(m) m$n_presence, integer(1)),
    n_background  = vapply(models, function(m) m$n_background, integer(1)),
    present_cells = vapply(models, function(m) m$present_cells, integer(1))
  )

  structure(
    list(
      richness = richness,
      species  = species_tbl,
      models   = models_tbl,
      fits     = models,
      ranges   = ranges,
      surfaces = surfaces,
      times    = times,
      focus    = focus,
      baseline = baseline,
      replicates = replicates,
      sets     = sets
    ),
    class = "richcast_series"
  )
}

#' Extract a richness surface from a series
#'
#' @param series A `richcast_series` from [run_hindcast_series()].
#' @param time Year CE, or `"present"`.
#' @return A `SpatRaster` of species counts.
#' @seealso [richness_grid()], [run_hindcast_series()]
#' @export
richness_surface <- function(series, time) {
  if (!inherits(series, "richcast_series")) {
    rc_abort("{.arg series} must come from {.fn run_hindcast_series}.")
  }
  key <- as.character(time)
  if (length(series$surfaces) == 0) {
    rc_abort(c(
      "This series holds no richness surfaces.",
      "i" = "Re-run with {.code keep_surfaces = TRUE}."
    ))
  }
  if (!key %in% names(series$surfaces)) {
    rc_abort(c(
      "No surface for time {.val {time}}.",
      "i" = "Available: {.val {names(series$surfaces)}}."
    ))
  }
  terra::unwrap(series$surfaces[[key]])
}

#' Richness surfaces as a tidy data frame
#'
#' Long-format cell values for every slice, ready for a faceted heatmap. One
#' row per cell per time slice, so a focus of a few tens of thousands of cells
#' across a dozen slices is comfortable and a continental one is not -- crop
#' the focus, or select `times`, before reaching for this on a large region.
#'
#' @param series A `richcast_series`.
#' @param times Which slices to include. Defaults to every hindcast slice;
#'   pass `"present"` to include the present-day baseline.
#' @param drop_na Drop cells outside the focus mask.
#' @return A tibble with `x`, `y`, `time` and `richness`.
#' @seealso [richness_surface()]
#' @examples
#' \dontrun{
#' grid <- richness_grid(res)
#' ggplot2::ggplot(grid, ggplot2::aes(x, y, fill = richness)) +
#'   ggplot2::geom_raster() +
#'   ggplot2::facet_wrap(~time)
#' }
#' @export
richness_grid <- function(series, times = NULL, drop_na = TRUE) {
  if (!inherits(series, "richcast_series")) {
    rc_abort("{.arg series} must come from {.fn run_hindcast_series}.")
  }
  if (length(series$surfaces) == 0) {
    rc_abort(c(
      "This series holds no richness surfaces.",
      "i" = "Re-run with {.code keep_surfaces = TRUE}."
    ))
  }
  keys <- as.character(times %||% series$times)
  missing <- setdiff(keys, names(series$surfaces))
  if (length(missing) > 0) {
    rc_abort("No surface for {.val {missing}}.")
  }

  out <- lapply(keys, function(k) {
    r <- terra::unwrap(series$surfaces[[k]])
    d <- terra::as.data.frame(r, xy = TRUE, na.rm = drop_na)
    names(d)[3] <- "richness"
    d$time <- k
    tibble::as_tibble(d)
  })
  out <- dplyr::bind_rows(out)
  # Keep chronology numeric where possible; "present" forces character.
  if (!any(keys == "present")) out$time <- as.numeric(out$time)
  out[, c("x", "y", "time", "richness")]
}

#' @export
print.richcast_series <- function(x, ...) {
  cli::cli_text("{.cls richcast_series}")
  cli::cli_text("  {length(x$fits)} species x {length(x$times)} slices ({min(x$times)}-{max(x$times)} CE)")
  cli::cli_text("  focus: {x$focus$label}")
  cli::cli_text("  {.code $richness} {.code $species} {.code $models} {.code $fits} {.code $ranges}")
  invisible(x)
}


# ==============================================================================
# Internals
# ==============================================================================

#' Choose which species are worth modelling for a given focus
#' @noRd
select_species <- function(db, focus, species, buffer, quiet) {
  if (!is.null(species)) {
    keep <- db$species %in% normalise_species(species)
    missing <- setdiff(normalise_species(species), db$species)
    if (length(missing) > 0) {
      cli::cli_warn("Not in database, skipped: {.val {missing}}.")
    }
    return(sf::st_drop_geometry(db[keep, "species"]))
  }

  bb <- sf::st_bbox(focus$geometry)
  wide <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = bb[["xmin"]] - buffer, xmax = bb[["xmax"]] + buffer,
      ymin = bb[["ymin"]] - buffer, ymax = bb[["ymax"]] + buffer),
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
      "{sum(hits)}/{nrow(db)} species intersect the focus (buffered {buffer} deg)."
    )
  }
  if (sum(hits) == 0) {
    rc_abort(c(
      "No species ranges intersect the focus.",
      "i" = "Check the focus coordinates are c(xmin, xmax, ymin, ymax) in degrees."
    ))
  }
  sf::st_drop_geometry(db[hits, "species"])
}

#' Run one pipeline step, converting failure into a skip or an abort
#' @noRd
try_step <- function(expr, what, species, on_error) {
  tryCatch(
    expr,
    error = function(e) {
      msg <- c(
        "Failed to {what} for {.val {species}}.",
        "x" = conditionMessage(e)
      )
      if (on_error == "stop") cli::cli_abort(msg) else cli::cli_warn(msg)
      NULL
    }
  )
}
