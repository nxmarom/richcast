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
#' Do not read a *narrow* interval as a stable answer. Replicates redraw the
#' sample, so they measure sampling variability -- and a range small enough to
#' be sampled exhaustively has none. Sampling 1000 points with replacement from
#' an 8-cell range returns the same 8 cells on every seed, so the replicate
#' spread collapses toward zero exactly where the estimate is least
#' trustworthy. Measured on synthetic ranges clipped to a fixed number of
#' cells, an 8-cell range reproduced its projected extent identically across
#' five seeds at slices where a 512-cell range varied by 10%. The interval was
#' tight because there was nothing left to resample, not because the answer was
#' firm.
#'
#' Screen on `range_cells` via `min_cells` instead; see the next section.
#'
#' Read these as **precision, not accuracy**. They describe how much the answer
#' moves when the sample is redrawn from the same range polygon under the same
#' reconstruction. They exclude the range polygon being wrong, the
#' reconstruction being wrong, niche conservatism, and the choice of model --
#' and in this pipeline the climate-product difference has been *larger* than
#' the replicate spread. Use [ensemble_series()] for that layer.
#'
#' @section Resolvability:
#'
#' Some species are too small for the climate grid to say anything about. A
#' range covering `k` cells gives the model `k` distinct climate vectors --
#' presences are drawn with replacement, so `nsample` cannot manufacture more
#' -- and every projected extent is then a small integer whose value turns on
#' which handful of cells the range happens to contain.
#'
#' `min_cells` names the floor. Species below it are **still fitted, projected
#' and reported**: they keep their rows in `$species`, `$models` and `$ranges`,
#' so nothing becomes uninspectable. They are held out of the richness
#' surfaces only, and listed in `$resolvability` with their range size and the
#' reason. Silent shrinkage of the assemblage is the failure this is meant to
#' avoid, so the omission is always visible in the object and in `print()`.
#'
#' `range_cells` counts cells, not area, because the grid is what limits the
#' inference: the same range is better resolved on a finer reconstruction. A
#' cutoff in cells therefore only means something alongside the resolution it
#' was measured at, and `$resolvability` records the `min_cells` used.
#'
#' `min_cells` is also stored on the series itself, so a saved result carries
#' the setting it was produced under. That matters because a series written
#' before the screen existed, one written with the screen disabled, and one
#' where every species passed all print as "nothing held out" -- and only the
#' last is a clean bill of health. Its absence identifies the first case;
#' `print()` reports which of the three applies.
#'
#' The default of 100 is the largest cutoff at which every species measured
#' below it was demonstrably unstable: across 31 species spanning 11 to 12549
#' cells, all 8 below 100 varied by more than 25% across sampling seeds at
#' their worst slice, and at 200 the rule stops holding. `min_cells = 0`
#' disables the screen. `?richcast-resolvability` records the measurements,
#' including what the screen does *not* fix.
#'
#' @section Two geographies:
#'
#' `$species` reports every range size twice, over two different regions, and
#' which one a figure should use is a real choice rather than a formality.
#'
#' `cells` is counted over the species' **study extent** -- its own range plus
#' the fitting buffer -- because that is the region the model was fitted and
#' projected over. It is the honest measure of what happened to the *species*.
#' For a widespread taxon it is a continental number: in the Tian Shan
#' vignette *Mus musculus* contributes 190225 cells at 1350 CE, of which 13215
#' lie in the study region.
#'
#' `focus_cells` is the same range clipped to `focus`, on the same grid and
#' with the same `touches` rule the richness surface uses. It is the measure
#' that answers "what happened *here*", and it is the one commensurable with
#' `$richness`.
#'
#' The two are counts on **different grids**: `cells` on the climate grid,
#' `focus_cells` on the `resolution` grid the richness surface is built at.
#' Their ratio is therefore not a fraction of range inside the focus, and
#' `focus_cells` can exceed `cells` where the richness grid is the finer of
#' the two. Compare each column against itself across slices -- the grid
#' cancels -- and use area, not counts, to compare one against the other.
#'
#' The distinction matters most for the species `focus` barely contains.
#' Species are selected by intersecting a `prefilter_buffer`-degree expansion
#' of the focus, so a run legitimately includes taxa that reach toward the
#' study region without entering it. Those species have rows in `$species`
#' and trajectories that move -- and `focus_cells` of zero at every slice.
#' Plotting `cells` without checking `focus_cells` puts them on the page as
#' though they were part of the local assemblage.
#'
#' Both are differenced the same two ways, giving `focus_delta_from_previous`
#' and `focus_delta_from_present` alongside the study-extent pair.
#'
#' @section Is the focus inside the niche:
#'
#' `focus_cells` says how much of the focus a species holds. It does not say
#' whether holding it means anything, and that is a separate question with a
#' separate answer.
#'
#' A projected range is a threshold applied to a continuous surface. Where the
#' focus sits well inside a species' modelled niche, the cell count moves when
#' the climate moves. Where it sits *at* the cutoff, the count moves when
#' anything moves: a shift of a few hundredths in suitability flips large areas
#' across the threshold, and the resulting trajectory is a property of the
#' cutoff rather than of the region. The two cases are indistinguishable in
#' `focus_cells` alone -- both give plausible numbers that change between
#' slices.
#'
#' `focus_suit_q90` is the ninetieth percentile of suitability inside the focus
#' before thresholding, and `focus_suit_margin` is that minus the species' own
#' threshold. A negative margin means the focus does not clear the cutoff even
#' at its ninetieth percentile, so whatever cells the species contributes there
#' are an edge effect. Screen on `focus_suit_margin > 0` before reading
#' regional trajectories as range responses.
#'
#' Volatility is not a usable substitute. In the Tian Shan vignette nine of 26
#' species have a negative margin, holding between 0.4% and 17% of the region
#' above threshold; a fold-change screen on `focus_cells` catches only three of
#' them, because a footprint can stay small without swinging. The six it misses
#' include the species that led the regional growth ranking.
#'
#' With `replicates > 1`, `focus_cells` is measured on the representative
#' (median-extent) replicate -- the one `$ranges` stores -- rather than being
#' a median of per-replicate focus counts.
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
#' @param min_cells Smallest range, in grid cells, that can support a
#'   trajectory. Species below it are still fitted, projected and reported, but
#'   are left out of the richness surfaces and listed in `res$resolvability`.
#'   `0` disables the screen. See the resolvability section.
#' @param keep_surfaces Retain the richness raster for every slice, so maps can
#'   be drawn afterwards. Stored wrapped, so the result still survives
#'   `saveRDS()`. Set `FALSE` for very large foci where only the summary
#'   statistics are wanted.
#' @param on_error `"warn"` skips a failing species and carries on;
#'   `"stop"` aborts the run.
#' @param quiet Suppress per-species progress.
#' @param ... Further arguments passed to [fit_sdm()].
#' @return A `richcast_series` with three tibbles -- `richness` (one row per
#'   time slice, summarising `focus` only), `species` (one row per species per
#'   slice: `cells` over the species' study extent, `focus_cells` over `focus`,
#'   each differenced against the previous slice and against the baseline, plus
#'   `focus_suit_q90` and `focus_suit_margin` saying whether the focus clears
#'   the species' threshold at all), and `models` (fit diagnostics) -- plus the
#'   fitted models and projected ranges. See the two-geographies section before
#'   plotting `cells`, and the niche section before reading `focus_cells` as a
#'   range response.
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
                                min_cells = 100,
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
  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 0) {
    rc_abort("{.arg min_cells} must be a single non-negative number.")
  }
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
        ranges[[as.character(tt)]][[sp]] <- ps$ranges[[ps$representative]]
        cell_val <- stats::median(ps$cells)
        counts[[length(counts) + 1L]] <- tibble::tibble(
          species = sp, time = tt, cells = cell_val,
          present_cells = base_cells,
          cells_min = min(ps$cells), cells_max = max(ps$cells),
          cells_sd = stats::sd(ps$cells),
          focus_suit_q90 = focus_suit_q90(terra::unwrap(ps$suit), focus)
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
          present_cells = base_cells,
          focus_suit_q90 = focus_suit_q90(terra::unwrap(proj$suit), focus)
        )
      }
    }
  }
  if (!quiet) cli::cli_progress_done(.envir = environment())

  if (length(models) == 0) {
    rc_abort("No species could be modelled; nothing to summarise.")
  }

  # --- Resolvability ------------------------------------------------------
  # A range covering a handful of grid cells cannot support a trajectory. The
  # model sees only as many distinct climate vectors as the range has cells,
  # so whether any of them clears the threshold at a given slice turns on which
  # cells the draw happened to hit. The resulting series is an artefact of the
  # sample, not a measurement with error around it, and stacking it into
  # richness moves the assemblage total on that basis.
  #
  # Such species are fitted, projected and reported exactly as before -- they
  # stay in `$species`, `$models` and `$ranges`, so the underlying numbers
  # remain inspectable -- and are held out of the richness surfaces only.
  resolvability <- tibble::tibble(
    species     = names(models),
    range_cells = vapply(models, function(m) m$range_cells %||% NA_integer_,
                         integer(1)),
    min_cells   = as.integer(min_cells)
  )
  resolvability$resolvable <- resolvability$range_cells >= min_cells
  resolvability$reason <- ifelse(
    resolvability$resolvable, NA_character_,
    sprintf("range covers %d cells, below min_cells = %d",
            resolvability$range_cells, as.integer(min_cells))
  )
  keep <- resolvability$species[resolvability$resolvable]

  if (length(keep) == 0) {
    rc_abort(c(
      "No species clears {.arg min_cells} = {min_cells}; richness would be empty.",
      "i" = "Largest range covers {max(resolvability$range_cells)} cells.",
      "i" = "Lower {.arg min_cells}, or prepare climate at a finer resolution."
    ))
  }
  n_dropped <- nrow(resolvability) - length(keep)
  if (n_dropped > 0 && !quiet) {
    cli::cli_alert_warning(c(
      "{n_dropped} species held out of richness as not resolvable at this grid: {.val {setdiff(resolvability$species, keep)}}."
    ))
    cli::cli_alert_info("See {.code $resolvability}; they remain in {.code $species} and {.code $models}.")
  }

  # --- Richness per slice -------------------------------------------------
  # The surfaces are kept, not just their summaries: a mean over a region is a
  # poor substitute for seeing where in that region the species actually are.
  surfaces <- list()
  richness <- lapply(times, function(tt) {
    slice_ranges <- ranges[[as.character(tt)]]
    slice_ranges <- slice_ranges[intersect(names(slice_ranges), keep)]

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
        avail <- rep_ranges[[key]][intersect(names(rep_ranges[[key]]), keep)]
        per_sp <- lapply(avail, function(x) x[[r]])
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
  present_ranges <- baselines[intersect(names(models), keep)]
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

  # --- Per-species occupancy of the focus ---------------------------------
  # `cells` is measured over each species' own study extent, so for a
  # widespread taxon it is a continental number and moves for continental
  # reasons. Richness is a focus-only quantity. Without a focus-restricted
  # per-species count the two cannot be read against each other, and a
  # trajectory panel silently answers a different question from the richness
  # curve beside it. Every species is counted, held-out ones included, so the
  # column means the same thing in every row.
  focus_counts <- dplyr::bind_rows(lapply(times, function(tt) {
    fc <- focus_cells_each(ranges[[as.character(tt)]], focus,
                           resolution = resolution)
    tibble::tibble(species = names(fc), time = tt, focus_cells = unname(fc))
  }))
  focus_base <- focus_cells_each(baselines, focus, resolution = resolution)
  thresholds <- vapply(models, function(m) m$threshold, numeric(1))

  # --- Species-level deltas ----------------------------------------------
  species_tbl <- dplyr::bind_rows(counts) |>
    dplyr::left_join(focus_counts, by = c("species", "time")) |>
    dplyr::mutate(
      # Distance from the cutoff, in suitability units. Negative means the
      # focus does not clear this species' own threshold even at its ninetieth
      # percentile, so whatever cells it does contribute are an edge effect.
      focus_suit_margin = .data$focus_suit_q90 -
        unname(thresholds[.data$species]),
      # A species with no suitable area anywhere has no range polygon, so it
      # is absent from `ranges` and from `baselines` rather than present with
      # an empty geometry. That is a zero here, not an unknown -- but a
      # baseline that failed to project is genuinely unknown, and
      # `present_cells` already records which is which.
      focus_cells = dplyr::coalesce(.data$focus_cells, 0L),
      focus_present_cells = ifelse(
        is.na(.data$present_cells), NA_integer_,
        dplyr::coalesce(unname(focus_base[.data$species]), 0L)
      )
    ) |>
    dplyr::arrange(.data$species, .data$time) |>   # numeric time: correct lag
    dplyr::group_by(.data$species) |>
    dplyr::mutate(
      delta_from_previous = .data$cells - dplyr::lag(.data$cells),
      delta_from_present  = .data$cells - .data$present_cells,
      focus_delta_from_previous =
        .data$focus_cells - dplyr::lag(.data$focus_cells),
      focus_delta_from_present =
        .data$focus_cells - .data$focus_present_cells
    ) |>
    dplyr::ungroup()

  models_tbl <- tibble::tibble(
    species       = names(models),
    auc           = vapply(models, function(m) m$auc, numeric(1)),
    threshold     = vapply(models, function(m) m$threshold, numeric(1)),
    n_presence    = vapply(models, function(m) m$n_presence, integer(1)),
    n_background  = vapply(models, function(m) m$n_background, integer(1)),
    range_cells   = vapply(models, function(m) m$range_cells %||% NA_integer_,
                           integer(1)),
    present_cells = vapply(models, function(m) m$present_cells, integer(1)),
    resolvable    = resolvability$resolvable[match(names(models),
                                                   resolvability$species)]
  )

  structure(
    list(
      richness = richness,
      species  = species_tbl,
      models   = models_tbl,
      resolvability = resolvability,
      min_cells = as.integer(min_cells),
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

#' Report whether the resolvability screen ran, and with what setting
#'
#' Four states are distinguishable and only one of them is a clean pass. A
#' series written before the screen existed has no `min_cells` field at all;
#' one written with `min_cells = 0` ran with the screen disabled; one with a
#' positive setting either cleared everyone or held some species out. Row
#' counts cannot tell these apart -- "nothing listed" is the printed form of
#' both a clean pass and an absent test -- so the status is keyed on the
#' presence of the field and reported unconditionally.
#'
#' @param x A `richcast_series`.
#' @return A one-line character description, invisibly.
#' @noRd
screen_status <- function(x) {
  if (is.null(x$min_cells)) {
    return("resolvability: NOT RECORDED (series predates the screen)")
  }
  if (x$min_cells == 0) {
    return("resolvability: disabled (min_cells = 0)")
  }
  n_out <- sum(!x$resolvability$resolvable)
  sprintf("resolvability: ran at min_cells = %d, %d of %d held out",
          x$min_cells, n_out, nrow(x$resolvability))
}

#' @export
print.richcast_series <- function(x, ...) {
  cli::cli_text("{.cls richcast_series}")
  cli::cli_text("  {length(x$fits)} species x {length(x$times)} slices ({min(x$times)}-{max(x$times)} CE)")
  cli::cli_text("  focus: {x$focus$label}")
  cli::cli_text("  {screen_status(x)}")
  n_out <- sum(!x$resolvability$resolvable)
  if (n_out > 0) {
    cli::cli_text("  {n_out} species not resolvable at this grid, held out of richness ({.code $resolvability})")
  }
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
