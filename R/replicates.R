# ==============================================================================
# Replicate ensembles
#
# A single fit is one draw. Redrawing the presence and background samples moves
# the answer a great deal -- for one Tian Shan marmot, six seeds at the default
# nsample = 100 spanned a 2.4-fold range in present-day extent and disagreed on
# the SIGN of the change at one slice. A point estimate with no interval around
# it hides that.
#
# What this quantifies is precision, not accuracy: how much the answer moves
# when the pseudo-occurrences are redrawn from the same range polygon under the
# same climate reconstruction. It says nothing about whether the polygon, the
# reconstruction, the niche-conservatism assumption, or maxnet is right. Those
# are larger, and `ensemble_series()` is the tool for the climate part.
# ==============================================================================

#' Fit replicate models for one species
#'
#' Refits the same species repeatedly, redrawing presence and background points
#' each time and re-optimising the threshold, so the spread across replicates
#' can be reported.
#'
#' @param db,species,climate As for [fit_sdm()].
#' @param replicates Number of replicate fits.
#' @param seed Base seed. Replicate `i` uses `seed + i - 1`, so a run is
#'   reproducible and replicate seeds never collide.
#' @param quiet Suppress per-replicate messages.
#' @param ... Further arguments passed to [fit_sdm()].
#' @return A `richcast_sdm_set`: a list of `richcast_sdm` objects plus the
#'   species name.
#' @seealso [project_replicates()], [run_hindcast_series()]
#' @export
fit_replicates <- function(db, species, climate, replicates = 20,
                           seed = 123, quiet = FALSE, ...) {

  if (!is.numeric(replicates) || length(replicates) != 1 || replicates < 1) {
    rc_abort("{.arg replicates} must be a single positive number.")
  }
  replicates <- as.integer(replicates)

  fits <- vector("list", replicates)
  for (i in seq_len(replicates)) {
    fits[[i]] <- fit_sdm(db, species, climate, seed = seed + i - 1L,
                         quiet = TRUE, ...)
  }

  nm <- fits[[1]]$species
  if (!quiet) {
    pc <- vapply(fits, function(f) f$present_cells, integer(1))
    cli::cli_alert_info(
      "[{nm}] {replicates} replicates: present-day cells {min(pc)}-{max(pc)} (CV {round(100 * stats::sd(pc) / mean(pc), 1)}%)."
    )
  }

  structure(list(species = nm, fits = fits, replicates = replicates,
                 seed = seed),
            class = "richcast_sdm_set")
}

#' @export
print.richcast_sdm_set <- function(x, ...) {
  pc <- vapply(x$fits, function(f) f$present_cells, integer(1))
  cli::cli_text("{.cls richcast_sdm_set} {.strong {x$species}}: {x$replicates} replicates")
  cli::cli_text("  present-day cells: {min(pc)}-{max(pc)} (median {stats::median(pc)})")
  invisible(x)
}

#' Project a replicate set onto one time slice
#'
#' Reads the slice once and scores every replicate against it. That ordering is
#' the point: raster I/O dominates, so this is far cheaper than calling
#' [project_sdm()] once per replicate.
#'
#' @param set A `richcast_sdm_set` from [fit_replicates()].
#' @param climate A climate source.
#' @param time Year CE, or `"present"`.
#' @param window Optional [gaussian_window()].
#' @param quiet Suppress messages.
#' @return A `richcast_projection_set`: `cells` (one per replicate) and
#'   `ranges` (one `sfc` or `NULL` per replicate).
#' @seealso [fit_replicates()]
#' @export
project_replicates <- function(set, climate, time, window = NULL,
                               quiet = FALSE) {

  if (!inherits(set, "richcast_sdm_set")) {
    rc_abort("{.arg set} must come from {.fn fit_replicates}.")
  }
  ref <- set$fits[[1]]

  # All replicates share the study extent and predictors, so one read serves
  # the whole set.
  past <- climate_for_projection(
    climate, time, ref$predictors, terra::ext(ref$study_extent),
    window = window, quiet = quiet
  )
  missing <- setdiff(ref$predictors, names(past))
  if (length(missing) > 0) {
    rc_abort(c(
      "[{set$species}] Predictor{?s} {.val {missing}} absent from the {time} slice.",
      "i" = "A model can only be projected onto the variables it was fitted with."
    ))
  }
  past <- past[[ref$predictors]]

  cells  <- integer(set$replicates)
  ranges <- vector("list", set$replicates)
  for (i in seq_len(set$replicates)) {
    f <- set$fits[[i]]
    b <- score_raster(past, f$model, f$threshold, f$land)$binary
    cells[i]    <- b$cells
    ranges[[i]] <- b$polygon
  }

  structure(list(species = set$species, time = time, cells = cells,
                 ranges = ranges, replicates = set$replicates,
                 window = window),
            class = "richcast_projection_set")
}

#' @export
print.richcast_projection_set <- function(x, ...) {
  cli::cli_text("{.cls richcast_projection_set} {.strong {x$species}} @ {x$time}")
  cli::cli_text("  suitable cells across {x$replicates} replicates: {min(x$cells)}-{max(x$cells)}")
  invisible(x)
}


# ==============================================================================
# Combining independent series
# ==============================================================================

#' Combine several series into an ensemble
#'
#' Summarises richness across runs that differ in something structural -- most
#' usefully the climate reconstruction. Where replicate intervals describe how
#' much redrawing the sample moves the answer, this describes how much the
#' choice of input does, and in this pipeline the latter has been the larger of
#' the two.
#'
#' The spread reported is the range and standard deviation across members at
#' each slice, not a confidence interval: a handful of climate products is not
#' a random sample from a population of products, so quantiles of them would
#' invite a reading they cannot support.
#'
#' @param series A named list of `richcast_series` objects, all run over the
#'   same time slices and focus.
#' @return A tibble with one row per time slice: the per-member mean richness,
#'   plus `ens_mean`, `ens_min`, `ens_max` and `ens_sd` across members.
#' @seealso [run_hindcast_series()]
#' @examples
#' \dontrun{
#' ensemble_series(list(chelsa = res_chelsa, worldclim = res_worldclim))
#' }
#' @export
ensemble_series <- function(series) {

  if (!is.list(series) || length(series) < 2) {
    rc_abort("{.arg series} must be a list of at least two {.cls richcast_series}.")
  }
  if (!all(vapply(series, inherits, logical(1), "richcast_series"))) {
    rc_abort("Every element of {.arg series} must be a {.cls richcast_series}.")
  }
  if (is.null(names(series)) || any(!nzchar(names(series)))) {
    rc_abort("{.arg series} must be named, so members can be told apart.")
  }

  times <- lapply(series, function(s) sort(s$times))
  if (length(unique(lapply(times, as.numeric))) > 1) {
    cli::cli_warn(c(
      "Members do not cover the same time slices.",
      "i" = "Only slices present in every member are summarised."
    ))
  }
  common <- Reduce(intersect, times)
  if (length(common) == 0) {
    rc_abort("Members share no time slices.")
  }

  long <- dplyr::bind_rows(lapply(names(series), function(nm) {
    series[[nm]]$richness |>
      dplyr::filter(!is.na(.data$time), .data$time %in% common) |>
      dplyr::transmute(member = nm, time = .data$time,
                       mean_richness = .data$mean_richness)
  }))

  wide <- long |>
    tidyr_pivot_wider(names_from = "member", values_from = "mean_richness")

  member_cols <- setdiff(names(wide), "time")
  m <- as.matrix(wide[member_cols])

  wide$ens_mean <- rowMeans(m, na.rm = TRUE)
  wide$ens_min  <- apply(m, 1, min, na.rm = TRUE)
  wide$ens_max  <- apply(m, 1, max, na.rm = TRUE)
  wide$ens_sd   <- apply(m, 1, stats::sd, na.rm = TRUE)
  wide$ens_range <- wide$ens_max - wide$ens_min

  dplyr::arrange(wide, .data$time)
}

#' Reshape long to wide without taking a tidyr dependency
#' @noRd
tidyr_pivot_wider <- function(data, names_from, values_from) {
  keys <- setdiff(names(data), c(names_from, values_from))
  out <- unique(data[keys])
  for (nm in unique(data[[names_from]])) {
    sub <- data[data[[names_from]] == nm, c(keys, values_from), drop = FALSE]
    names(sub)[names(sub) == values_from] <- nm
    out <- dplyr::left_join(out, sub, by = keys)
  }
  tibble::as_tibble(out)
}
