# ==============================================================================
# Species distribution models
#
# The central design decision: fitting and projection are separate functions.
# A model depends only on present-day climate and the species' range, so it is
# fitted ONCE and then projected onto as many time slices as wanted. The source
# pipeline refitted inside its time loop, repeating identical work once per
# slice and overwriting each slice's diagnostics with the next.
#
# Rasters are stored wrapped (terra::wrap) and vectors as sf, so a fitted model
# survives saveRDS(). Bare SpatRaster/SpatVector objects hold external pointers
# and reload as invalid handles -- silently, until first use.
# ==============================================================================

#' Fit a species distribution model
#'
#' Fits a MaxEnt-style model (via \pkg{maxnet}) for one species against
#' present-day climate, and returns an object that [project_sdm()] can push
#' onto any number of palaeoclimate slices.
#'
#' Presence and background points are drawn from within and outside the
#' species' range polygon respectively, both restricted to a study extent
#' buffered around the range. When `db` carries point geometry (from
#' [gbif_occurrences()]), the occurrences are used directly and only the
#' background is sampled.
#'
#' @param db A `richcast_db` from [build_taxon_db()].
#' @param species Species name, or a row index into `db`.
#' @param climate A climate source from [climate_dir()] or [pastclim_climate()].
#' @param predictors Character vector of climate variables. `NULL` uses every
#'   variable the present-day slice provides.
#' @param nsample Number of presence and background points to draw. Measured
#'   across 32 rodent species under the default `p10` threshold, raising this
#'   from 100 to 1000 cut the worst-case coefficient of variation in modelled
#'   range size from 28.6% to 5.5%, at negligible cost.
#' @param buffer Study-extent buffer around the range, as a proportion of the
#'   range's bounding-box diagonal.
#' @param buffer_range Minimum and maximum buffer in degrees, clamping `buffer`.
#' @param extent Optional fixed `c(xmin, xmax, ymin, ymax)` study extent,
#'   overriding `buffer`. Fitting on the species' full range is usually
#'   preferable to fitting on the study region, since a truncated extent
#'   truncates the estimated niche.
#' @param threshold Binarisation rule. `"tss"` maximises the true skill
#'   statistic; `"p10"` uses the tenth percentile of predictions at training
#'   presences; `"mtp"` uses the minimum prediction at a training presence; or
#'   supply a number in `(0, 1)` for a fixed cutoff.
#'
#'   The choice matters more than it looks. `"tss"` is referenced to the
#'   background, so a wider background makes discrimination easier and pushes
#'   the cutoff up: measured across eight ungulates, widening the background
#'   buffer moved the median TSS threshold from 0.78 to 0.86 and the median
#'   predicted range from 18 cells to 2. `"p10"` and `"mtp"` are referenced to
#'   the presences instead, so they do not tighten as the background grows.
#' @param seed Random seed for point sampling.
#' @param land Optional `sf`/`sfc` land outline used to mask predictions and
#'   restrict the background to land. Defaults to Natural Earth at medium
#'   resolution.
#' @param quiet Suppress progress messages.
#' @return A `richcast_sdm`.
#' @seealso [project_sdm()], [run_hindcast_series()]
#' @export
fit_sdm <- function(db,
                    species,
                    climate,
                    predictors = NULL,
                    nsample = 1000,
                    buffer = 0.15,
                    buffer_range = c(1, 8),
                    extent = NULL,
                    threshold = "p10",
                    seed = 123,
                    land = NULL,
                    quiet = FALSE) {

  rlang::check_installed("maxnet", "to fit species distribution models.")
  say <- function(...) if (!quiet) cli::cli_alert_info(...)

  row <- db_row(db, species)
  sp_name <- row$species
  sp_geom <- sf::st_geometry(row)

  is_point <- as.character(sf::st_geometry_type(sp_geom, by_geometry = FALSE)) %in%
    c("POINT", "MULTIPOINT")

  # --- Study extent ------------------------------------------------------
  study_ext <- if (!is.null(extent)) {
    terra::ext(extent[1], extent[2], extent[3], extent[4])
  } else {
    buffered_extent(sp_geom, pct = buffer,
                    min_deg = buffer_range[1], max_deg = buffer_range[2],
                    quiet = quiet)
  }

  # --- Present-day climate ----------------------------------------------
  say("[{sp_name}] Loading present-day climate")
  present <- climate_at(climate, "present", predictors, study_ext)
  predictors <- names(present)
  if (length(predictors) == 0) {
    rc_abort("No climate variables available for the present-day slice.")
  }

  land_geom <- land %||% land_outline("medium")
  land_vec <- terra::vect(sf::st_sf(geometry = sf::st_geometry(land_geom)))
  sp_vec <- terra::vect(sf::st_sf(geometry = sp_geom))
  study_vec <- terra::vect(study_ext, crs = "EPSG:4326")

  # --- Presence and background samples -----------------------------------
  # Erasing the range from the study extent only makes sense for a polygon.
  # Points have no area to erase, and terra::erase() returns zero features
  # rather than the untouched extent, so routing occurrences through the same
  # call aborted every fit on "no background area left" -- a diagnosis exactly
  # backwards from the truth, since a point range fills nothing.
  #
  # With occurrences the background is the study extent itself. That is the
  # usual presence-background convention: the background describes what was
  # available, and a cell containing an occurrence was available too.
  bg_vec <- if (is_point) {
    terra::intersect(study_vec, land_vec)
  } else {
    erased <- terra::erase(study_vec, sp_vec)
    if (is.null(erased) || nrow(erased) == 0) {
      rc_abort(c(
        "[{sp_name}] No background area left after erasing the range.",
        "i" = "The range fills the study extent; widen {.arg buffer}."
      ))
    }
    terra::intersect(erased, land_vec)
  }
  if (is.null(bg_vec) || nrow(bg_vec) == 0) {
    rc_abort(c(
      "[{sp_name}] No background area left inside the study extent.",
      "i" = "The extent may fall entirely off the {.arg land} outline."
    ))
  }

  set.seed(seed)
  if (is_point) {
    # nrow(row) is 1 for a combined MULTIPOINT however many occurrences it
    # holds, so it cannot report the count; the geometry has to be counted.
    n_occ <- nrow(terra::geom(sp_vec))
    say("[{sp_name}] Using {n_occ} occurrence point{?s} directly")
    pres_df <- stats::na.omit(terra::extract(present, sp_vec, ID = FALSE))
  } else {
    pres_df <- sample_cells(present, sp_vec, nsample)
  }
  bg_df <- sample_cells(present, bg_vec, nsample)

  if (nrow(pres_df) < 5 || nrow(bg_df) < 5) {
    rc_abort(c(
      "[{sp_name}] Too few usable samples: {nrow(pres_df)} presence, {nrow(bg_df)} background.",
      "i" = "The range may fall entirely on cells with missing climate data."
    ))
  }
  say("[{sp_name}] Sampled {nrow(pres_df)} presence, {nrow(bg_df)} background")

  # A range narrower than nsample cells cannot yield nsample distinct samples.
  # terra warns about this per call, which buries the signal under noise across
  # a batch; report it once per species, naming the species, and carry the
  # realised counts in the returned object so they can be audited afterwards.
  if (nrow(pres_df) < nsample && !quiet) {
    cli::cli_alert_warning(
      "[{sp_name}] Range yielded {nrow(pres_df)} of {nsample} presence samples; the range is small relative to the grid."
    )
  }

  response <- c(rep(1, nrow(pres_df)), rep(0, nrow(bg_df)))
  covars <- rbind(pres_df[predictors], bg_df[predictors])

  # --- Fit ---------------------------------------------------------------
  model <- maxnet::maxnet(p = response, data = covars, regmult = 1)
  fitted <- as.numeric(stats::predict(model, newdata = covars, type = "cloglog"))

  auc <- if (requireNamespace("modEvA", quietly = TRUE)) {
    modEvA::AUC(obs = response, pred = fitted, plot = FALSE)$AUC
  } else {
    NA_real_
  }

  cutoff <- resolve_threshold(threshold, response, fitted, sp_name, quiet)
  say("[{sp_name}] AUC {round(auc, 3)}, threshold {round(cutoff, 3)}")

  # --- Present-day prediction --------------------------------------------
  suit <- terra::predict(present, model, type = "cloglog", na.rm = TRUE)
  suit <- terra::mask(suit, land_vec)
  binary <- binarise(suit, cutoff)

  structure(
    list(
      species        = sp_name,
      model          = model,
      predictors     = predictors,
      threshold      = cutoff,
      threshold_rule = if (is.character(threshold)) threshold else "fixed",
      auc            = auc,
      importance     = variable_importance(model),
      n_presence     = nrow(pres_df),
      n_background   = nrow(bg_df),
      study_extent   = as.vector(study_ext),
      present_cells  = binary$cells,
      present_range  = binary$polygon,
      present_suit   = terra::wrap(suit),
      land           = sf::st_geometry(land_geom),
      seed           = seed
    ),
    class = "richcast_sdm"
  )
}

#' Project a fitted model onto a palaeoclimate slice
#'
#' @param sdm A `richcast_sdm` from [fit_sdm()].
#' @param climate A climate source.
#' @param time Year CE to project onto.
#' @param window Optional [gaussian_window()] for time-averaging.
#' @param quiet Suppress progress messages.
#' @return A `richcast_projection`: the suitability surface, the binarised
#'   range as `sf`, and the suitable-cell count.
#' @seealso [fit_sdm()]
#' @export
project_sdm <- function(sdm, climate, time, window = NULL, quiet = FALSE) {

  if (!inherits(sdm, "richcast_sdm")) {
    rc_abort("{.arg sdm} must come from {.fn fit_sdm}.")
  }
  say <- function(...) if (!quiet) cli::cli_alert_info(...)

  study_ext <- terra::ext(sdm$study_extent)
  say("[{sdm$species}] Projecting onto {time} CE")

  past <- climate_for_projection(
    climate, time, sdm$predictors, study_ext, window = window, quiet = quiet
  )

  missing <- setdiff(sdm$predictors, names(past))
  if (length(missing) > 0) {
    rc_abort(c(
      "[{sdm$species}] Predictor{?s} {.val {missing}} absent from the {time} CE slice.",
      "i" = "A model can only be projected onto the variables it was fitted with."
    ))
  }
  past <- past[[sdm$predictors]]

  scored <- score_raster(past, sdm$model, sdm$threshold, sdm$land)
  suit <- scored$suit
  binary <- scored$binary

  if (binary$cells == 0 && !quiet) {
    cli::cli_alert_warning(
      "[{sdm$species}] No cells above threshold at {time} CE; range is empty."
    )
  }

  structure(
    list(
      species   = sdm$species,
      time      = time,
      threshold = sdm$threshold,
      cells     = binary$cells,
      range     = binary$polygon,
      suit      = terra::wrap(suit),
      window    = window
    ),
    class = "richcast_projection"
  )
}

#' Extract the suitability surface from a fitted model or projection
#'
#' Rasters are stored wrapped so that objects survive `saveRDS()`; this
#' unwraps one back into a usable `SpatRaster`.
#'
#' @param x A `richcast_sdm` or `richcast_projection`.
#' @return A `SpatRaster`.
#' @export
suitability <- function(x) {
  r <- switch(
    class(x)[1],
    richcast_sdm        = x$present_suit,
    richcast_projection = x$suit,
    rc_abort("{.arg x} must be a {.cls richcast_sdm} or {.cls richcast_projection}.")
  )
  terra::unwrap(r)
}

#' @export
print.richcast_sdm <- function(x, ...) {
  cli::cli_text("{.cls richcast_sdm} {.strong {x$species}}")
  cli::cli_text("  {length(x$predictors)} predictors: {.val {x$predictors}}")
  cli::cli_text("  AUC {round(x$auc, 3)} | threshold {round(x$threshold, 3)} ({x$threshold_rule})")
  cli::cli_text("  present-day suitable cells: {x$present_cells}")
  invisible(x)
}

#' @export
print.richcast_projection <- function(x, ...) {
  cli::cli_text("{.cls richcast_projection} {.strong {x$species}} @ {x$time} CE")
  cli::cli_text("  suitable cells: {x$cells}")
  invisible(x)
}


# ==============================================================================
# Internals
# ==============================================================================

#' Pull one species row from a taxon database
#' @noRd
db_row <- function(db, species) {
  if (is.numeric(species)) {
    if (species < 1 || species > nrow(db)) {
      rc_abort("Row index {species} is out of range (db has {nrow(db)} rows).")
    }
    return(db[species, ])
  }
  key <- normalise_species(species)
  hit <- which(db$species == key)
  if (length(hit) == 0) {
    rc_abort(c(
      "Species {.val {species}} not found in the database.",
      "i" = "Names are normalised by {.fn normalise_species}; looked for {.val {key}}."
    ))
  }
  if (length(hit) > 1) {
    rows <- db[hit, ]
    types <- as.character(sf::st_geometry_type(rows))

    # A point source keeps one row per occurrence on purpose:
    # build_taxon_db() skips the dissolve for points, so several rows for one
    # species is the normal shape of occurrence data, not a defect to warn
    # about. Taking the first would have fitted the model to a single
    # occurrence, and the advice to rebuild with dissolve = TRUE could not have
    # helped, since that path never runs for points.
    if (all(types %in% c("POINT", "MULTIPOINT"))) {
      out <- rows[1, ]
      sf::st_geometry(out) <- sf::st_combine(sf::st_geometry(rows))
      return(out)
    }

    cli::cli_warn(c(
      "{.val {key}} has {length(hit)} rows; using the first.",
      "i" = "Rebuild with {.code build_taxon_db(dissolve = TRUE)} to union them."
    ))
  }
  db[hit[1], ]
}

#' Predict one model onto an already-loaded raster and binarise
#'
#' Factored out so a set of replicate models can be scored against a single
#' raster read. Raster I/O dominates projection cost, so loading a slice once
#' and predicting N models against it is close to N times cheaper than calling
#' [project_sdm()] N times.
#'
#' @param past A `SpatRaster` of predictors, already cropped and ordered.
#' @param model A fitted maxnet model.
#' @param threshold Binarisation cutoff.
#' @param land_geom `sfc` land outline used to mask.
#' @return A list with `suit` (SpatRaster) and `binary` (from [binarise()]).
#' @noRd
score_raster <- function(past, model, threshold, land_geom) {
  suit <- terra::predict(past, model, type = "cloglog", na.rm = TRUE)
  suit <- terra::mask(suit, terra::vect(sf::st_sf(geometry = land_geom)))
  list(suit = suit, binary = binarise(suit, threshold))
}

#' Sample climate values from the cells covered by a polygon
#'
#' terra emits "[spatSample] fewer cells returned than requested" whenever the
#' mask holds fewer cells than asked for. That is expected for narrow-ranged
#' species and is reported by the caller in species-labelled form, so the raw
#' warning is suppressed here rather than repeated once per call.
#' @noRd
sample_cells <- function(r, poly, n) {
  as.data.frame(
    withCallingHandlers(
      terra::spatSample(
        terra::mask(r, poly), n,
        replace = TRUE, na.rm = TRUE, as.points = TRUE
      ),
      warning = function(w) {
        if (grepl("fewer cells returned", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
}

#' Study extent as a proportion of the range's bounding-box diagonal
#'
#' Narrow endemics get a tight background, wide-ranging species a broad one.
#' @noRd
buffered_extent <- function(geom, pct = 0.15, min_deg = 1, max_deg = 8,
                            quiet = FALSE) {
  bb <- sf::st_bbox(geom)
  diag <- sqrt((bb[["xmax"]] - bb[["xmin"]])^2 + (bb[["ymax"]] - bb[["ymin"]])^2)
  buf <- min(max(pct * diag, min_deg), max_deg)
  if (!quiet) {
    cli::cli_alert_info(
      "Range diagonal {round(diag, 1)} deg; buffer {round(buf, 2)} deg."
    )
  }
  terra::ext(
    bb[["xmin"]] - buf, bb[["xmax"]] + buf,
    bb[["ymin"]] - buf, bb[["ymax"]] + buf
  )
}

#' Resolve the binarisation threshold
#' @noRd
resolve_threshold <- function(threshold, obs, pred, sp_name, quiet) {
  if (is.numeric(threshold)) {
    if (threshold <= 0 || threshold >= 1) {
      rc_abort("A fixed {.arg threshold} must lie strictly between 0 and 1.")
    }
    return(threshold)
  }

  presence_pred <- pred[obs == 1]

  # Presence-referenced rules. Unlike TSS these ask "what score do known
  # occurrences achieve?" rather than "what score best separates presence from
  # background?", so they do not tighten as the background widens.
  if (identical(threshold, "p10")) {
    # Tenth-percentile training presence: tolerate 10% omission, the usual
    # choice where occurrence data carry locational error.
    return(unname(stats::quantile(presence_pred, 0.10, na.rm = TRUE)))
  }
  if (identical(threshold, "mtp")) {
    # Minimum training presence: the most permissive rule, admitting anywhere
    # at least as suitable as the worst known occurrence.
    return(min(presence_pred, na.rm = TRUE))
  }

  if (!identical(threshold, "tss")) {
    rc_abort('{.arg threshold} must be "tss", "p10", "mtp", or a number in (0, 1).')
  }
  if (!requireNamespace("modEvA", quietly = TRUE)) {
    cli::cli_warn(
      "{.pkg modEvA} not installed; falling back to a fixed threshold of 0.5."
    )
    return(0.5)
  }
  opt <- modEvA::optiThresh(obs = obs, pred = pred, measures = "TSS", plot = FALSE)
  opt$optimals.each$value[1]
}

#' Binarise a suitability surface and vectorise the suitable area
#'
#' Uses `ifel()` to produce 1/NA rather than 1/0, so `as.polygons()` returns
#' only the suitable area instead of a two-class coverage of the whole extent.
#' The polygon is returned as sf, which serialises.
#' @noRd
binarise <- function(suit, cutoff) {
  above <- suit > cutoff
  cells <- terra::global(above, "sum", na.rm = TRUE)[1, 1]
  cells <- if (is.na(cells)) 0 else as.integer(cells)

  polygon <- NULL
  if (cells > 0) {
    p <- terra::as.polygons(terra::ifel(above, 1, NA), dissolve = TRUE, na.rm = TRUE)
    if (!is.null(p) && nrow(p) > 0) {
      polygon <- sf::st_geometry(sf::st_as_sf(p))
    }
  }
  list(cells = cells, polygon = polygon)
}

#' Summed absolute lambda per variable, as a crude importance measure
#' @noRd
variable_importance <- function(model) {
  lam <- model$betas
  if (is.null(lam) || length(lam) == 0) return(numeric(0))

  base_var <- function(nm) {
    if (grepl("^hinge\\(", nm)) sub("^hinge\\(([^,]+),.*", "\\1", nm)
    else if (grepl("^I\\(", nm)) sub("^I\\(([^\\^]+).*", "\\1", nm)
    else nm
  }
  contrib <- tapply(abs(lam), vapply(names(lam), base_var, character(1)), sum)
  sort(contrib, decreasing = TRUE)
}
