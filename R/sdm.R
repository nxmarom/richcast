# ==============================================================================
# Species distribution models
#
# One model per species: an equally weighted ensemble of a random forest
# (ranger) and MaxEnt (maxnet), trained on pseudo-presences drawn from inside
# the species' range polygon against background drawn from the rest of its
# study extent.
#
# Fitting and projection are separate functions. A model depends only on
# present-day climate and the species' range, so it is fitted ONCE and then
# projected onto as many time slices as wanted.
#
# Rasters are stored wrapped (terra::wrap) and vectors as sf, so a fitted model
# survives saveRDS(). Bare SpatRaster/SpatVector objects hold external pointers
# and reload as invalid handles -- silently, until first use.
# ==============================================================================

#' Fit a species distribution model
#'
#' Fits an ensemble of a random forest (via \pkg{ranger}) and MaxEnt (via
#' \pkg{maxnet}) for one species against present-day climate, and returns an
#' object that [project_sdm()] can push onto any number of palaeoclimate
#' slices.
#'
#' # Training data
#'
#' The study extent is the range polygon's bounding box, widened on every side
#' by `buffer` times its diagonal (30% by default). Inside it, `n_presence`
#' pseudo-presences are drawn from grid cells inside the range polygon and
#' `n_background` background points from land cells outside it, both without
#' replacement. A range covering fewer cells than `n_presence` contributes
#' every cell it has; `range_cells` on the result says how many that was.
#'
#' # Members and ensemble
#'
#' * **MaxEnt**: `maxnet::maxnet()` with default features and `regmult = 1`,
#'   predicted on the cloglog scale.
#' * **Random forest**: a probability forest from `ranger::ranger()`, with each
#'   tree drawn from an equal number of presences and background points
#'   (balanced down-sampling), which keeps the 1:10 class imbalance from
#'   flattening predicted probabilities toward zero.
#' * **Ensemble**: the unweighted mean of the two, so it is also on `[0, 1]`.
#'
#' The presence/absence threshold is applied to the ensemble only. The default
#' `"p10"` is the tenth percentile of ensemble predictions at the training
#' presences: it tolerates 10% omission, and being referenced to the presences
#' it does not tighten as the background widens.
#'
#' # Evaluation
#'
#' Before the final fit, a stratified `test_frac` of the presences and
#' background is held out, both members are trained on the rest, and each
#' member and the ensemble are scored on the held-out points:
#'
#' * **AUC**, the probability that a held-out presence scores above a held-out
#'   background point.
#' * **Continuous Boyce index** (Hirzel et al. 2006): the Spearman correlation
#'   between suitability and the predicted-to-expected ratio of held-out
#'   presences, over a moving window. Near 1 means presences concentrate where
#'   the model says suitability is high; near 0 means no better than random.
#'
#' The final members are then refitted on all the points. Set `test_frac = 0`
#' to skip evaluation.
#'
#' @param db A `richcast_db` from [build_taxon_db()].
#' @param species Species name, or a row index into `db`.
#' @param climate A climate source from [climate_dir()] or [pastclim_climate()].
#' @param predictors Character vector of climate variables. `NULL` uses every
#'   variable the present-day slice provides.
#' @param n_presence Number of pseudo-presences drawn inside the range.
#' @param n_background Number of background points drawn outside it.
#' @param buffer Study-extent buffer around the range's bounding box, as a
#'   proportion of the box's diagonal.
#' @param threshold `"p10"`, or a fixed number in `(0, 1)`.
#' @param test_frac Share of presences and background held out for evaluation.
#' @param num_trees Trees in the random forest.
#' @param seed Random seed for sampling, the hold-out split and the forest.
#' @param land Optional `sf`/`sfc` land outline used to mask predictions and
#'   restrict the background to land. Defaults to Natural Earth at medium
#'   resolution.
#' @param quiet Suppress progress messages.
#' @return A `richcast_sdm`. `metrics` holds AUC and Boyce for `maxent`, `rf`
#'   and `ensemble`; `range_cells` is the number of grid cells with climate the
#'   range covers.
#' @references Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A.
#'   (2006). Evaluating the ability of habitat suitability models to predict
#'   species presences. *Ecological Modelling*, 199, 142-152.
#' @seealso [project_sdm()], [run_hindcast_series()]
#' @export
fit_sdm <- function(db,
                    species,
                    climate,
                    predictors = NULL,
                    n_presence = 100,
                    n_background = 1000,
                    buffer = 0.3,
                    threshold = "p10",
                    test_frac = 0.25,
                    num_trees = 500,
                    seed = 123,
                    land = NULL,
                    quiet = FALSE) {

  rlang::check_installed(c("maxnet", "ranger"),
                         "to fit species distribution models.")
  say <- function(...) if (!quiet) cli::cli_alert_info(...)
  if (!is.numeric(test_frac) || length(test_frac) != 1 ||
      test_frac < 0 || test_frac >= 1) {
    rc_abort("{.arg test_frac} must be a single number in [0, 1).")
  }

  row <- db_row(db, species)
  sp_name <- row$species
  sp_geom <- sf::st_geometry(row)

  # --- Study extent and present-day climate -------------------------------
  study_ext <- study_extent(sp_geom, buffer)
  say("[{sp_name}] Loading present-day climate")
  present <- climate_at(climate, "present", predictors, study_ext)
  predictors <- names(present)
  if (length(predictors) == 0) {
    rc_abort("No climate variables available for the present-day slice.")
  }

  land_geom <- land %||% land_outline("medium")
  land_vec <- terra::vect(sf::st_sf(geometry = sf::st_geometry(land_geom)))
  # Range maps that are valid for sf can still break GEOS inside terra's
  # erase() ("unable to assign free hole to a shell"), so repair them in
  # terra's own terms first.
  sp_vec <- terra::makeValid(terra::vect(sf::st_sf(geometry = sp_geom)))
  study_vec <- terra::vect(study_ext, crs = "EPSG:4326")

  bg_vec <- terra::erase(study_vec, sp_vec)
  if (!is.null(bg_vec) && nrow(bg_vec) > 0) bg_vec <- terra::intersect(bg_vec, land_vec)
  if (is.null(bg_vec) || nrow(bg_vec) == 0) {
    rc_abort(c(
      "[{sp_name}] No background area left in the study extent.",
      "i" = "The range fills the extent, or the extent is off the {.arg land} outline; widen {.arg buffer}."
    ))
  }

  # --- Pseudo-presences and background ------------------------------------
  set.seed(seed)
  range_cells <- count_cells(present, sp_vec)
  pres_df <- sample_cells(present, sp_vec, n_presence)
  bg_df <- sample_cells(present, bg_vec, n_background)
  if (nrow(pres_df) < 5 || nrow(bg_df) < 5) {
    rc_abort(c(
      "[{sp_name}] Too few usable samples: {nrow(pres_df)} presence, {nrow(bg_df)} background.",
      "i" = "The range may fall on cells with missing climate data, or be smaller than the grid."
    ))
  }
  if (!quiet && nrow(pres_df) < n_presence) {
    cli::cli_alert_warning(
      "[{sp_name}] Range covers only {range_cells} grid cell{?s}; using {nrow(pres_df)} pseudo-presence{?s}."
    )
  }
  say("[{sp_name}] {nrow(pres_df)} pseudo-presences, {nrow(bg_df)} background")

  response <- c(rep(1L, nrow(pres_df)), rep(0L, nrow(bg_df)))
  covars <- rbind(pres_df[predictors], bg_df[predictors])

  # --- Hold-out evaluation ------------------------------------------------
  metrics <- if (test_frac > 0) {
    evaluate_holdout(response, covars, test_frac, num_trees, seed)
  } else {
    tibble::tibble(model = c("maxent", "rf", "ensemble"),
                   auc = NA_real_, boyce = NA_real_)
  }

  # --- Final fit on every point -------------------------------------------
  members <- fit_members(response, covars, num_trees, seed)
  fitted <- predict_members(members, covars)
  cutoff <- resolve_threshold(threshold, response, fitted[, "ensemble"])

  ens <- metrics[metrics$model == "ensemble", ]
  say("[{sp_name}] ensemble AUC {round(ens$auc, 3)}, Boyce {round(ens$boyce, 3)}, threshold {round(cutoff, 3)}")

  # --- Present-day prediction ---------------------------------------------
  suit <- predict_surface(present, members, land_geom)

  structure(
    list(
      species        = sp_name,
      members        = members,
      predictors     = predictors,
      threshold      = cutoff,
      threshold_rule = if (is.character(threshold)) threshold else "fixed",
      metrics        = metrics,
      n_presence     = nrow(pres_df),
      n_background   = nrow(bg_df),
      range_cells    = range_cells,
      study_extent   = as.vector(study_ext),
      present_cells  = count_above(suit, cutoff),
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
#' @param time Year CE to project onto, or `"present"`.
#' @param quiet Suppress progress messages.
#' @return A `richcast_projection`: the ensemble suitability surface over the
#'   species' study extent, the threshold, and the count of cells above it.
#' @seealso [fit_sdm()]
#' @export
project_sdm <- function(sdm, climate, time, quiet = FALSE) {

  if (!inherits(sdm, "richcast_sdm")) {
    rc_abort("{.arg sdm} must come from {.fn fit_sdm}.")
  }
  if (!quiet) cli::cli_alert_info("[{sdm$species}] Projecting onto {time} CE")

  past <- climate_at(climate, time, sdm$predictors, terra::ext(sdm$study_extent))
  missing <- setdiff(sdm$predictors, names(past))
  if (length(missing) > 0) {
    rc_abort(c(
      "[{sdm$species}] Predictor{?s} {.val {missing}} absent from the {time} CE slice.",
      "i" = "A model can only be projected onto the variables it was fitted with."
    ))
  }
  suit <- predict_surface(past[[sdm$predictors]], sdm$members, sdm$land)
  cells <- count_above(suit, sdm$threshold)

  if (cells == 0 && !quiet) {
    cli::cli_alert_warning(
      "[{sdm$species}] No cells above threshold at {time} CE; range is empty."
    )
  }

  structure(
    list(
      species   = sdm$species,
      time      = time,
      threshold = sdm$threshold,
      cells     = cells,
      suit      = terra::wrap(suit)
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
#' @param binary Return presence (1) / absence (0) at the model's threshold
#'   instead of continuous suitability.
#' @return A `SpatRaster`.
#' @export
suitability <- function(x, binary = FALSE) {
  r <- switch(
    class(x)[1],
    richcast_sdm        = x$present_suit,
    richcast_projection = x$suit,
    rc_abort("{.arg x} must be a {.cls richcast_sdm} or {.cls richcast_projection}.")
  )
  r <- terra::unwrap(r)
  if (binary) r <- terra::ifel(r > x$threshold, 1, 0)
  r
}

#' @export
print.richcast_sdm <- function(x, ...) {
  m <- x$metrics
  cli::cli_text("{.cls richcast_sdm} {.strong {x$species}}")
  cli::cli_text("  {length(x$predictors)} predictors: {.val {x$predictors}}")
  cli::cli_text("  {x$n_presence} pseudo-presences ({x$range_cells} range cell{?s}), {x$n_background} background")
  for (i in seq_len(nrow(m))) {
    cli::cli_text("  {m$model[i]}: AUC {round(m$auc[i], 3)} | Boyce {round(m$boyce[i], 3)}")
  }
  cli::cli_text("  threshold {round(x$threshold, 3)} ({x$threshold_rule}); present-day cells: {x$present_cells}")
  invisible(x)
}

#' @export
print.richcast_projection <- function(x, ...) {
  cli::cli_text("{.cls richcast_projection} {.strong {x$species}} @ {x$time}")
  cli::cli_text("  suitable cells: {x$cells}")
  invisible(x)
}


# ==============================================================================
# Model internals
# ==============================================================================

#' Fit both ensemble members
#' @noRd
fit_members <- function(response, covars, num_trees, seed) {
  covars <- as.data.frame(covars)
  maxent <- maxnet::maxnet(p = response, data = covars, regmult = 1)

  # Balanced down-sampling: every tree sees as many background points as
  # presences. Class order follows the factor levels, "0" then "1".
  n_min <- min(sum(response == 1), sum(response == 0))
  frac <- rep(n_min / length(response), 2)
  rf <- ranger::ranger(
    x = covars, y = factor(response, levels = c(0, 1)),
    probability = TRUE, num.trees = num_trees,
    sample.fraction = frac, replace = TRUE,
    seed = seed, num.threads = 1
  )
  rf$predictions <- NULL  # out-of-bag predictions, not needed afterwards
  list(maxent = maxent, rf = rf)
}

#' Predict members and ensemble for a data frame of predictors
#' @return A matrix with columns `maxent`, `rf` and `ensemble`.
#' @noRd
predict_members <- function(members, newdata) {
  newdata <- as.data.frame(newdata)
  n <- nrow(newdata)
  if (n == 0) {
    return(matrix(numeric(0), ncol = 3,
                  dimnames = list(NULL, c("maxent", "rf", "ensemble"))))
  }
  mx <- as.numeric(stats::predict(members$maxent, newdata = newdata,
                                  type = "cloglog"))
  rf <- stats::predict(members$rf, data = newdata,
                       num.threads = 1)$predictions[, "1"]
  cbind(maxent = mx, rf = rf, ensemble = (mx + rf) / 2)
}

#' Ensemble suitability over a predictor raster, masked to land
#' @noRd
predict_surface <- function(r, members, land_geom) {
  suit <- terra::predict(
    r, members,
    fun = function(model, data, ...) predict_members(model, data)[, "ensemble"],
    na.rm = TRUE
  )
  names(suit) <- "suitability"
  terra::mask(suit, terra::vect(sf::st_sf(geometry = land_geom)))
}

#' Score members on a stratified hold-out
#' @noRd
evaluate_holdout <- function(response, covars, test_frac, num_trees, seed) {
  set.seed(seed)
  pick <- function(idx) idx[sample.int(length(idx), max(1, round(length(idx) * test_frac)))]
  test <- c(pick(which(response == 1)), pick(which(response == 0)))

  members <- fit_members(response[-test], covars[-test, , drop = FALSE],
                         num_trees, seed)
  pred <- predict_members(members, covars[test, , drop = FALSE])
  obs <- response[test]

  models <- colnames(pred)
  tibble::tibble(
    model = models,
    auc   = vapply(models, function(m) auc_score(pred[obs == 1, m], pred[obs == 0, m]),
                   numeric(1), USE.NAMES = FALSE),
    boyce = vapply(models, function(m) boyce_index(pred[obs == 1, m], pred[, m]),
                   numeric(1), USE.NAMES = FALSE)
  )
}

#' Area under the ROC curve, by the Mann-Whitney statistic
#'
#' Ties count half, so a model that scores everything equally gets 0.5.
#' @param pres,bg Predictions at presences and at background points.
#' @noRd
auc_score <- function(pres, bg) {
  pres <- pres[!is.na(pres)]
  bg <- bg[!is.na(bg)]
  n1 <- length(pres)
  n0 <- length(bg)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pres, bg))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Continuous Boyce index
#'
#' Hirzel et al. (2006), as implemented in `ecospat::ecospat.boyce()`: a
#' window a tenth of the suitability range wide slides across it in `res`
#' steps; in each, the share of presences is divided by the share of all
#' evaluated points (the expected share under a random model), and the index is
#' the Spearman correlation of that ratio with window position. Windows with no
#' evaluated points are dropped, as are runs of an unchanged ratio, which
#' would otherwise inflate the correlation with ties.
#'
#' @param pres Predictions at presences.
#' @param fit Predictions at every evaluated point (presences and background),
#'   standing in for the suitability available across the study extent.
#' @param res Number of window positions.
#' @return A number in `[-1, 1]`, or `NA` if too few windows are informative.
#' @noRd
boyce_index <- function(pres, fit, res = 100) {
  pres <- pres[!is.na(pres)]
  fit <- fit[!is.na(fit)]
  if (length(pres) == 0 || length(fit) == 0) return(NA_real_)
  lo <- min(fit)
  hi <- max(fit)
  if (hi <= lo) return(NA_real_)
  width <- (hi - lo) / 10
  starts <- seq(lo, hi - width, length.out = res)

  ratio <- vapply(starts, function(s) {
    expected <- mean(fit >= s & fit <= s + width)
    if (expected == 0) return(NA_real_)
    mean(pres >= s & pres <= s + width) / expected
  }, numeric(1))
  mid <- starts + width / 2

  ok <- !is.na(ratio)
  ratio <- ratio[ok]
  mid <- mid[ok]
  keep <- c(TRUE, diff(ratio) != 0)
  ratio <- ratio[keep]
  mid <- mid[keep]
  if (length(ratio) < 3 || stats::sd(ratio) == 0) return(NA_real_)
  stats::cor(ratio, mid, method = "spearman")
}

#' Resolve the presence/absence threshold
#' @noRd
resolve_threshold <- function(threshold, obs, pred) {
  if (is.numeric(threshold)) {
    if (length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      rc_abort("A fixed {.arg threshold} must lie strictly between 0 and 1.")
    }
    return(threshold)
  }
  if (!identical(threshold, "p10")) {
    rc_abort('{.arg threshold} must be "p10" or a number in (0, 1).')
  }
  # Tenth-percentile training presence: tolerate 10% omission.
  unname(stats::quantile(pred[obs == 1], 0.10, na.rm = TRUE))
}


# ==============================================================================
# Spatial internals
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
    cli::cli_warn(c(
      "{.val {key}} has {length(hit)} rows; using the first.",
      "i" = "Rebuild with {.code build_taxon_db(dissolve = TRUE)} to union them."
    ))
  }
  db[hit[1], ]
}

#' Study extent: the range's bounding box plus a share of its diagonal
#'
#' Narrow endemics get a tight background, wide-ranging species a broad one.
#' Clamped to the globe so a high-latitude range does not ask for latitude 95.
#' @noRd
study_extent <- function(geom, buffer = 0.3) {
  bb <- sf::st_bbox(geom)
  diag <- sqrt((bb[["xmax"]] - bb[["xmin"]])^2 + (bb[["ymax"]] - bb[["ymin"]])^2)
  buf <- buffer * diag
  terra::ext(
    max(bb[["xmin"]] - buf, -180), min(bb[["xmax"]] + buf, 180),
    max(bb[["ymin"]] - buf, -90),  min(bb[["ymax"]] + buf, 90)
  )
}

#' Count the grid cells with climate that a polygon covers
#' @noRd
count_cells <- function(r, poly) {
  masked <- terra::mask(r[[1]], poly)
  n <- terra::global(!is.na(masked), "sum", na.rm = TRUE)[1, 1]
  if (is.na(n)) 0L else as.integer(n)
}

#' Count cells above a threshold
#' @noRd
count_above <- function(suit, cutoff) {
  n <- terra::global(suit > cutoff, "sum", na.rm = TRUE)[1, 1]
  if (is.na(n)) 0L else as.integer(n)
}

#' Sample climate values from distinct cells covered by a polygon
#'
#' Without replacement, so a range smaller than `n` cells contributes each of
#' its cells once rather than the same few repeatedly. terra warns "fewer cells
#' returned than requested" in that case; the caller reports it in
#' species-labelled form, so the raw warning is muffled here.
#' @noRd
sample_cells <- function(r, poly, n) {
  as.data.frame(
    withCallingHandlers(
      terra::spatSample(terra::mask(r, poly), n, method = "random",
                        replace = FALSE, na.rm = TRUE),
      warning = function(w) {
        if (grepl("fewer|requested|available", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
}
