# ==============================================================================
# Fossil presences
#
# A present-day range is what is left of a species' niche after the Holocene;
# where a species has since been extirpated, the range no longer shows the
# climates it once lived in. Dated fossil occurrences put those climates back:
# each is paired with the palaeoclimate of its own time slice and added to the
# model's presences, with time-matched background, so the model learns the
# niche rather than the difference between the present and the past.
# ==============================================================================

#' Read a table of dated fossil occurrences
#'
#' Puts fossil occurrences into the long form [fit_sdm()] and
#' [run_hindcast_series()] use: one row per dated unit and time slice, with
#' the probability that the unit falls in that slice. Called by both, so a
#' table only needs to come through here to be checked.
#'
#' Two input forms are accepted, each with `species`, `lon` and `lat`:
#'
#' * **Wide**: a `slices` column of `"ka_bp:probability"` pairs separated by
#'   semicolons, e.g. `"34:0.8;36:0.2"`, one row per dated unit.
#' * **Long**: `ka_bp` and `weight` columns, one row per unit and slice.
#'   Rows with the same `unit_id` (or, without one, the same species and
#'   coordinates) are one unit.
#'
#' # Input format
#'
#' One dated unit is one dated find, layer or cluster of layers at a site.
#' Its age uncertainty is given as the probability of each model slice
#' (`ka_bp`, in thousands of years before present, on the climate's slice
#' grid). In the **wide** form, each row is a unit:
#'
#' | species | unit_id | locality | lon | lat | slices |
#' |---|---|---|---|---|---|
#' | Cervus_elaphus | KEB-D | Kebara Cave | 34.94 | 32.56 | 36:0.3;38:0.4;40:0.3 |
#' | Cervus_elaphus | AMU-B2 | Amud Cave | 35.50 | 32.87 | 50:0.2;52:0.5;54:0.3 |
#' | Capreolus_capreolus | KEB-D | Kebara Cave | 34.94 | 32.56 | 36:0.3;38:0.4;40:0.3 |
#' | Capreolus_capreolus | QAF-IX | Qafzeh | 35.30 | 32.68 | 32:1 |
#'
#' The same table in the **long** form has one row per unit and slice:
#'
#' | species | unit_id | locality | lon | lat | ka_bp | weight |
#' |---|---|---|---|---|---|---|
#' | Cervus_elaphus | KEB-D | Kebara Cave | 34.94 | 32.56 | 36 | 0.3 |
#' | Cervus_elaphus | KEB-D | Kebara Cave | 34.94 | 32.56 | 38 | 0.4 |
#' | Cervus_elaphus | KEB-D | Kebara Cave | 34.94 | 32.56 | 40 | 0.3 |
#' | Cervus_elaphus | AMU-B2 | Amud Cave | 35.50 | 32.87 | 50 | 0.2 |
#' | ... | | | | | | |
#'
#' `species` must match the names in the taxon database (underscores or
#' spaces), and must name the fitted species, not a merged taxon. A unit
#' known to fall in a single slice gets one pair with probability 1 (Qafzeh
#' above); weights need not sum to 1. A layer dated by a range rather than by
#' dates can be spread evenly over the slices the range covers. The values in
#' these tables are made up.
#'
#' A unit's probabilities are rescaled to sum to 1. An optional `block`
#' column, `"lon_lat"` of a box centre, gives the 1 x 1 degree box searched
#' for the nearest cell with climate when the site's own cell has none in a
#' slice; without it, the eight cells around the site are searched. Other
#' columns, such as `locality`, are kept.
#'
#' # Data sources
#'
#' richcast ships no fossil data. A natural source is the ROCEEH Out of
#' Africa Database (ROAD; Kandel et al. 2023), which can be queried from R
#' with the \pkg{roadDB} package (Kanaeva et al. 2026); the dated red deer
#' and roe deer occurrences used for the Middle East analysis were compiled
#' from ROAD that way. ROAD content is published under CC BY-SA 4.0: cite
#' Kandel et al. (2023) whenever you use it, and share any table derived
#' from it under the same licence. The table in the example below is made
#' up.
#'
#' @param x A data frame, or a path to a CSV file.
#' @return A tibble with `species`, `unit_id`, `lon`, `lat`, `block`, `ka_bp`
#'   and `prob`, plus any other input columns.
#' @references Kandel, A. W., Sommer, C., Kanaeva, Z., Bolus, M., Bruch, A.
#'   A., Groth, C., Haidle, M. N., Hertler, C., Heß, J., Malina, M., Märker,
#'   M., Hochschild, V., Mosbrugger, V., Schrenk, F., & Conard, N. J. (2023).
#'   The ROCEEH Out of Africa Database (ROAD): A large-scale research database
#'   serves as an indispensable tool for human evolutionary studies. *PLOS
#'   ONE*, 18(8), e0289513. \doi{10.1371/journal.pone.0289513}
#'
#'   Kanaeva, Z., Borre Pedersen, J., Sommer, C., & Streicher, T. P. (2026).
#'   *roadDB: Access Data from the ROCEEH Out of Africa Database (ROAD)*. R
#'   package version 0.2.0. \doi{10.32614/CRAN.package.roadDB}
#' @seealso [fit_sdm()] ("Fossil presences")
#' @examples
#' # Wide form: one row per dated unit (made-up values)
#' wide <- data.frame(
#'   species  = c("Cervus_elaphus", "Cervus_elaphus", "Capreolus_capreolus"),
#'   unit_id  = c("KEB-D", "AMU-B2", "QAF-IX"),
#'   locality = c("Kebara Cave", "Amud Cave", "Qafzeh"),
#'   lon      = c(34.94, 35.50, 35.30),
#'   lat      = c(32.56, 32.87, 32.68),
#'   slices   = c("36:0.3;38:0.4;40:0.3", "50:0.2;52:0.5;54:0.3", "32:1")
#' )
#' fossil_presences(wide)
#'
#' # Long form: one row per unit and slice (the Kebara and Qafzeh units above)
#' long <- data.frame(
#'   species  = c(rep("Cervus_elaphus", 3), "Capreolus_capreolus"),
#'   unit_id  = c("KEB-D", "KEB-D", "KEB-D", "QAF-IX"),
#'   lon      = c(34.94, 34.94, 34.94, 35.30),
#'   lat      = c(32.56, 32.56, 32.56, 32.68),
#'   ka_bp    = c(36, 38, 40, 32),
#'   weight   = c(0.3, 0.4, 0.3, 1)
#' )
#' fossil_presences(long)
#'
#' # A CSV in either form works too:
#' # fossils <- fossil_presences("my_fossils.csv")
#' # fit_sdm(db, "Cervus_elaphus", clim, fossils = fossils)
#' @export
fossil_presences <- function(x) {
  if (is.character(x) && length(x) == 1) x <- utils::read.csv(x, stringsAsFactors = FALSE)
  if (inherits(x, "richcast_fossils")) return(x)
  if (!is.data.frame(x)) rc_abort("{.arg fossils} must be a data frame or a CSV path.")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  need_cols(x, c("species", "lon", "lat"))
  if (!is.numeric(x$lon) || !is.numeric(x$lat) || anyNA(x$lon) || anyNA(x$lat)) {
    rc_abort("{.arg fossils} needs numeric, non-missing {.field lon} and {.field lat}.")
  }
  x$species <- normalise_species(x$species)
  if (is.null(x$unit_id)) {
    x$unit_id <- if (!is.null(x$slices)) sprintf("F%03d", seq_len(nrow(x))) else
      paste(x$species, x$lon, x$lat, sep = "@")
  }
  x$unit_id <- as.character(x$unit_id)
  if (is.null(x$block)) x$block <- NA_character_

  if (!is.null(x$slices)) {
    parts <- strsplit(gsub("\\s", "", as.character(x$slices)), ";", fixed = TRUE)
    rows <- lapply(seq_len(nrow(x)), function(i) {
      kv <- do.call(rbind, strsplit(parts[[i]], ":", fixed = TRUE))
      if (is.null(kv) || ncol(kv) != 2) {
        rc_abort("Unit {.val {x$unit_id[i]}}: {.field slices} must read {.val ka:prob;ka:prob}.")
      }
      base <- x[rep(i, nrow(kv)), setdiff(names(x), "slices"), drop = FALSE]
      base$ka_bp <- as.numeric(kv[, 1])
      base$weight <- as.numeric(kv[, 2])
      base
    })
    x <- do.call(rbind, rows)
  }
  if (is.null(x$weight) && !is.null(x$prob)) x$weight <- x$prob
  need_cols(x, c("ka_bp", "weight"))
  if (anyNA(x$ka_bp) || anyNA(x$weight) || any(x$weight < 0)) {
    rc_abort("{.arg fossils}: every slice needs a numeric age and a non-negative probability.")
  }
  x <- x[x$weight > 0, , drop = FALSE]
  key <- paste(x$species, x$unit_id, sep = "\r")
  x$prob <- x$weight / stats::ave(x$weight, key, FUN = sum)
  x$weight <- NULL
  rownames(x) <- NULL
  first <- c("species", "unit_id", "lon", "lat", "block", "ka_bp", "prob")
  out <- tibble::as_tibble(x[c(first, setdiff(names(x), first))])
  class(out) <- c("richcast_fossils", class(out))
  out
}

#' Draw fossil presences and their time-matched background
#'
#' `n` draws are spread evenly over the units (each unit gets `n %/% units`,
#' and a random few one more); each draw picks a slice by the unit's slice
#' probabilities and takes that slice's climate at the site. Background,
#' `n_background` points in all, is drawn uniformly over land in the study
#' extent from the same slices, in proportion to the draws in each.
#' @return A list of `presence` and `background` data frames of predictors,
#'   and `draws`, one row per presence with its unit and slice.
#' @noRd
draw_fossils <- function(fossils, climate, predictors, study_ext, land_vec,
                         n, n_background, climate_times = NULL) {
  units <- unique(fossils$unit_id)
  per_unit <- rep(n %/% length(units), length(units))
  extra <- n - sum(per_unit)
  if (extra > 0) {
    bump <- sample.int(length(units), extra)
    per_unit[bump] <- per_unit[bump] + 1L
  }
  draws <- do.call(rbind, lapply(seq_along(units), function(i) {
    if (per_unit[i] == 0) return(NULL)
    u <- fossils[fossils$unit_id == units[i], , drop = FALSE]
    pick <- if (nrow(u) == 1) rep(1L, per_unit[i]) else
      sample.int(nrow(u), per_unit[i], replace = TRUE, prob = u$prob)
    u[pick, c("unit_id", "lon", "lat", "block", "ka_bp"), drop = FALSE]
  }))
  draws$time <- ka_to_time(draws$ka_bp)
  if (!is.null(climate_times)) {
    gone <- setdiff(draws$time[draws$time != "present"], as.character(climate_times))
    if (length(gone) > 0) {
      rc_abort("No climate slice for fossil age{?s} {.val {unique(draws$ka_bp[draws$time %in% gone])}} ka.")
    }
  }

  slices <- table(draws$time)
  bg_n <- round(n_background * as.numeric(slices) / sum(slices))
  names(bg_n) <- names(slices)
  study_vec <- terra::vect(study_ext, crs = "EPSG:4326")
  bg_area <- terra::intersect(study_vec, land_vec)

  pres <- vector("list", length(slices))
  bg <- vector("list", length(slices))
  for (j in seq_along(slices)) {
    tt <- names(slices)[j]
    r <- climate_at(climate, if (tt == "present") "present" else as.numeric(tt),
                    predictors, study_ext)[[predictors]]
    here <- which(draws$time == tt)
    pres[[j]] <- site_climate(r, draws[here, , drop = FALSE])
    if (bg_n[[tt]] > 0) bg[[j]] <- sample_cells(r, bg_area, bg_n[[tt]])
  }
  pres <- do.call(rbind, pres)
  ok <- stats::complete.cases(pres[predictors])
  if (!all(ok)) {
    cli::cli_warn("{sum(!ok)} fossil draw{?s} had no climate within reach of the site and {?was/were} dropped.")
  }
  list(presence = pres[ok, predictors, drop = FALSE],
       background = do.call(rbind, bg)[predictors],
       draws = pres[ok, setdiff(names(pres), predictors), drop = FALSE])
}

#' Climate at fossil sites, falling back to the nearest cell with climate
#'
#' Within the unit's `block` box (centre +/- 0.5 degrees), or else the eight
#' cells around the site.
#' @noRd
site_climate <- function(r, draws) {
  xy <- cbind(draws$lon, draws$lat)
  vals <- terra::extract(r, xy)
  vals <- vals[names(r)]
  miss <- which(!stats::complete.cases(vals))
  for (i in miss) {
    b <- draws$block[i]
    if (!is.na(b) && nzchar(b)) {
      cxy <- as.numeric(strsplit(b, "_", fixed = TRUE)[[1]])
      box <- terra::ext(cxy[1] - 0.5, cxy[1] + 0.5, cxy[2] - 0.5, cxy[2] + 0.5)
    } else {
      res <- terra::res(r)
      box <- terra::ext(xy[i, 1] - 1.5 * res[1], xy[i, 1] + 1.5 * res[1],
                        xy[i, 2] - 1.5 * res[2], xy[i, 2] + 1.5 * res[2])
    }
    cells <- terra::cells(r[[1]], box)
    if (length(cells) == 0) next
    cv <- terra::extract(r, cells)
    keep <- stats::complete.cases(cv)
    if (!any(keep)) next
    cxy <- terra::xyFromCell(r, cells[keep])
    near <- which.min((cxy[, 1] - xy[i, 1])^2 + (cxy[, 2] - xy[i, 2])^2)
    vals[i, ] <- cv[keep, , drop = FALSE][near, ]
  }
  cbind(as.data.frame(vals), draws[c("unit_id", "ka_bp", "time")])
}

#' A slice age in ka BP as the time key climate sources use
#' @noRd
ka_to_time <- function(ka) {
  ifelse(ka == 0, "present", format(bp_to_ce(-ka * 1000), scientific = FALSE, trim = TRUE))
}

#' Abort unless a fossil table has the columns it needs
#' @noRd
need_cols <- function(x, cols) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    rc_abort(c(
      "{.arg fossils} lacks column{?s} {.field {missing}}.",
      "i" = "Give either {.field slices} or {.field ka_bp} and {.field weight}, alongside {.field species}, {.field lon} and {.field lat}."
    ))
  }
}
