# ==============================================================================
# The whole pipeline, over a series of time slices
# ==============================================================================

#' Hindcast a whole assemblage across a series of time slices
#'
#' Selects the species near `region` with [species_near()], fits one ensemble
#' model per species with [fit_sdm()], projects each onto every requested time
#' slice, and stacks the results into richness surfaces for the region.
#'
#' Each species is fitted **once**, against present-day climate over its own
#' study extent, then projected repeatedly. Training is independent of the
#' region: the region decides only which species are trained and where
#' richness is summarised.
#'
#' Each slice gets two surfaces on the climate grid cropped to the region:
#' `richness`, the number of species above their own threshold in each cell,
#' and `expected`, the sum of the species' ensemble suitabilities -- a
#' threshold-free estimate of how many species a cell supports.
#'
#' @param db A `richcast_db` from [build_taxon_db()].
#' @param climate A climate source from [climate_dir()] or [pastclim_climate()].
#' @param times Numeric vector of years CE.
#' @param region The modelling region, from [region()].
#' @param species The taxa to model. `NULL` (the default) uses the
#'   Pleistocene zooarchaeological list for the preset `region`
#'   ([zooarch_taxa()]); otherwise a character vector of your own. Either way
#'   each species is narrowed by [species_near()]: a species more than 10
#'   degrees from `region` is never trained.
#' @param merge Taxa that count as one identification, as a named list:
#'   `list(Dama_sp = c("Dama_dama", "Dama_mesopotamica"))`. A merged taxon
#'   named in `species` has each member modelled separately, on its own range,
#'   and the outputs combined: the taxon is present in a cell when any member
#'   is, and its suitability there is that of the member furthest above its
#'   own threshold. Modelling the union of the ranges instead lets the larger
#'   range swamp the smaller one. `NULL` (the default) uses the merges the
#'   region's zooarchaeological list defines.
#' @param keep_surfaces Retain the richness surfaces, so maps can be drawn
#'   afterwards. Stored wrapped, so the result survives `saveRDS()`.
#' @param on_error `"warn"` skips a failing species and carries on; `"stop"`
#'   aborts the run.
#' @param quiet Suppress per-species progress.
#' @param ... Further arguments passed to [fit_sdm()], e.g. `predictors` or
#'   `n_presence`.
#' @return A `richcast_series`:
#'   * `richness`: one row per slice (`present` first) with mean, median and
#'     maximum richness over the region, and `mean_expected`.
#'   * `species`: one row per taxon per slice. `cells` counts suitable cells
#'     over the study extent (summed over members for a merged taxon),
#'     `region_cells` over the region, each with its change from the present.
#'   * `models`: one row per fitted species, with the taxon it belongs to, the
#'     fit diagnostics and the AUC and Boyce index of each ensemble member and
#'     of the ensemble.
#'   * `fits`: the fitted models, for [project_sdm()] and [richness_at()].
#'   * `taxa`: which fitted species make up each taxon.
#'   * `surfaces`: the richness surfaces; see [richness_surface()].
#' @seealso [fit_sdm()], [richness_at()], [richness_surface()]
#' @export
run_hindcast_series <- function(db,
                                climate,
                                times,
                                region,
                                species = NULL,
                                merge = NULL,
                                keep_surfaces = TRUE,
                                on_error = c("warn", "stop"),
                                quiet = FALSE,
                                ...) {

  on_error <- match.arg(on_error)
  check_region(region)
  if (!is.numeric(times) || length(times) == 0) {
    rc_abort("{.arg times} must be a non-empty numeric vector of years CE.")
  }
  # Chronological from the outset, so every lag runs forwards in time.
  times <- sort(unique(as.numeric(times)))
  keys <- c("present", as.character(times))

  # The species list defaults to the region's zooarchaeological record, and
  # the 10-degree rule is not optional: a list can only narrow it.
  species <- normalise_species(species %||% default_taxa(region))
  merge <- merge %||% zooarch_merges(region)
  names(merge) <- normalise_species(names(merge))
  taxa <- stats::setNames(lapply(species, function(t) {
    normalise_species(merge[[t]] %||% t)
  }), species)

  wanted <- unique(unlist(taxa))
  near <- species_near(db, region, quiet = quiet)
  missing <- setdiff(wanted, db$species)
  if (length(missing) > 0) {
    cli::cli_warn("Not in database, skipped: {.val {missing}}.")
  }
  far <- setdiff(intersect(wanted, db$species), near)
  if (length(far) > 0) {
    cli::cli_warn(c(
      "More than {near_distance} deg from {.emph {region$label}}, skipped: {.val {far}}.",
      "i" = "Species are only trained for regions near their present range."
    ))
  }
  taxa <- lapply(taxa, intersect, near)
  taxa <- taxa[lengths(taxa) > 0]
  if (length(taxa) == 0) {
    rc_abort(c(
      "No species to model.",
      "i" = "No requested range lies within {near_distance} deg of {.emph {region$label}}."
    ))
  }

  fits <- list()
  fitted_taxa <- list()
  stacks <- NULL
  counts <- list()

  if (!quiet) {
    cli::cli_progress_bar("Fitting and projecting", total = length(taxa),
                          .envir = environment())
  }

  for (tx in names(taxa)) {
    if (!quiet) cli::cli_progress_update(.envir = environment())

    # Each member is fitted and projected on its own range; the taxon's
    # projections are combined only when stacked.
    member_proj <- list()
    for (sp in taxa[[tx]]) {
      fitted <- try_step(
        fit_sdm(db, sp, climate, quiet = quiet, ...),
        what = "fit", species = sp, on_error = on_error
      )
      if (is.null(fitted)) next

      # The richness grid is the climate grid cropped to the region, read
      # once the first model says which variables exist.
      if (is.null(stacks)) {
        template <- climate_at(climate, "present", fitted$predictors[1],
                               region_ext(region))
        stacks <- stats::setNames(lapply(keys, function(k) new_stack(template)),
                                  keys)
      }

      projections <- list(present = list(suit = suitability(fitted),
                                         cells = fitted$present_cells))
      for (tt in times) {
        proj <- try_step(
          project_sdm(fitted, climate, tt, quiet = quiet),
          what = paste("project onto", tt), species = sp, on_error = on_error
        )
        if (is.null(proj)) next
        projections[[as.character(tt)]] <- list(suit = suitability(proj),
                                                cells = proj$cells)
      }
      fits[[sp]] <- fitted
      member_proj[[sp]] <- projections
    }
    if (length(member_proj) == 0) next
    fitted_taxa[[tx]] <- names(member_proj)

    for (k in keys) {
      have <- Filter(function(p) !is.null(p[[k]]), member_proj)
      if (length(have) == 0) next
      stacks[[k]] <- stack_add(
        stacks[[k]],
        lapply(have, function(p) p[[k]]$suit),
        vapply(names(have), function(sp) fits[[sp]]$threshold, numeric(1))
      )
      counts[[length(counts) + 1L]] <- tibble::tibble(
        species = tx,
        time = if (k == "present") NA_real_ else as.numeric(k),
        cells = sum(vapply(have, function(p) as.integer(p[[k]]$cells), integer(1))),
        region_cells = attr(stacks[[k]], "region_cells")
      )
    }
  }
  if (!quiet) cli::cli_progress_done(.envir = environment())

  if (length(fits) == 0) {
    rc_abort("No species could be modelled; nothing to summarise.")
  }

  # --- Richness per slice -------------------------------------------------
  surfaces <- list()
  richness <- lapply(keys, function(k) {
    s <- stack_finish(stacks[[k]])
    if (stacks[[k]]$n == 0) {
      cli::cli_warn(c(
        "No species contributed at {k}.",
        "x" = "Richness for this slice is missing because nothing was projected, not because the assemblage was empty.",
        "i" = "Check that climate rasters exist for {k}."
      ))
    }
    if (keep_surfaces) surfaces[[k]] <<- terra::wrap(s)
    richness_row(s, period = if (k == "present") "present" else paste0("hindcast_", k),
                 time = if (k == "present") NA_real_ else as.numeric(k),
                 n_species = stacks[[k]]$n)
  })
  richness <- dplyr::bind_rows(richness)

  # --- Per-taxon change ---------------------------------------------------
  all_counts <- dplyr::bind_rows(counts)
  base <- all_counts[is.na(all_counts$time), ]
  species_tbl <- all_counts[!is.na(all_counts$time), ] |>
    dplyr::left_join(
      dplyr::select(base, "species", present_cells = "cells",
                    region_present_cells = "region_cells"),
      by = "species"
    ) |>
    dplyr::mutate(
      delta_from_present = .data$cells - .data$present_cells,
      region_delta_from_present = .data$region_cells - .data$region_present_cells
    ) |>
    dplyr::arrange(.data$species, .data$time)

  models <- models_table(fits)
  taxon_of <- unlist(lapply(names(fitted_taxa), function(tx) {
    stats::setNames(rep(tx, length(fitted_taxa[[tx]])), fitted_taxa[[tx]])
  }))
  models <- tibble::add_column(models, taxon = unname(taxon_of[models$species]),
                               .after = "species")

  structure(
    list(
      richness = richness,
      species  = species_tbl,
      models   = models,
      fits     = fits,
      taxa     = fitted_taxa,
      surfaces = surfaces,
      times    = times,
      region   = region
    ),
    class = "richcast_series"
  )
}

#' One row per fitted model, with metrics spread wide
#' @noRd
models_table <- function(fits) {
  rows <- lapply(fits, function(f) {
    m <- f$metrics
    out <- tibble::tibble(
      species       = f$species,
      range_cells   = f$range_cells,
      n_presence    = f$n_presence,
      n_background  = f$n_background,
      threshold     = f$threshold,
      present_cells = f$present_cells
    )
    for (i in seq_len(nrow(m))) {
      out[[paste0("auc_", m$model[i])]] <- m$auc[i]
      out[[paste0("boyce_", m$model[i])]] <- m$boyce[i]
    }
    out
  })
  dplyr::bind_rows(rows)
}

#' Extract a richness surface from a series
#'
#' @param series A `richcast_series` from [run_hindcast_series()].
#' @param time Year CE, or `"present"`.
#' @param layer `"richness"` for the count of species above threshold,
#'   `"expected"` for the sum of suitabilities, or `"both"`.
#' @return A `SpatRaster`.
#' @seealso [richness_grid()], [run_hindcast_series()]
#' @export
richness_surface <- function(series, time, layer = c("richness", "expected", "both")) {
  layer <- match.arg(layer)
  check_surfaces(series)
  key <- as.character(time)
  if (!key %in% names(series$surfaces)) {
    rc_abort(c(
      "No surface for time {.val {time}}.",
      "i" = "Available: {.val {names(series$surfaces)}}."
    ))
  }
  r <- terra::unwrap(series$surfaces[[key]])
  if (layer == "both") r else r[[layer]]
}

#' Richness surfaces as a tidy data frame
#'
#' Long-format cell values for every slice, ready for a faceted map.
#'
#' @param series A `richcast_series`.
#' @param times Which slices to include. Defaults to every hindcast slice;
#'   include `"present"` to add the present day.
#' @param drop_na Drop cells no model covered.
#' @return A tibble with `x`, `y`, `time`, `richness` and `expected`.
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
  check_surfaces(series)
  keys <- as.character(times %||% series$times)
  missing <- setdiff(keys, names(series$surfaces))
  if (length(missing) > 0) {
    rc_abort("No surface for {.val {missing}}.")
  }
  out <- dplyr::bind_rows(lapply(keys, function(k) {
    d <- terra::as.data.frame(terra::unwrap(series$surfaces[[k]]), xy = TRUE,
                              na.rm = drop_na)
    d$time <- k
    tibble::as_tibble(d)
  }))
  # Keep chronology numeric where possible; "present" forces character.
  if (!any(keys == "present")) out$time <- as.numeric(out$time)
  out[, c("x", "y", "time", "richness", "expected")]
}

#' Predict richness, and the species expected, at coordinates
#'
#' Evaluates every fitted model at the given points, so the answer is exact at
#' the climate grid's resolution and needs no stored surfaces. A species counts
#' towards `richness` where its ensemble suitability is above its own
#' threshold; `expected_richness` sums the suitabilities instead. The
#' `species` table lists every modelled species at each point, most suitable
#' first, which is the expected species list read straight off the suitability
#' surface.
#'
#' A point outside a species' study extent is outside anything its model can
#' speak to, so that species has `NA` suitability there and is not counted.
#'
#' @param x A `richcast_series`, a single `richcast_sdm`, or a list of them.
#' @param lon,lat Coordinates in decimal degrees, recycled to a common length.
#' @param time Year(s) CE, and/or `"present"`.
#' @param climate The climate source the models were fitted with.
#' @return A `richcast_point` with two tibbles: `richness` (one row per point
#'   per time) and `species` (one row per point, time and species).
#' @seealso [run_hindcast_series()]
#' @examples
#' \dontrun{
#' pt <- richness_at(res, lon = 35.5, lat = 33, time = c("present", 850),
#'                   climate = clim)
#' pt$richness
#' subset(pt$species, present)
#' }
#' @export
richness_at <- function(x, lon, lat, time = "present", climate) {
  fits <- as_fits(x)
  if (!is.numeric(lon) || !is.numeric(lat) || length(lon) == 0 || length(lat) == 0) {
    rc_abort("{.arg lon} and {.arg lat} must be numeric.")
  }
  n <- max(length(lon), length(lat))
  lon <- rep_len(lon, n)
  lat <- rep_len(lat, n)
  xy <- cbind(lon, lat)

  vars <- unique(unlist(lapply(fits, function(f) f$predictors)))
  pad <- 1
  ext <- terra::ext(min(lon) - pad, max(lon) + pad, min(lat) - pad, max(lat) + pad)

  pts <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"),
                      crs = 4326)
  on_land <- with_planar_fallback(
    function() {
      suppressMessages(lengths(sf::st_intersects(pts, fits[[1]]$land)) > 0)
    },
    what = "land test", quiet = TRUE
  )

  sp_rows <- list()
  for (tt in as.character(time)) {
    clim <- climate_at(climate, if (tt == "present") "present" else as.numeric(tt),
                       vars, ext)
    vals <- terra::extract(clim, xy)
    for (f in fits) {
      se <- f$study_extent
      inside <- lon >= se[1] & lon <= se[2] & lat >= se[3] & lat <= se[4] & on_land
      d <- vals[f$predictors]
      ok <- inside & stats::complete.cases(d)
      suit <- rep(NA_real_, n)
      if (any(ok)) suit[ok] <- predict_members(f$members, d[ok, , drop = FALSE])[, "ensemble"]
      sp_rows[[length(sp_rows) + 1L]] <- tibble::tibble(
        lon = lon, lat = lat, time = tt,
        species = f$species, suitability = suit, threshold = f$threshold,
        present = !is.na(suit) & suit > f$threshold
      )
    }
  }
  species_tbl <- dplyr::bind_rows(sp_rows) |>
    collapse_taxa(taxa_of(x), c("lon", "lat", "time")) |>
    dplyr::arrange(.data$time, .data$lon, .data$lat,
                   dplyr::desc(.data$suitability))

  richness_tbl <- species_tbl |>
    dplyr::group_by(.data$lon, .data$lat, .data$time) |>
    dplyr::summarise(
      richness = sum(.data$present),
      expected_richness = sum(.data$suitability, na.rm = TRUE),
      species_modelled = sum(!is.na(.data$suitability)),
      .groups = "drop"
    )

  structure(list(richness = richness_tbl, species = species_tbl),
            class = "richcast_point")
}

#' Predict richness, and the species expected, over a small focus area
#'
#' The area counterpart of [richness_at()]. A single climate cell is often
#' too small a unit to read an assemblage from, so this pools a handful of
#' cells -- say the 2 x 2 block covering one landscape -- and asks which
#' species the models place anywhere in it.
#'
#' Every grid cell whose centre lies inside `area` is evaluated. A species is
#' **present** in the area when its ensemble suitability clears its threshold
#' in at least one of those cells, and its `suitability` is the highest value
#' among them, so the two agree: present exactly when `suitability >
#' threshold`. `richness` counts the species present; `expected_richness`
#' sums their suitabilities. Cells without climate (sea, or ice) are skipped,
#' so an area on a coast can hold fewer cells in some slices than others --
#' `cells` records how many were evaluated.
#'
#' @param x A `richcast_series`, a single `richcast_sdm`, or a list of them.
#'   For a series, merged taxa are reported once, as in [run_hindcast_series()].
#' @param area The focus area, a [region()] box.
#' @param time Year(s) CE, and/or `"present"`.
#' @param climate The climate source the models were fitted with.
#' @return A `richcast_area` with two tibbles: `richness` (one row per time)
#'   and `species` (one row per time and species, with `suitability`,
#'   `threshold`, `present` and `cells_present`).
#' @seealso [richness_at()], [run_hindcast_series()]
#' @examples
#' \dontrun{
#' focus <- region(c(34.5, 35.5, 31.5, 32.5), label = "focus area")
#' fa <- richness_in(res, focus, time = c("present", -4050), climate = clim)
#' fa$richness
#' }
#' @export
richness_in <- function(x, area, time = "present", climate) {
  sp_rows <- list()
  area_rows <- list()
  for (tt in as.character(time)) {
    g <- suitability_grid(x, area, tt, climate)
    area_rows[[tt]] <- tibble::tibble(time = tt, cells = attr(g, "cells"))
    sp_rows[[tt]] <- g |>
      dplyr::group_by(.data$species, .data$threshold) |>
      dplyr::summarise(
        cells_present = sum(.data$suitability > .data$threshold, na.rm = TRUE),
        suitability = if (all(is.na(.data$suitability))) NA_real_ else
          max(.data$suitability, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(time = tt,
                    present = !is.na(.data$suitability) &
                      .data$suitability > .data$threshold) |>
      dplyr::select("time", "species", "suitability", "threshold", "present",
                    "cells_present")
  }
  species_tbl <- dplyr::bind_rows(sp_rows) |>
    collapse_taxa(taxa_of(x), "time") |>
    dplyr::arrange(.data$time, dplyr::desc(.data$suitability))

  richness_tbl <- species_tbl |>
    dplyr::group_by(.data$time) |>
    dplyr::summarise(
      richness = sum(.data$present),
      expected_richness = sum(.data$suitability, na.rm = TRUE),
      species_modelled = sum(!is.na(.data$suitability)),
      .groups = "drop"
    ) |>
    dplyr::left_join(dplyr::bind_rows(area_rows), by = "time")

  structure(list(richness = richness_tbl, species = species_tbl, area = area),
            class = "richcast_area")
}

#' Suitability of every fitted species in every cell of an area
#'
#' The cell-level table behind [richness_in()]: each grid cell whose centre
#' lies inside `area`, evaluated by every fitted model whose study extent
#' covers it. Cells without climate are skipped. Merged taxa are not
#' collapsed here -- each fitted species keeps its own rows -- so thresholds
#' other than the fitted ones can be applied afterwards.
#'
#' @inheritParams richness_in
#' @param time A single year CE, or `"present"`.
#' @return A tibble with `x`, `y`, `species`, `suitability` and `threshold`,
#'   one row per cell and fitted species. The number of cells with climate is
#'   attached as the `cells` attribute.
#' @seealso [richness_in()], [presence_thresholds()]
#' @export
suitability_grid <- function(x, area, time, climate) {
  fits <- as_fits(x)
  check_region(area, arg = "area")
  if (length(time) != 1) rc_abort("{.arg time} must be a single slice.")
  vars <- unique(unlist(lapply(fits, function(f) f$predictors)))
  tt <- as.character(time)

  clim <- climate_at(climate, if (tt == "present") "present" else as.numeric(tt),
                     vars, NULL)
  cells <- terra::cells(clim[[1]], region_ext(area))
  xy <- terra::xyFromCell(clim, cells)
  b <- area$box
  keep <- xy[, 1] > b[1] & xy[, 1] < b[2] & xy[, 2] > b[3] & xy[, 2] < b[4]
  xy <- xy[keep, , drop = FALSE]
  vals <- as.data.frame(terra::extract(clim, xy))
  has_climate <- stats::complete.cases(vals[vars])

  rows <- lapply(fits, function(f) {
    se <- f$study_extent
    inside <- xy[, 1] >= se[1] & xy[, 1] <= se[2] &
      xy[, 2] >= se[3] & xy[, 2] <= se[4]
    d <- vals[f$predictors]
    ok <- inside & stats::complete.cases(d)
    suit <- rep(NA_real_, nrow(xy))
    if (any(ok)) {
      suit[ok] <- predict_members(f$members, d[ok, , drop = FALSE])[, "ensemble"]
    }
    tibble::tibble(x = xy[, 1], y = xy[, 2], species = f$species,
                   suitability = suit, threshold = f$threshold)
  })
  out <- dplyr::bind_rows(rows)
  attr(out, "cells") <- sum(has_climate)
  out
}

#' Alternative presence thresholds for each fitted model
#'
#' richcast binarises with the p10 threshold, which tolerates 10% omission of
#' the training presences. Two alternatives, for reading presence as a band
#' rather than a line, computed from each model's present-day suitability
#' surface without refitting:
#'
#' * `min_presence`: the lowest suitability the model assigns anywhere inside
#'   the species' own range polygon. The most permissive rule; suitability
#'   between it and p10 reads as possible but marginal presence. It is set by
#'   the single worst cell in the range, so it is sensitive to ragged range
#'   edges.
#' * `tss`: the cutoff that maximises the true skill statistic (sensitivity +
#'   specificity - 1) when cells inside the range polygon are treated as
#'   presences and the rest of the study extent as absences. It depends on how
#'   wide the study extent is, and can fall above or below p10.
#'
#' @param x A `richcast_series`, a single `richcast_sdm`, or a list of them.
#' @param db The `richcast_db` the models were fitted from.
#' @return A tibble with `species`, `threshold` (the fitted p10 threshold),
#'   `min_presence` and `tss`.
#' @seealso [suitability_grid()]
#' @export
presence_thresholds <- function(x, db) {
  fits <- as_fits(x)
  dplyr::bind_rows(lapply(fits, function(f) {
    row <- db_row(db, f$species)
    poly <- terra::makeValid(terra::vect(sf::st_sf(geometry = sf::st_geometry(row))))
    suit <- suitability(f)
    pres <- terra::values(terra::mask(suit, poly), mat = FALSE, na.rm = TRUE)
    abs <- terra::values(terra::mask(suit, poly, inverse = TRUE), mat = FALSE,
                         na.rm = TRUE)
    tibble::tibble(species = f$species, threshold = f$threshold,
                   min_presence = if (length(pres)) min(pres) else NA_real_,
                   tss = tss_cutoff(pres, abs))
  }))
}

#' The cutoff maximising sensitivity + specificity - 1
#' @noRd
tss_cutoff <- function(pres, abs) {
  if (length(pres) == 0 || length(abs) == 0) return(NA_real_)
  cand <- unique(stats::quantile(c(pres, abs), seq(0, 1, length.out = 501),
                                 names = FALSE))
  tss <- vapply(cand, function(c) mean(pres >= c) + mean(abs < c) - 1,
                numeric(1))
  cand[which.max(tss)]
}

#' @export
print.richcast_area <- function(x, ...) {
  cli::cli_text("{.cls richcast_area} {.emph {x$area$label}}")
  r <- x$richness
  for (i in seq_len(nrow(r))) {
    sp <- x$species$species[x$species$present & x$species$time == r$time[i]]
    cli::cli_text(
      "  {r$time[i]}: richness {r$richness[i]}, expected {round(r$expected_richness[i], 2)} ({r$cells[i]} cell{?s})"
    )
    if (length(sp)) cli::cli_text("    {.emph {sp}}")
  }
  cli::cli_text("  {.code $richness} {.code $species}")
  invisible(x)
}

#' Which fitted species make up each taxon, if the object records it
#' @noRd
taxa_of <- function(x) {
  if (inherits(x, "richcast_series")) x$taxa else NULL
}

#' Collapse member species rows to one row per merged taxon
#'
#' Within each group the member furthest above its own threshold stands for
#' the taxon, so `present` stays exactly `suitability > threshold`. Rows of
#' species that are their own taxon pass through unchanged.
#' @noRd
collapse_taxa <- function(tbl, taxa, by) {
  if (is.null(taxa) ||
      all(vapply(names(taxa), function(t) identical(taxa[[t]], t), logical(1)))) {
    return(tbl)
  }
  lookup <- unlist(lapply(names(taxa), function(tx) {
    stats::setNames(rep(tx, length(taxa[[tx]])), taxa[[tx]])
  }))
  tbl$species <- dplyr::coalesce(unname(lookup[tbl$species]), tbl$species)
  tbl$margin <- tbl$suitability - tbl$threshold
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "species")))) |>
    dplyr::arrange(dplyr::desc(.data$margin), .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(-"margin")
}

#' Fitted models from a series, a single model, or a list of models
#' @noRd
as_fits <- function(x) {
  if (inherits(x, "richcast_series")) {
    x$fits
  } else if (inherits(x, "richcast_sdm")) {
    stats::setNames(list(x), x$species)
  } else if (is.list(x) && length(x) > 0 &&
             all(vapply(x, inherits, logical(1), "richcast_sdm"))) {
    x
  } else {
    rc_abort("{.arg x} must be a {.cls richcast_series}, a {.cls richcast_sdm}, or a list of them.")
  }
}

#' @export
print.richcast_point <- function(x, ...) {
  cli::cli_text("{.cls richcast_point}")
  r <- x$richness
  for (i in seq_len(nrow(r))) {
    sp <- x$species$species[x$species$present &
                              x$species$lon == r$lon[i] &
                              x$species$lat == r$lat[i] &
                              x$species$time == r$time[i]]
    cli::cli_text(
      "  ({r$lon[i]}, {r$lat[i]}) @ {r$time[i]}: richness {r$richness[i]}, expected {round(r$expected_richness[i], 2)}"
    )
    if (length(sp)) cli::cli_text("    {.emph {sp}}")
  }
  cli::cli_text("  {.code $richness} {.code $species}")
  invisible(x)
}

#' @export
print.richcast_series <- function(x, ...) {
  cli::cli_text("{.cls richcast_series}")
  cli::cli_text("  {length(x$fits)} species x {length(x$times)} slices ({min(x$times)} to {max(x$times)} CE)")
  cli::cli_text("  region: {x$region$label}")
  if (nrow(x$models) > 0) {
    cli::cli_text(
      "  median ensemble AUC {round(stats::median(x$models$auc_ensemble, na.rm = TRUE), 3)}, Boyce {round(stats::median(x$models$boyce_ensemble, na.rm = TRUE), 3)}"
    )
  }
  cli::cli_text("  {.code $richness} {.code $species} {.code $models} {.code $fits} {.code $surfaces}")
  invisible(x)
}


# ==============================================================================
# Internals
# ==============================================================================

#' @noRd
check_surfaces <- function(series) {
  if (!inherits(series, "richcast_series")) {
    rc_abort("{.arg series} must come from {.fn run_hindcast_series}.")
  }
  if (length(series$surfaces) == 0) {
    rc_abort(c(
      "This series holds no richness surfaces.",
      "i" = "Re-run with {.code keep_surfaces = TRUE}."
    ))
  }
  invisible(TRUE)
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
