# ==============================================================================
# Range sources
#
# A range source is a lazy description of where occurrence geometry comes from.
# It is resolved to an sf object by resolve_ranges(). Keeping resolution lazy
# means build_taxon_db() can report what it is about to read before spending
# minutes on a 170 MB shapefile.
# ==============================================================================

new_range_source <- function(type, ..., resolver) {
  structure(
    list(type = type, params = list(...), resolver = resolver),
    class = "richcast_range_source"
  )
}

#' @export
print.richcast_range_source <- function(x, ...) {
  cli::cli_text("{.cls richcast_range_source} <{x$type}>")
  for (nm in names(x$params)) {
    val <- x$params[[nm]]
    if (is.atomic(val) && length(val) <= 4) {
      cli::cli_text("  {nm}: {.val {val}}")
    } else {
      cli::cli_text("  {nm}: <{class(val)[1]}[{length(val)}]>")
    }
  }
  invisible(x)
}

#' Range source: IUCN Red List range polygons
#'
#' `iucn_folder()` reads every shapefile (and geopackage) found anywhere under
#' a folder, so a set of Red List downloads unpacked side by side -- say
#' `bovid_IUCN/`, `cervid_IUCN/` and `equid_IUCN/` -- is read as one source.
#' `iucn_shapefile()` reads a single file.
#'
#' The data are **not** redistributed by this package and are not bundled with
#' it: the IUCN Red List Terms of Use (v3, section 4) prohibit redistribution of
#' Red List data, whole or in part, including within derivative works. Download
#' the ranges yourself from
#' <https://www.iucnredlist.org/resources/spatial-data-download>, accept the
#' terms, and cite the version you used.
#'
#' # Which polygons count as today's range
#'
#' Red List files carry historical and uncertain range alongside the current
#' one. By default only polygons coded as extant (`PRESENCE` 1-3: extant,
#' probably extant, possibly extant) and native or reintroduced (`ORIGIN` 1-2)
#' are kept, since an extinct or introduced patch is not part of the climate
#' niche a model should learn. Set either argument to `NULL` to keep every
#' code. Files without the column are not filtered on it.
#'
#' @param path For `iucn_folder()`, a folder searched recursively. For
#'   `iucn_shapefile()`, a `.shp` or `.gpkg` file.
#' @param species_col Name of the binomial column. IUCN exports call this
#'   `SCI_NAME`.
#' @param species Optional character vector of binomials. When given, only
#'   these are read, via an OGR attribute query; names match in both
#'   `Genus species` and `Genus_species` form.
#' @param presence,origin IUCN `PRESENCE` and `ORIGIN` codes to keep, or `NULL`
#'   for all. See details.
#' @param layer Optional layer name, for multi-layer sources.
#' @return A `richcast_range_source`.
#' @seealso [build_taxon_db()], [sf_polygons()]
#' @examples
#' src <- iucn_folder("UngulatePolygons")
#' print(src)
#'
#' # Only the taxa you actually intend to model:
#' iucn_shapefile(
#'   "~/iucn_artiodactyla/data_0.shp",
#'   species = c("Gazella gazella", "Sus scrofa", "Capra ibex")
#' )
#' @export
iucn_folder <- function(path, species_col = "SCI_NAME", species = NULL,
                        presence = 1:3, origin = 1:2) {
  new_range_source(
    "iucn_folder",
    path = path, species_col = species_col, species = species,
    presence = presence, origin = origin,
    resolver = function(params) {
      if (!dir.exists(params$path)) {
        rc_abort(c(
          "Range folder not found at {.path {params$path}}.",
          "i" = iucn_download_hint()
        ))
      }
      files <- list.files(params$path, pattern = "\\.(shp|gpkg)$",
                          recursive = TRUE, full.names = TRUE,
                          ignore.case = TRUE)
      if (length(files) == 0) {
        rc_abort("No {.file .shp} or {.file .gpkg} files under {.path {params$path}}.")
      }
      cli::cli_progress_step(
        "Reading {length(files)} range file{?s} from {.path {params$path}}"
      )
      parts <- lapply(files, function(f) {
        x <- read_iucn_file(f, params$species_col, params$species,
                            layer = NULL, must_match = FALSE)
        x <- filter_iucn_codes(x, params$presence, params$origin)
        # Downloads for different groups do not always share columns; the
        # species column and the geometry are all that is needed downstream.
        x <- x[, params$species_col]
        sf::st_geometry(x) <- "geometry"
        if (is.na(sf::st_crs(x)) || sf::st_crs(x) == sf::st_crs(4326)) x
        else sf::st_transform(x, 4326)
      })
      x <- do.call(rbind, parts)
      if (nrow(x) == 0) {
        rc_abort(c(
          "No features matched.",
          "i" = "Check the names against the {.val {params$species_col}} column, and the {.arg presence}/{.arg origin} filters."
        ))
      }
      x
    }
  )
}

#' @rdname iucn_folder
#' @export
iucn_shapefile <- function(path, species_col = "SCI_NAME", species = NULL,
                           presence = 1:3, origin = 1:2, layer = NULL) {
  new_range_source(
    "iucn_shapefile",
    path = path, species_col = species_col, species = species,
    presence = presence, origin = origin, layer = layer,
    resolver = function(params) {
      if (!file.exists(params$path)) {
        rc_abort(c(
          "Range shapefile not found at {.path {params$path}}.",
          "i" = iucn_download_hint()
        ))
      }
      cli::cli_progress_step("Reading {.path {basename(params$path)}}")
      x <- read_iucn_file(params$path, params$species_col, params$species,
                          params$layer, must_match = TRUE)
      filter_iucn_codes(x, params$presence, params$origin)
    }
  )
}

#' @noRd
iucn_download_hint <- function() {
  paste(
    "IUCN Red List ranges are not bundled with richcast and must be",
    "downloaded separately from",
    "{.url https://www.iucnredlist.org/resources/spatial-data-download}."
  )
}

#' Read one IUCN file, optionally only some species
#' @noRd
read_iucn_file <- function(path, species_col, species, layer, must_match) {
  lyr <- layer %||% tools::file_path_sans_ext(basename(path))
  x <- if (is.null(species)) {
    if (is.null(layer)) sf::st_read(path, quiet = TRUE)
    else sf::st_read(path, layer = lyr, quiet = TRUE)
  } else {
    read_species_subset(path, lyr, species_col, species)
  }
  check_col(x, species_col)
  if (must_match && nrow(x) == 0) {
    rc_abort(c(
      "No features matched.",
      "i" = "Check the names against the {.val {species_col}} column of the source."
    ))
  }
  x
}

#' Keep only the IUCN presence and origin codes asked for
#' @noRd
filter_iucn_codes <- function(x, presence, origin) {
  if (!is.null(presence) && "PRESENCE" %in% names(x)) {
    x <- x[x$PRESENCE %in% presence, ]
  }
  if (!is.null(origin) && "ORIGIN" %in% names(x)) {
    x <- x[x$ORIGIN %in% origin, ]
  }
  x
}

#' Read only the requested species from a vector source
#'
#' Uses an OGR SQL predicate so unwanted features are never materialised in R.
#' On an unindexed shapefile this is not faster to *read* -- OGR scans every
#' feature either way, and evaluating the predicate costs a little extra -- but
#' it keeps memory down and, more importantly, spares the dissolve: unioning 8
#' species is a different proposition from unioning 322.
#'
#' Names are matched in both spaced and underscored form, since range databases
#' and trait tables disagree about which they use.
#'
#' @param path,layer Source and layer name.
#' @param species_col Column holding binomials in the source.
#' @param species Character vector of species to keep.
#' @return An `sf` object.
#' @noRd
read_species_subset <- function(path, layer, species_col, species) {

  # Match either spelling, and escape quotes rather than trusting the input.
  variants <- unique(c(gsub("_", " ", species), gsub(" ", "_", species)))
  variants <- variants[!is.na(variants) & nzchar(variants)]
  quoted <- paste0("'", gsub("'", "''", variants), "'")

  query <- sprintf('SELECT * FROM "%s" WHERE "%s" IN (%s)',
                   layer, species_col, paste(quoted, collapse = ", "))

  out <- tryCatch(
    sf::st_read(path, query = query, quiet = TRUE),
    error = function(e) NULL
  )

  if (is.null(out)) {
    # Not every driver supports SQL. Falling back is slower and heavier, but a
    # working slow path beats an error the user cannot act on.
    cli::cli_alert_warning(
      "Attribute query unsupported by this driver; reading in full and filtering."
    )
    out <- sf::st_read(path, quiet = TRUE)
    check_col(out, species_col)
    out <- out[normalise_species(out[[species_col]]) %in%
                 normalise_species(species), ]
  }
  out
}

#' Range source: arbitrary sf polygons
#'
#' The general-purpose backend. Use it for any polygon range map you already
#' hold in R -- a national atlas, expert-drawn ranges, alpha hulls, or the
#' output of someone else's SDM.
#'
#' @param x An `sf` object with polygon or multipolygon geometry.
#' @param species_col Name of the column holding species names.
#' @return A `richcast_range_source`.
#' @examples
#' ranges <- sf::st_sf(
#'   species = c("Genus_alpha", "Genus_beta"),
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(cbind(c(0, 2, 2, 0, 0), c(0, 0, 2, 2, 0)))),
#'     sf::st_polygon(list(cbind(c(1, 3, 3, 1, 1), c(1, 1, 3, 3, 1)))),
#'     crs = 4326
#'   )
#' )
#' sf_polygons(ranges)
#' @export
sf_polygons <- function(x, species_col = "species") {
  if (!inherits(x, "sf")) {
    rc_abort("{.arg x} must be an {.cls sf} object, not {.cls {class(x)[1]}}.")
  }
  check_col(x, species_col)
  new_range_source(
    "sf_polygons",
    x = x, species_col = species_col,
    resolver = function(params) params$x
  )
}

#' Union all features belonging to the same species
#'
#' Published range maps are full of self-intersecting rings and duplicated
#' vertices, which sf's default s2 backend rejects outright ("Loop N is not
#' valid"). Real IUCN downloads trip this reliably, so a spherical attempt is
#' followed by a planar retry rather than an error. Planar union is a slight
#' approximation at continental scale and across the antimeridian, but the
#' alternative is refusing to read the data at all.
#'
#' @param x An `sf` object with a `species` column.
#' @param quiet Suppress messages.
#' @return An `sf` object with one row per species.
#' @noRd
dissolve_by_species <- function(x, quiet = FALSE) {

  x <- suppressWarnings(sf::st_make_valid(x))

  out <- with_planar_fallback(
    function() {
      suppressWarnings(sf::st_make_valid(x)) |>
        dplyr::group_by(.data$species) |>
        dplyr::summarise(.groups = "drop")
    },
    what = "union", quiet = quiet
  )

  suppressWarnings(sf::st_make_valid(out))
}

#' Resolve a range source to an sf object
#' @param src A `richcast_range_source`.
#' @return An `sf` object.
#' @noRd
resolve_ranges <- function(src) {
  if (!inherits(src, "richcast_range_source")) {
    rc_abort(c(
      "{.arg ranges} must be a range source.",
      "i" = "Build one with {.fn iucn_folder}, {.fn iucn_shapefile} or {.fn sf_polygons}, or pass a folder path."
    ))
  }
  src$resolver(src$params)
}


# ==============================================================================
# Assembling the taxon database
# ==============================================================================

#' Assemble a taxon database from ranges and (optionally) traits
#'
#' Joins a resolved range source to a trait table and returns the object every
#' other richcast function consumes.
#'
#' Two behaviours are worth knowing about, because both silently corrupted
#' results in the pipeline this package was extracted from:
#'
#' * **Multi-row species are dissolved.** Range databases routinely store one
#'   species across many features -- IUCN splits by subspecies, seasonality and
#'   disjunct patches, so a global rodent download holds ~3100 features for
#'   ~2350 species. Treating each feature as a species means a model is trained
#'   against a single fragment while the rest of the species' true range is
#'   handed to the background sample, i.e. the model learns to discriminate the
#'   species from itself. Set `dissolve = FALSE` only if you genuinely want
#'   per-feature modelling.
#' * **The trait join is reported, not silent.** Species present in the ranges
#'   but absent from `traits` keep their geometry and receive `NA` traits. The
#'   count is printed, because an inner-join-shaped drop is easy to mistake for
#'   a real biological filter.
#'
#' @param ranges A range source from [iucn_folder()], [iucn_shapefile()] or
#'   [sf_polygons()]; a bare `sf` object; or a folder path, which is read with
#'   [iucn_folder()] defaults.
#' @param traits Optional data frame of species attributes, e.g.
#'   [rodent_traits]. Joined on normalised species name.
#' @param trait_species_col Name of the species column in `traits`.
#' @param dissolve Logical. Union multiple features of the same species into a
#'   single geometry. Defaults to `TRUE`.
#' @param crs Target coordinate reference system. Defaults to EPSG:4326.
#' @param quiet Suppress progress messages.
#' @return A `richcast_db`: an `sf` object with a `species` column first,
#'   trait columns next, and geometry last.
#' @seealso [iucn_folder()], [rodent_traits]
#' @examples
#' ranges <- sf::st_sf(
#'   species = c("Genus_alpha", "Genus_alpha", "Genus_beta"),
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
#'     sf::st_polygon(list(cbind(c(2, 3, 3, 2, 2), c(0, 0, 1, 1, 0)))),
#'     sf::st_polygon(list(cbind(c(1, 3, 3, 1, 1), c(1, 1, 3, 3, 1)))),
#'     crs = 4326
#'   )
#' )
#' # The two Genus_alpha features collapse into one range.
#' build_taxon_db(sf_polygons(ranges), quiet = TRUE)
#' @export
build_taxon_db <- function(ranges,
                           traits = NULL,
                           trait_species_col = "species",
                           dissolve = TRUE,
                           crs = 4326,
                           quiet = FALSE) {

  say <- function(...) if (!quiet) cli::cli_alert_info(...)

  src <- if (inherits(ranges, "sf")) {
    sf_polygons(ranges)
  } else if (is.character(ranges) && length(ranges) == 1) {
    iucn_folder(ranges)
  } else {
    ranges
  }
  x <- resolve_ranges(src)
  species_col <- src$params$species_col %||% "species"

  n_in <- nrow(x)

  # --- Canonical species column -----------------------------------------
  x <- x[!is.na(x[[species_col]]), ]
  x$species <- normalise_species(x[[species_col]])
  if (species_col != "species") x[[species_col]] <- NULL
  x <- x[, c("species", setdiff(names(x), c("species", attr(x, "sf_column"))))]

  # --- CRS ---------------------------------------------------------------
  if (is.na(sf::st_crs(x))) {
    cli::cli_warn("Range source has no CRS; assuming EPSG:{crs}.")
    sf::st_crs(x) <- crs
  } else if (sf::st_crs(x) != sf::st_crs(crs)) {
    say("Reprojecting ranges to EPSG:{crs}.")
    x <- sf::st_transform(x, crs)
  }

  if (!all(sf::st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"))) {
    rc_abort(c(
      "Ranges must be polygons.",
      "i" = "richcast samples pseudo-presences from inside each range polygon."
    ))
  }

  # --- Dissolve multi-feature species ------------------------------------
  if (dissolve) {
    n_species <- dplyr::n_distinct(x$species)
    if (n_species < n_in) {
      say("Dissolving {n_in} feature{?s} into {n_species} species range{?s}.")
      x <- dissolve_by_species(x, quiet = quiet)
    }
  }

  # --- Trait join --------------------------------------------------------
  if (!is.null(traits)) {
    check_col(traits, trait_species_col, arg = "trait_species_col")
    tr <- tibble::as_tibble(traits)
    tr[[trait_species_col]] <- normalise_species(tr[[trait_species_col]])
    tr <- dplyr::distinct(tr, .data[[trait_species_col]], .keep_all = TRUE)
    names(tr)[names(tr) == trait_species_col] <- "species"

    x <- dplyr::left_join(x, tr, by = "species")

    matched <- sum(x$species %in% tr$species)
    if (!quiet) {
      cli::cli_alert_info(
        "Trait join: {matched}/{nrow(x)} species matched, {nrow(x) - matched} with NA traits."
      )
      if (matched < nrow(x)) {
        cli::cli_alert_warning(
          "Filtering on a trait column will silently drop the {nrow(x) - matched} unmatched species."
        )
      }
    }
  }

  x <- sf::st_sf(x)
  x <- x[, c("species", setdiff(names(x), c("species", attr(x, "sf_column"))))]
  class(x) <- c("richcast_db", class(x))
  x
}

#' @export
print.richcast_db <- function(x, ...) {
  # Count by name rather than by arithmetic on ncol(): st_drop_geometry() keeps
  # the class but removes a column, which made a trait-less database report
  # "-1 trait columns".
  n_traits <- length(setdiff(names(x), c("species", attr(x, "sf_column"))))
  cli::cli_text(
    "{.cls richcast_db}: {nrow(x)} species, {n_traits} trait column{?s}"
  )
  NextMethod()
}

#' Filter a taxon database
#'
#' A thin, self-documenting wrapper around [dplyr::filter()] that reports how
#' many species each expression removed. Useful because trait filters often
#' remove far more (or far less) than intended -- for example `S_index > 0` on
#' [rodent_traits] looks like a risk threshold but every matched species scores
#' above zero, so it only ever removes species missing from the trait table.
#'
#' @param db A `richcast_db`.
#' @param ... Filter expressions, passed to [dplyr::filter()].
#' @param quiet Suppress the report.
#' @return A filtered `richcast_db`.
#' @examples
#' db <- build_taxon_db(
#'   sf_polygons(sf::st_sf(
#'     species = c("Genus_alpha", "Genus_beta"),
#'     mass = c(10, 500),
#'     geometry = sf::st_sfc(
#'       sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
#'       sf::st_polygon(list(cbind(c(2, 3, 3, 2, 2), c(0, 0, 1, 1, 0)))),
#'       crs = 4326
#'     )
#'   )),
#'   quiet = TRUE
#' )
#' filter_taxa(db, mass < 100)
#' @export
filter_taxa <- function(db, ..., quiet = FALSE) {
  before <- nrow(db)

  out <- withCallingHandlers(
    dplyr::filter(db, ...),
    error = function(e) {
      # A filter on a column that is not there almost always means the trait
      # table was never joined, so say that rather than leaving the user with
      # dplyr's "object not found".
      traits <- setdiff(names(db), c("species", attr(db, "sf_column")))
      rc_abort(c(
        "Filter failed: {conditionMessage(e)}",
        if (length(traits) == 0) c(
          "x" = "This database has no trait columns.",
          "i" = "Pass {.arg traits} to {.fn build_taxon_db} to join a trait table."
        ) else c(
          "i" = "Available trait columns: {.val {traits}}."
        )
      ))
    }
  )
  if (!quiet) {
    cli::cli_alert_info(
      "Filter kept {nrow(out)}/{before} species ({before - nrow(out)} removed)."
    )
  }
  out
}
