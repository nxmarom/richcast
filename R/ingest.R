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

#' Range source: an IUCN Red List range shapefile
#'
#' Points at a Red List spatial download on your own disk. The data are **not**
#' redistributed by this package and are not bundled with it: the IUCN Red List
#' Terms of Use (v3, section 4) prohibit redistribution of Red List data, whole
#' or in part, including within derivative works. Download the ranges yourself
#' from <https://www.iucnredlist.org/resources/spatial-data-download>, accept
#' the terms, and cite the version you used.
#'
#' @param path Path to the shapefile (`.shp`) or geopackage.
#' @param species_col Name of the binomial column. IUCN exports normally call
#'   this `SCI_NAME`.
#' @param species Optional character vector of binomials. When given, only
#'   these are read, via an OGR attribute query, and names are matched in both
#'   `Genus species` and `Genus_species` form. Worth using when you want a
#'   handful of taxa from a continental download: the read itself is no faster
#'   on an unindexed shapefile, but nothing unwanted is materialised and the
#'   dissolve shrinks accordingly.
#' @param layer Optional layer name, for multi-layer sources. Defaults to the
#'   file name without its extension.
#' @return A `richcast_range_source`.
#' @seealso [sf_polygons()], [gbif_occurrences()], [build_taxon_db()]
#' @examples
#' src <- iucn_shapefile("~/iucn_rodentia/data_0.shp")
#' print(src)
#'
#' # Only the taxa you actually intend to model:
#' iucn_shapefile(
#'   "~/iucn_artiodactyla/data_0.shp",
#'   species = c("Gazella gazella", "Sus scrofa", "Capra ibex")
#' )
#' @export
iucn_shapefile <- function(path, species_col = "SCI_NAME", species = NULL,
                           layer = NULL) {
  new_range_source(
    "iucn_shapefile",
    path = path, species_col = species_col, species = species, layer = layer,
    resolver = function(params) {
      if (!file.exists(params$path)) {
        rc_abort(c(
          "Range shapefile not found at {.path {params$path}}.",
          "i" = paste(
            "IUCN Red List ranges are not bundled with richcast and must be",
            "downloaded separately from",
            "{.url https://www.iucnredlist.org/resources/spatial-data-download}."
          )
        ))
      }
      lyr <- params$layer %||% tools::file_path_sans_ext(basename(params$path))

      if (is.null(params$species)) {
        cli::cli_progress_step("Reading {.path {basename(params$path)}}")
        x <- if (is.null(params$layer)) {
          sf::st_read(params$path, quiet = TRUE)
        } else {
          sf::st_read(params$path, layer = lyr, quiet = TRUE)
        }
      } else {
        cli::cli_progress_step(
          "Reading {length(params$species)} species from {.path {basename(params$path)}}"
        )
        x <- read_species_subset(params$path, lyr, params$species_col,
                                 params$species)
      }

      check_col(x, params$species_col)
      if (nrow(x) == 0) {
        rc_abort(c(
          "No features matched.",
          "i" = "Check the names against the {.val {params$species_col}} column of the source."
        ))
      }
      x
    }
  )
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

#' Range source: GBIF occurrence records
#'
#' Fetches point occurrences from GBIF. Unlike the polygon backends this needs
#' no manual download and no licence negotiation, so it is the quickest way to
#' get a working pipeline for a new taxon. It is also the methodologically
#' cleaner input: [fit_sdm()] samples pseudo-occurrences from within polygon
#' ranges, whereas real occurrence points can be used directly.
#'
#' Requires the \pkg{rgbif} package.
#'
#' @param species Character vector of binomials to fetch.
#' @param limit Maximum records per species. GBIF caps a single query at 100000.
#' @param ... Further arguments passed to [rgbif::occ_search()], e.g.
#'   `year = "1950,2000"` or `country = "KG"`.
#' @return A `richcast_range_source`.
#' @examples
#' src <- gbif_occurrences(c("Marmota_baibacina", "Marmota_bobak"), limit = 500)
#' print(src)
#' @export
gbif_occurrences <- function(species, limit = 5000, ...) {
  new_range_source(
    "gbif_occurrences",
    species = species, limit = limit, dots = list(...),
    resolver = function(params) {
      rlang::check_installed("rgbif", "to fetch GBIF occurrences.")
      sp <- gsub("_", " ", params$species)
      cli::cli_progress_bar("Querying GBIF", total = length(sp))
      out <- lapply(seq_along(sp), function(i) {
        cli::cli_progress_update(id = NULL)
        res <- do.call(rgbif::occ_search, c(
          list(
            scientificName = sp[i],
            hasCoordinate = TRUE,
            hasGeospatialIssue = FALSE,
            limit = params$limit
          ),
          params$dots
        ))
        d <- res$data
        if (is.null(d) || nrow(d) == 0) {
          cli::cli_warn("No GBIF records for {.val {sp[i]}}.")
          return(NULL)
        }
        data.frame(
          species = params$species[i],
          decimalLongitude = d$decimalLongitude,
          decimalLatitude = d$decimalLatitude
        )
      })
      out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
      if (is.null(out) || nrow(out) == 0) {
        rc_abort("GBIF returned no usable records for any requested species.")
      }
      sf::st_as_sf(
        out,
        coords = c("decimalLongitude", "decimalLatitude"),
        crs = 4326
      )
    }
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
      "i" = "Build one with {.fn iucn_shapefile}, {.fn sf_polygons}, or {.fn gbif_occurrences}."
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
#' @param ranges A range source from [iucn_shapefile()], [sf_polygons()] or
#'   [gbif_occurrences()].
#' @param traits Optional data frame of species attributes, e.g.
#'   [rodent_traits]. Joined on normalised species name.
#' @param trait_species_col Name of the species column in `traits`.
#' @param dissolve Logical. Union multiple features of the same species into a
#'   single geometry. Defaults to `TRUE`; ignored for point sources.
#' @param crs Target coordinate reference system. Defaults to EPSG:4326.
#' @param quiet Suppress progress messages.
#' @return A `richcast_db`: an `sf` object with a `species` column first,
#'   trait columns next, and geometry last.
#' @seealso [iucn_shapefile()], [gbif_occurrences()], [rodent_traits]
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

  src <- if (inherits(ranges, "sf")) sf_polygons(ranges) else ranges
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

  is_point <- all(sf::st_geometry_type(x) %in% c("POINT", "MULTIPOINT"))

  # --- Dissolve multi-feature species ------------------------------------
  if (dissolve && !is_point) {
    n_species <- dplyr::n_distinct(x$species)
    if (n_species < n_in) {
      say("Dissolving {n_in} feature{?s} into {n_species} species range{?s}.")
      x <- dissolve_by_species(x, quiet = quiet)
    }
  } else if (is_point) {
    say("Point source: {n_in} record{?s} across {dplyr::n_distinct(x$species)} species.")
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
