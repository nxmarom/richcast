#' Normalise species names to a canonical form
#'
#' Range databases and trait tables rarely agree on binomial formatting, and a
#' silent mismatch here shows up much later as an unexplained drop in species
#' count. This applies one canonical form to both sides of the join:
#' whitespace-trimmed, internal whitespace collapsed to a single underscore,
#' and sentence case (`Genus_species`).
#'
#' @param x Character vector of species names.
#' @return A character vector of the same length.
#' @examples
#' normalise_species(c("Marmota  baibacina", "marmota_bobak", " Rattus rattus "))
#' @export
normalise_species <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]_]+", "_", x)
  # Sentence case: capitalise the genus, lower-case the rest.
  ifelse(
    is.na(x) | !nzchar(x),
    NA_character_,
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  )
}

#' Abort with a richcast-flavoured error
#'
#' `.envir` is passed explicitly so that `{}` interpolation resolves against
#' the calling function's locals rather than this wrapper's frame.
#' @noRd
rc_abort <- function(message, ..., .envir = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = .envir, call = .envir)
}

#' Convert between years BP and years CE
#'
#' richcast works in years CE throughout, because a single signed axis orders
#' correctly whether a study spans the Holocene or the last millennium. Deep
#' time is conventionally quoted in years before present, though, and "present"
#' in that convention means 1950 CE -- so palaeoclimate reconstructions land on
#' years like 1950, 950, -50, -1050 rather than round thousands.
#'
#' These convert, so the offset is written down once instead of being
#' rediscovered as an off-by-1950 error.
#'
#' @param bp,ce Numeric vectors of years before present / years CE.
#' @param present Year CE that BP counts back from. 1950 by convention.
#' @return A numeric vector.
#' @examples
#' # The Last Glacial Maximum, roughly 20 ka BP:
#' bp_to_ce(-20000)
#'
#' # Beyer2020's millennial slices, as richcast wants them:
#' bp_to_ce(seq(0, -10000, by = -1000))
#'
#' ce_to_bp(1950)
#' @export
bp_to_ce <- function(bp, present = 1950) {
  present + bp
}

#' @rdname bp_to_ce
#' @export
ce_to_bp <- function(ce, present = 1950) {
  ce - present
}

#' Run a geometry operation, retrying in planar mode if s2 rejects the input
#'
#' Published range maps routinely contain rings that are valid enough for
#' planar tools but self-intersecting when interpreted on the sphere, and sf's
#' s2 backend refuses them outright ("Loop N is not valid"). `st_make_valid()`
#' does not reliably repair this, because planar validity and spherical
#' validity are different properties. Rather than failing on data that every
#' GIS will happily draw, spherical evaluation is attempted first and planar
#' evaluation used as a fallback.
#'
#' Takes a function rather than an expression so the retry genuinely
#' re-evaluates.
#'
#' @param f A zero-argument function performing the geometry operation.
#' @param what Short description used in the warning.
#' @param quiet Suppress the warning.
#' @noRd
with_planar_fallback <- function(f, what = "geometry operation", quiet = FALSE) {
  out <- tryCatch(f(), error = function(e) NULL)
  if (!is.null(out)) return(out)

  if (!quiet) {
    cli::cli_alert_warning(
      "Spherical {what} failed on invalid geometry; retrying in planar mode."
    )
  }
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  f()
}

#' Check that a column exists in a data frame, with a helpful message
#' @noRd
check_col <- function(data, col, arg = "species_col") {
  if (!col %in% names(data)) {
    rc_abort(c(
      "Column {.val {col}} not found.",
      "i" = "Available columns: {.val {names(data)}}.",
      "x" = "Set {.arg {arg}} to one of these."
    ))
  }
  invisible(TRUE)
}
