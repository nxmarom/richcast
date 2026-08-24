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

#' Validate a c(xmin, xmax, ymin, ymax) extent
#'
#' The commonest mistake is supplying bbox order, `c(xmin, ymin, xmax, ymax)`,
#' which is what `sf::st_bbox()` prints and what most GIS tools show. That
#' misordering usually lands a longitude in a latitude slot, so the check tests
#' for it explicitly rather than leaving terra to report an "invalid extent".
#'
#' @param x Numeric vector of length 4.
#' @param arg Argument name to quote in the message.
#' @return `x`, invisibly.
#' @noRd
check_extent <- function(x, arg = "extent") {
  if (!is.numeric(x) || length(x) != 4) {
    rc_abort("{.arg {arg}} must be numeric of length 4: c(xmin, xmax, ymin, ymax).")
  }
  if (anyNA(x)) rc_abort("{.arg {arg}} must not contain missing values.")

  # If reading the input as bbox order c(xmin, ymin, xmax, ymax) yields a valid
  # extent, that is almost certainly what was meant -- so suggest the fix
  # rather than only reporting the failure.
  swapped <- x[c(1, 3, 2, 4)]
  looks_like_bbox <- swapped[1] < swapped[2] && swapped[3] < swapped[4] &&
    swapped[3] >= -90 && swapped[4] <= 90

  if (x[1] >= x[2] || x[3] >= x[4]) {
    rc_abort(c(
      "{.arg {arg}} must be c(xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax.",
      "x" = "Got xmin={x[1]}, xmax={x[2]}, ymin={x[3]}, ymax={x[4]}.",
      if (looks_like_bbox)
        c("i" = "This looks like bbox order c(xmin, ymin, xmax, ymax). Did you mean c({swapped[1]}, {swapped[2]}, {swapped[3]}, {swapped[4]})?")
    ))
  }
  if (x[3] < -90 || x[4] > 90) {
    rc_abort(c(
      "{.arg {arg}} has latitudes outside [-90, 90]: ymin={x[3]}, ymax={x[4]}.",
      "i" = "The order is c(xmin, xmax, ymin, ymax), not c(xmin, ymin, xmax, ymax)."
    ))
  }
  invisible(x)
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
  # A non-scalar here is almost always R's partial matching quietly redirecting
  # a `species =` argument onto `species_col =` against an older installed
  # version that has no `species` argument. Caught here it names the cause;
  # left alone it surfaces as "the condition has length > 1" from the %in%.
  if (!is.character(col) || length(col) != 1L || is.na(col)) {
    rc_abort(c(
      "{.arg {arg}} must be a single column name.",
      "x" = "Got {.cls {class(col)[1]}} of length {length(col)}.",
      if (length(col) > 1) c(
        "i" = "Passing several names looks like a species filter. That is {.arg species}, not {.arg {arg}} -- and R partial-matches {.code species=} onto {.arg species_col} on versions without it, so check that richcast is up to date."
      )
    ))
  }
  if (!col %in% names(data)) {
    rc_abort(c(
      "Column {.val {col}} not found.",
      "i" = "Available columns: {.val {names(data)}}.",
      "x" = "Set {.arg {arg}} to one of these."
    ))
  }
  invisible(TRUE)
}
