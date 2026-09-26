# ==============================================================================
# Pleistocene zooarchaeological species lists
#
# A hindcast should only place species where the Pleistocene record says they
# could plausibly have been. The default species set for a preset region is a
# list compiled from zooarchaeological and palaeontological reports from that
# region, stored in inst/extdata/zooarch_taxa.csv with the evidence and source
# for every entry, so each inclusion can be checked and argued with.
# ==============================================================================

#' Pleistocene zooarchaeological species lists
#'
#' The species richcast models by default for a preset [region()]: taxa
#' reported from Pleistocene zooarchaeological or palaeontological sites in
#' that region, restricted to species that have present-day IUCN range maps.
#' Extinct taxa (aurochs, *Equus hydruntinus*) and taxa outside the bundled
#' families have no range polygon to train on and are left out.
#'
#' # Middle East
#'
#' Sixteen bovids, cervids and equids from the Levant, the Zagros, Syria and
#' Arabia, compiled from:
#'
#' * Stewart, M., Louys, J., Price, G. J., Drake, N. A., Groucutt, H. S., &
#'   Petraglia, M. D. (2019). Middle and Late Pleistocene mammal fossils of
#'   Arabia and surrounding regions: implications for biogeography and hominin
#'   dispersals. *Quaternary International*, 515, 12-29.
#' * Stiner, M. C. (2005). *The Faunas of Hayonim Cave, Israel*. Peabody
#'   Museum, Harvard.
#' * Yeshurun, R., Bar-Oz, G., & Weinstein-Evron, M. (2007). Modern hunting
#'   behavior in the early Middle Paleolithic: faunal remains from Misliya
#'   Cave, Mount Carmel, Israel. *Journal of Human Evolution*, 53, 656-677.
#' * Stiner, M. C., Barkai, R., & Gopher, A. (2009). Cooperative hunting and
#'   meat sharing 400-200 kya at Qesem Cave, Israel. *PNAS*, 106,
#'   13207-13212.
#' * Tsahar, E., Izhaki, I., Lev-Yadun, S., & Bar-Oz, G. (2009). Distribution
#'   and extinction of ungulates during the Holocene of the southern Levant.
#'   *PLoS ONE*, 4, e5316.
#' * Mata-Gonzalez, M., Starkovich, B. M., Zeidi, M., & Conard, N. J. (2023).
#'   Evidence of diverse animal exploitation during the Middle Paleolithic at
#'   Ghar-e Boof (southern Zagros). *Scientific Reports*, 13, 19006.
#'
#' Older names are mapped to current IUCN taxonomy: *Capra ibex* from the
#' Negev to *Capra nubiana*, *Ovis orientalis* to *Ovis gmelini*, *Equus
#' caballus* to *Equus ferus*.
#'
#' Lists for the other presets are not yet compiled; pass your own species to
#' [run_hindcast_series()] for those regions.
#'
#' @param region A [region()] or a preset name. `NULL` returns every list.
#' @return A tibble with `region`, `species`, `evidence` and `source`.
#' @examples
#' zooarch_taxa("middle_east")
#' @export
zooarch_taxa <- function(region = NULL) {
  path <- system.file("extdata", "zooarch_taxa.csv", package = "richcast")
  tbl <- tibble::as_tibble(utils::read.csv(path, stringsAsFactors = FALSE))
  if (is.null(region)) return(tbl)
  key <- if (inherits(region, "richcast_region")) region$preset else
    region(region)$preset
  tbl[tbl$region %in% key, ]
}

#' The default species list for a region, or an error saying there is none
#' @noRd
default_taxa <- function(region) {
  tbl <- zooarch_taxa(region)
  if (nrow(tbl) == 0) {
    have <- unique(zooarch_taxa()$region)
    rc_abort(c(
      "No Pleistocene zooarchaeological species list is bundled for {.emph {region$label}}.",
      "i" = "Lists exist for: {.val {have}}.",
      "i" = "Pass your own list as {.arg species}."
    ))
  }
  tbl$species
}
