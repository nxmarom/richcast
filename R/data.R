#' Rodent ecological and zoonotic trait table
#'
#' Species-level attributes for 2000+ rodent species, from Supplementary Data 1
#' of Ecke et al. (2022). Bundled as a worked example of the trait table that
#' [build_taxon_db()] joins to a range source; there is nothing rodent-specific
#' about the rest of the package.
#'
#' @format A tibble with one row per species. Columns of most interest:
#' \describe{
#'   \item{species}{Binomial in `Genus_species` form (normalised by
#'     [normalise_species()]).}
#'   \item{S_index}{Synanthropy index. Note that every species present in this
#'     table scores above zero, so `S_index > 0` is a presence-in-table filter
#'     rather than a threshold.}
#'   \item{Synanthropic}{Binary synanthropy flag.}
#'   \item{Hunted}{Binary hunting-pressure flag.}
#'   \item{Reservoir}{Binary zoonotic-reservoir flag.}
#'   \item{NAm}{North America occurrence flag. Renamed from the source column
#'     `NA`, which collides with R's missing-value sentinel.}
#' }
#' Remaining columns carry habitat use, zoonosis counts by pathogen group,
#' population-dynamics measures, continental occurrence and MSW05 taxonomy.
#' See the source publication for full definitions.
#'
#' @source Ecke, F., Han, B. A., Hornfeldt, B., Khalil, H., Magnusson, M.,
#'   Singh, N. J., & Ostfeld, R. S. (2022). Population fluctuations and
#'   synanthropy explain transmission risk in rodent-borne zoonoses.
#'   *Nature Communications*, 13, 7532.
#'   \doi{10.1038/s41467-022-35273-7}. Redistributed under CC BY 4.0.
#'
#' @examples
#' head(rodent_traits[, c("species", "S_index", "Synanthropic", "Hunted")])
#' attr(rodent_traits, "source")
"rodent_traits"
