# ==============================================================================
# data-raw/rodent_traits.R
#
# Builds data/rodent_traits.rda from Supplementary Data 1 of:
#
#   Ecke, F., Han, B. A., Hornfeldt, B., Khalil, H., Magnusson, M.,
#   Singh, N. J., & Ostfeld, R. S. (2022). Population fluctuations and
#   synanthropy explain transmission risk in rodent-borne zoonoses.
#   Nature Communications, 13, 7532. https://doi.org/10.1038/s41467-022-35273-7
#
# Nature Communications publishes under CC BY 4.0, so this table is
# redistributable with attribution. See inst/CITATION and the `source` /
# `license` attributes attached below.
#
# NOTE: no IUCN Red List data is included here or anywhere else in the package.
# IUCN range polygons may not be redistributed (Terms of Use v3, section 4);
# users supply their own via iucn_shapefile().
# ==============================================================================

library(dplyr)
library(readr)

src <- "../Inputs/rodent_ecology.csv"  # relative to package root

raw <- readr::read_csv(src, show_col_types = FALSE, na = c("", "NA"))

# The source table has a column literally named "NA" (North America, alongside
# EU/AS/AF/AU/SA). Left alone it collides with the NA sentinel in almost every
# tidy verb, so it is renamed before anything else touches the data.
if ("NA" %in% names(raw)) {
  names(raw)[names(raw) == "NA"] <- "NAm"
}

rodent_traits <- raw |>
  dplyr::rename(species = "Rodent_species") |>
  dplyr::mutate(species = normalise_species(species)) |>
  dplyr::distinct(species, .keep_all = TRUE) |>
  dplyr::relocate("species", "S_index", "Synanthropic", "Hunted", "Reservoir") |>
  tibble::as_tibble()

attr(rodent_traits, "source") <- paste(
  "Ecke et al. (2022) Nature Communications 13:7532, Supplementary Data 1.",
  "https://doi.org/10.1038/s41467-022-35273-7"
)
attr(rodent_traits, "license") <- "CC BY 4.0"

usethis::use_data(rodent_traits, overwrite = TRUE, compress = "xz")
