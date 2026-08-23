# Small synthetic fixtures so tests need no network and no third-party data.

square <- function(xmin, ymin, size = 1) {
  sf::st_polygon(list(cbind(
    c(xmin, xmin + size, xmin + size, xmin, xmin),
    c(ymin, ymin, ymin + size, ymin + size, ymin)
  )))
}

# Genus_alpha is deliberately split across two disjunct features, mimicking the
# way IUCN stores subspecies and seasonal ranges.
fixture_ranges <- function() {
  sf::st_sf(
    SCI_NAME = c("Genus alpha", "Genus alpha", "Genus beta", "Genus gamma"),
    geometry = sf::st_sfc(
      square(0, 0), square(5, 0), square(2, 2), square(8, 8),
      crs = 4326
    )
  )
}

fixture_traits <- function() {
  tibble::tibble(
    species = c("Genus_alpha", "Genus_beta"),
    mass = c(10, 500),
    score = c(0.8, 0.2)
  )
}

# --- Fixtures for end-to-end model runs -------------------------------------
# Land covering the whole test extent, so the Natural Earth mask does not
# erase a synthetic study area that sits in the Gulf of Guinea.
fake_land <- function() sf::st_sfc(square(-1, -1, 12), crs = 4326)

# A climate gradient on disk, one directory per slice, each slice nudged so
# projected ranges differ between them.
structured_climate <- function(dir, times, present = "present_1985") {
  grid <- function(shift) {
    r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10,
                     ymin = 0, ymax = 10, crs = "EPSG:4326")
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
    terra::values(r) <- xy[, 1] * 3 + xy[, 2] + shift
    r
  }
  labels <- c(present, sprintf("time_%04d", times))
  for (i in seq_along(labels)) {
    d <- file.path(dir, labels[i])
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    b1 <- grid(i * 0.5)
    b12 <- grid(-i * 0.5) * 2
    names(b1) <- "bio01"
    names(b12) <- "bio12"
    terra::writeRaster(b1, file.path(d, "bio01.tif"), overwrite = TRUE)
    terra::writeRaster(b12, file.path(d, "bio12.tif"), overwrite = TRUE)
  }
  climate_dir(dir, present = present)
}

# Two species at opposite ends of that gradient.
two_species_db <- function() {
  build_taxon_db(
    sf_polygons(sf::st_sf(
      species = c("Genus_low", "Genus_high"),
      geometry = sf::st_sfc(square(0.5, 0.5, 3), square(6, 6, 3), crs = 4326)
    )),
    quiet = TRUE
  )
}
