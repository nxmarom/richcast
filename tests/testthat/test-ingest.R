test_that("normalise_species canonicalises binomials", {
  expect_equal(
    normalise_species(c("Marmota  baibacina", "marmota_bobak", " Rattus rattus ")),
    c("Marmota_baibacina", "Marmota_bobak", "Rattus_rattus")
  )
  expect_equal(normalise_species(c(NA, "")), c(NA_character_, NA_character_))
})

test_that("multi-feature species are dissolved into one range", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    quiet = TRUE
  )
  expect_s3_class(db, "richcast_db")
  expect_equal(nrow(db), 3)
  expect_setequal(db$species, c("Genus_alpha", "Genus_beta", "Genus_gamma"))

  # The dissolved Genus_alpha covers both disjunct squares, so its bounding
  # box must span them. Getting this wrong is the bug the dissolve step exists
  # to prevent: half the range would otherwise land in the background sample.
  bb <- sf::st_bbox(db[db$species == "Genus_alpha", ])
  expect_equal(as.numeric(bb["xmin"]), 0)
  expect_equal(as.numeric(bb["xmax"]), 6)
})

test_that("dissolve = FALSE preserves every feature", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    dissolve = FALSE, quiet = TRUE
  )
  expect_equal(nrow(db), 4)
})

test_that("trait join keeps unmatched species and warns about them", {
  src <- sf_polygons(fixture_ranges(), species_col = "SCI_NAME")
  expect_message(
    db <- build_taxon_db(src, traits = fixture_traits()),
    "Genus_gamma|NA traits|unmatched",
    all = FALSE
  )
  db <- build_taxon_db(src, traits = fixture_traits(), quiet = TRUE)
  expect_equal(nrow(db), 3)
  expect_true("mass" %in% names(db))
  expect_true(is.na(db$mass[db$species == "Genus_gamma"]))
})

test_that("species column ends up first and geometry last", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    traits = fixture_traits(), quiet = TRUE
  )
  expect_equal(names(db)[1], "species")
  expect_equal(names(db)[ncol(db)], attr(db, "sf_column"))
})

test_that("filter_taxa reports what it removed", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    traits = fixture_traits(), quiet = TRUE
  )
  expect_message(filter_taxa(db, mass < 100), "1/3")
  expect_equal(nrow(filter_taxa(db, mass < 100, quiet = TRUE)), 1)
})

test_that("bad inputs fail with actionable messages", {
  expect_error(build_taxon_db(42), "range source")
  expect_error(build_taxon_db("/nonexistent/folder"), "folder not found")
  expect_error(sf_polygons(data.frame(a = 1)), "sf")
  expect_error(
    sf_polygons(fixture_ranges(), species_col = "nope"),
    "not found"
  )
  expect_error(iucn_shapefile("/nonexistent/data_0.shp") |> resolve_ranges(),
               "not found")
})

test_that("CRS is normalised to EPSG:4326", {
  r <- sf::st_transform(fixture_ranges(), 3857)
  db <- build_taxon_db(sf_polygons(r, species_col = "SCI_NAME"), quiet = TRUE)
  expect_equal(sf::st_crs(db), sf::st_crs(4326))
})

test_that("bundled rodent_traits is well formed and attributed", {
  expect_s3_class(rodent_traits, "tbl_df")
  expect_true(all(c("species", "S_index", "Synanthropic") %in% names(rodent_traits)))
  expect_false("NA" %in% names(rodent_traits))
  expect_equal(anyDuplicated(rodent_traits$species), 0L)
  expect_match(attr(rodent_traits, "source"), "Ecke")
  expect_equal(attr(rodent_traits, "license"), "CC BY 4.0")
})

# --- Traits are optional ------------------------------------------------
# Not every taxon has an Ecke-style trait table. A range source alone must be
# enough to get all the way through the pipeline.

test_that("a database builds with no traits at all", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    quiet = TRUE
  )
  expect_s3_class(db, "richcast_db")
  expect_equal(names(db), c("species", attr(db, "sf_column")))
  expect_equal(nrow(db), 3)
  # Dissolving still happens; it does not depend on traits.
  expect_setequal(db$species, c("Genus_alpha", "Genus_beta", "Genus_gamma"))
})

test_that("filtering a trait-less database explains why the column is missing", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    quiet = TRUE
  )
  expect_error(filter_taxa(db, S_index > 0, quiet = TRUE), "no trait columns")
  expect_error(filter_taxa(db, S_index > 0, quiet = TRUE), "build_taxon_db")
})

test_that("filtering an unknown trait lists the traits that do exist", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    traits = fixture_traits(), quiet = TRUE
  )
  expect_error(filter_taxa(db, nonesuch > 0, quiet = TRUE), "Available trait columns")
})

# --- Reading only the species you want --------------------------------------

test_that("a species subset is read without materialising the rest", {
  # Round-trip through a real geopackage so the OGR query path is exercised,
  # not just mocked.
  skip_if_not(rlang::is_installed("sf"))
  f <- withr::local_tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(
    sf::st_sf(SCI_NAME = c("Genus alpha", "Genus alpha", "Genus beta",
                           "Genus gamma"),
              geometry = sf::st_sfc(square(0, 0), square(5, 0), square(2, 2),
                                    square(8, 8), crs = 4326)),
    f, quiet = TRUE
  ))

  src <- iucn_shapefile(f, species = c("Genus alpha", "Genus beta"))
  db <- build_taxon_db(src, quiet = TRUE)
  expect_setequal(db$species, c("Genus_alpha", "Genus_beta"))
  expect_equal(nrow(db), 2)

  # The two Genus_alpha features are still dissolved into one range.
  bb <- sf::st_bbox(db[db$species == "Genus_alpha", ])
  expect_equal(as.numeric(bb["xmax"]), 6)
})

test_that("underscored names match a spaced source and vice versa", {
  f <- withr::local_tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(
    sf::st_sf(SCI_NAME = c("Genus alpha", "Genus beta"),
              geometry = sf::st_sfc(square(0, 0), square(2, 2), crs = 4326)),
    f, quiet = TRUE
  ))
  db <- build_taxon_db(iucn_shapefile(f, species = "Genus_alpha"), quiet = TRUE)
  expect_equal(db$species, "Genus_alpha")
})

test_that("asking for species that are not there fails clearly", {
  f <- withr::local_tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(
    sf::st_sf(SCI_NAME = "Genus alpha",
              geometry = sf::st_sfc(square(0, 0), crs = 4326)),
    f, quiet = TRUE
  ))
  expect_error(
    build_taxon_db(iucn_shapefile(f, species = "Nothing here"), quiet = TRUE),
    "No features matched"
  )
})

test_that("names containing a quote are escaped, not injected", {
  # Defensive: binomials should never contain quotes, but the query is built
  # by string interpolation, so the escaping needs to hold anyway.
  f <- withr::local_tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(
    sf::st_sf(SCI_NAME = c("Genus alpha", "Genus o'brien"),
              geometry = sf::st_sfc(square(0, 0), square(2, 2), crs = 4326)),
    f, quiet = TRUE
  ))
  db <- build_taxon_db(iucn_shapefile(f, species = "Genus o'brien"),
                       quiet = TRUE)
  expect_equal(nrow(db), 1)
})

test_that("printing never reports a negative trait count", {
  db <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    quiet = TRUE
  )
  expect_message(print(db), "0 trait columns")
  # st_drop_geometry() keeps the class but drops a column; the count must not
  # go negative.
  expect_message(print(sf::st_drop_geometry(db)), "0 trait columns")

  with_traits <- build_taxon_db(
    sf_polygons(fixture_ranges(), species_col = "SCI_NAME"),
    traits = fixture_traits(), quiet = TRUE
  )
  expect_message(print(with_traits), "2 trait columns")
})

test_that("a vector where a column name belongs is diagnosed, not left to %in%", {
  # R partial-matches `species=` onto `species_col=` on versions lacking a
  # `species` argument, which used to surface as "the condition has length > 1".
  expect_error(
    sf_polygons(fixture_ranges(), species_col = c("Genus alpha", "Genus beta")),
    "single column name"
  )
  expect_error(
    sf_polygons(fixture_ranges(), species_col = c("Genus alpha", "Genus beta")),
    "richcast is up to date"
  )
  expect_error(sf_polygons(fixture_ranges(), species_col = NA), "single column name")
  expect_error(sf_polygons(fixture_ranges(), species_col = 1), "single column name")
})

# --- Reading a folder of IUCN downloads ---------------------------------------

write_iucn <- function(dir, name, df) {
  d <- file.path(dir, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(sf::st_write(df, file.path(d, "data_0.shp"), quiet = TRUE))
}

test_that("iucn_folder reads every shapefile under a folder as one source", {
  dir <- withr::local_tempdir()
  write_iucn(dir, "bovid_IUCN", sf::st_sf(
    SCI_NAME = c("Genus alpha", "Genus alpha"), PRESENCE = 1, ORIGIN = 1,
    EXTRA = "only here",
    geometry = sf::st_sfc(square(0, 0), square(5, 0), crs = 4326)
  ))
  write_iucn(dir, "cervid_IUCN", sf::st_sf(
    SCI_NAME = "Genus beta", PRESENCE = 1, ORIGIN = 1,
    geometry = sf::st_sfc(square(2, 2), crs = 4326)
  ))

  db <- build_taxon_db(iucn_folder(dir), quiet = TRUE)
  expect_setequal(db$species, c("Genus_alpha", "Genus_beta"))
  bb <- sf::st_bbox(db[db$species == "Genus_alpha", ])
  expect_equal(as.numeric(bb["xmax"]), 6)

  # A bare folder path means the same thing.
  expect_setequal(build_taxon_db(dir, quiet = TRUE)$species, db$species)
})

test_that("extinct and introduced polygons are left out by default", {
  dir <- withr::local_tempdir()
  write_iucn(dir, "g", sf::st_sf(
    SCI_NAME = c("Genus alpha", "Genus alpha", "Genus alpha"),
    PRESENCE = c(1, 5, 1), ORIGIN = c(1, 1, 3),
    geometry = sf::st_sfc(square(0, 0), square(5, 0), square(0, 5), crs = 4326)
  ))
  db <- build_taxon_db(iucn_folder(dir), quiet = TRUE)
  bb <- sf::st_bbox(db)
  expect_equal(as.numeric(bb[c("xmax", "ymax")]), c(1, 1))

  all_codes <- build_taxon_db(iucn_folder(dir, presence = NULL, origin = NULL),
                              quiet = TRUE)
  expect_equal(as.numeric(sf::st_bbox(all_codes)[c("xmax", "ymax")]), c(6, 6))
})

test_that("an empty folder fails clearly", {
  dir <- withr::local_tempdir()
  expect_error(build_taxon_db(iucn_folder(dir), quiet = TRUE), "No .*files")
})

test_that("point geometry is refused", {
  pts <- sf::st_sf(species = "Genus_alpha",
                   geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326))
  expect_error(build_taxon_db(pts, quiet = TRUE), "must be polygons")
})
