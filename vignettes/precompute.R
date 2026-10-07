# Re-render the precomputed vignette(s).
#
# middle-east.Rmd.orig needs IUCN range polygons and prepared Beyer2020
# climate slices, neither of which can ship with the package. It is knitted
# here, in a local analysis folder holding them, and the knitted markdown and
# figures are copied back into vignettes/ as the vignette source R CMD build
# sees. Run from the package root, with the analysis folder laid out as:
#
#   analysis/middle_east_120ka/
#     UngulatePolygons/   IUCN downloads (or a symlink to them)
#     beyer/              climate slices from prepare_climate()
#     runs/               saved model runs (created on first render)

analysis <- "analysis/middle_east_120ka"
if (!file.exists(file.path(analysis, "UngulatePolygons"))) {
  file.symlink(normalizePath("UngulatePolygons"),
               file.path(analysis, "UngulatePolygons"))
}

# The render must use this checkout of richcast, not an older installed one.
lib <- file.path(analysis, "lib")
dir.create(lib, showWarnings = FALSE)
install.packages(".", lib = lib, repos = NULL, type = "source", quiet = TRUE)
.libPaths(c(normalizePath(lib), .libPaths()))

file.copy("vignettes/middle-east.Rmd.orig", analysis, overwrite = TRUE)
old <- setwd(analysis)
knitr::knit("middle-east.Rmd.orig", "middle-east.Rmd")
setwd(old)

file.copy(file.path(analysis, "middle-east.Rmd"), "vignettes", overwrite = TRUE)
figs <- list.files(analysis, pattern = "^middle-east-.*\\.png$", full.names = TRUE)
file.copy(figs, "vignettes", overwrite = TRUE)

# --- richcast-quickstart.Rmd.orig ------------------------------------------------------
# Knitted in its own folder, laid out as the vignette describes:
#
#   analysis/quickstart/
#     geo_data/   IUCN downloads (here a symlink to UngulatePolygons)
#     beyer/      climate slices (here a symlink to the Middle East folder's)
#     lib/        this checkout of richcast
qs <- "analysis/quickstart"
dir.create(qs, showWarnings = FALSE)
if (!file.exists(file.path(qs, "geo_data"))) {
  file.symlink(normalizePath("UngulatePolygons"), file.path(qs, "geo_data"))
}
if (!file.exists(file.path(qs, "beyer"))) {
  file.symlink(normalizePath(file.path(analysis, "beyer")), file.path(qs, "beyer"))
}
qlib <- file.path(qs, "lib")
dir.create(qlib, showWarnings = FALSE)
install.packages(".", lib = qlib, repos = NULL, type = "source", quiet = TRUE)
.libPaths(c(normalizePath(qlib), .libPaths()))
file.copy("vignettes/richcast-quickstart.Rmd.orig", qs, overwrite = TRUE)
old <- setwd(qs)
knitr::knit("richcast-quickstart.Rmd.orig", "richcast-quickstart.Rmd")
setwd(old)
file.copy(file.path(qs, "richcast-quickstart.Rmd"), "vignettes", overwrite = TRUE)
file.copy(list.files(qs, pattern = "^richcast-quickstart-.*\\.png$", full.names = TRUE), "vignettes", overwrite = TRUE)
