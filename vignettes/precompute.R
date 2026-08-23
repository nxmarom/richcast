# ==============================================================================
# vignettes/precompute.R
#
# The Tian Shan vignette runs against a global IUCN Red List download and ~170
# MB of prepared climate slices, neither of which can be shipped or downloaded
# at build time. So it is precomputed: tianshan.Rmd.orig holds the real code,
# is knitted once on a machine that has the data, and the resulting
# tianshan.Rmd -- with outputs and figures baked in -- is what ships.
#
# Re-run this from the package root after changing tianshan.Rmd.orig:
#
#   Rscript vignettes/precompute.R
#
# Expect roughly 10 minutes. Requires:
#   ../Inputs/geo_data/iucn_rodentia/data_0.shp
#   ../Inputs/Climate/eurasia_slices/
# ==============================================================================

stopifnot(file.exists("DESCRIPTION"))

# Wrapped in a function so on.exit() has a frame to attach to. At top level in
# a sourced script it fires immediately, reverting the directory before knit
# runs.
precompute <- function() {
  old <- setwd("vignettes")
  on.exit(setwd(old), add = TRUE)

  # Figures go straight into vignettes/ with a per-vignette prefix, NOT into a
  # subdirectory. R CMD check treats a vignettes/figure/ directory as leftover
  # knitr debris and NOTEs about it, which for a precomputed vignette is a
  # false positive -- those files are the deliverable, not scratch.
  # error = FALSE so a failing chunk aborts the precompute. knitr's default is
  # to capture the error, print it into the output and carry on -- which for a
  # precomputed vignette means shipping a broken one, with stale figures left
  # over from the previous run still sitting on disk looking plausible.
  knitr::opts_chunk$set(fig.path = "tianshan-", error = FALSE)
  knitr::knit("tianshan.Rmd.orig", output = "tianshan.Rmd")
}

precompute()

cat("\nPrecomputed vignettes/tianshan.Rmd\n")
cat("Figures written to vignettes/figure/\n")
