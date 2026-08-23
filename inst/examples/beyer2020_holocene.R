# ==============================================================================
# Beyer2020: Holocene rodent richness in the Tian Shan
#
# A worked configuration for a deep-time reconstruction. Beyer et al. (2020)
# covers 0 to -120,000 BP at 0.5 degrees, globally, so it reaches far further
# back than CHELSA-TraCE21k and needs no aggregation.
#
# Three things differ from a WorldClim/CHELSA run, and all three are
# configuration rather than code:
#
#   1. Time is quoted BP but richcast works in CE, and "present" in the BP
#      convention means 1950. Use bp_to_ce() rather than writing the offset
#      out by hand.
#   2. Beyer2020 publishes its own 0 BP slice, so fit and projection can stay
#      within one product instead of crossing from WorldClim.
#   3. It is already at 0.5 degrees, so both aggregation factors are 1.
# ==============================================================================

library(richcast)
library(pastclim)

pastclim::set_data_path("~/pastclim_data/")

VARS <- c("bio01", "bio04", "bio05", "bio06",
          "bio12", "bio15", "bio16", "bio17")

# Beyer2020 steps every 1000 yr back to 22 ka BP, then every 2000 yr. Staying
# on a 2000-yr grid keeps the spacing regular across the whole range.
TIMES_BP <- seq(-2000, -20000, by = -2000)
TIMES_CE <- bp_to_ce(TIMES_BP)          # 0 BP == 1950 CE

# --- 1. Stage the slices (once; resumable) ----------------------------------
clim <- prepare_climate(
  path            = "climate/beyer",
  vars            = VARS,
  times           = TIMES_CE,
  extent          = c(60, 95, 32, 52),   # Tian Shan plus buffer for fitting
  dataset_present = "Beyer2020",         # fit on Beyer's own present slice
  dataset_past    = "Beyer2020",
  present_time    = 1950,
  agg_present     = 1,                   # already 0.5 deg
  agg_past        = 1,
  present_dir     = "beyer_1950"
)

check_climate_grids(clim, VARS)

# --- 2. Taxa ----------------------------------------------------------------
db <- build_taxon_db(
  ranges = iucn_shapefile("~/iucn_rodentia/data_0.shp"),
  traits = rodent_traits
)
db <- filter_taxa(db, !is.na(S_index))

# --- 3. Run -----------------------------------------------------------------
# Richness resolution should not out-run the climate: at 0.5 degree input,
# a 0.1 degree richness grid invents detail the model never had.
res <- run_hindcast_series(
  db, clim,
  times      = TIMES_CE,
  focus      = focus_box(c(68, 87, 39, 46), label = "Tian Shan"),
  predictors = VARS,
  resolution = 0.25,
  on_error   = "warn"
)

# --- 4. Report in BP, which is how deep time is read ------------------------
res$richness |>
  dplyr::mutate(bp = ifelse(period == "present", NA, ce_to_bp(time))) |>
  dplyr::select(period, time_ce = time, bp, mean_richness, max_richness)

# Maps per slice
richness_grid(res)
richness_surface(res, bp_to_ce(-6000))

# --- Notes ------------------------------------------------------------------
# Chronological smoothing: gaussian_window(step = 2000) matches the coarser
# half of Beyer2020. step = 1000 works only back to 22 ka BP; beyond that it
# silently lands between slices, and richcast drops the missing ones and
# renormalises rather than erroring.
#
# Global runs: Beyer2020 covers -180/180, -60/90, so extent = c(-180, 180,
# -60, 90) and focus_global() are viable here in a way they are not with
# Eurasia-clipped CHELSA slices. Expect every trait-matched species to be
# modelled rather than the ~32 that reach the Tian Shan.
#
# Missing variables: Beyer2020 has no bio02 or bio03. The eight above are all
# present.
