# ==============================================================================
# The evidence behind min_cells
#
# A documentation-only topic. The default cutoff governs which species enter an
# assemblage total, so the measurements behind it belong in the manual rather
# than in a commit message.
# ==============================================================================

#' Resolvability: why `min_cells` defaults to what it does
#'
#' @description
#' [run_hindcast_series()] holds species whose ranges cover fewer than
#' `min_cells` grid cells out of the richness surfaces. This topic records the
#' measurements the default rests on, so the number can be argued with.
#'
#' @section What goes wrong:
#'
#' `fit_sdm()` draws presences **with replacement** from the cells the range
#' covers. A range covering `k` cells therefore supplies exactly `k` distinct
#' climate vectors, whatever `nsample` is asked for. Across the eight Levantine
#' ungulates at 0.5 degrees:
#'
#' | species | range cells | 1000 samples give | AUC |
#' |---|---|---|---|
#' | *Dama mesopotamica* | 11 | each cell about 91 times | 0.987 |
#' | *Gazella gazella* | 19 | each cell about 53 times | 0.905 |
#' | *Capra ibex* | 79 | each cell about 13 times | 0.935 |
#' | *Capra aegagrus* | 747 | 544 distinct cells | 0.919 |
#' | *Sus scrofa* | 12549 | 966 distinct cells | 0.908 |
#'
#' Note the AUCs. *Dama* carries the highest in the assemblage on 11 distinct
#' points: separating 11 observations from 1000 background points is easy, and
#' the statistic says nothing about whether the resulting surface means
#' anything. AUC runs the wrong way here and must not be used as this screen.
#'
#' @section What it does to a trajectory:
#'
#' Ten seeds per species, six slices spanning 2-60 ka BP, default settings.
#' `traj_cv` is the coefficient of variation of the projected extent across
#' seeds, taken per slice and then medianed:
#'
#' | species | range cells | present-day CV | traj_cv |
#' |---|---|---|---|
#' | *Dama mesopotamica* | 11 | 14.9% | 58.6% |
#' | *Gazella gazella* | 19 | 2.9% | 95.8% |
#' | *Capra ibex* | 79 | 10.0% | 24.9% |
#' | *Capra aegagrus* | 747 | 2.9% | 6.8% |
#' | *Alcelaphus buselaphus* | 1243 | 1.7% | 9.1% |
#' | *Cervus elaphus* | 2479 | 3.0% | 7.2% |
#' | *Capreolus capreolus* | 3633 | 3.0% | 3.8% |
#' | *Sus scrofa* | 12549 | 3.2% | 5.3% |
#'
#' The present-day column is the trap. *Gazella*'s present-day estimate is
#' among the steadiest in the table (2.9%) while its hindcast trajectory is the
#' least stable (95.8%). Screening on how firm the present-day range looks
#' would have passed it.
#'
#' @section Why not replicate spread:
#'
#' Because it inverts. Replicates redraw the sample, so they measure sampling
#' variability -- and a range small enough to be sampled exhaustively has none.
#'
#' Clipping one donor range (*Cervus elaphus*) to compact blocks of a fixed
#' cell count, holding climate, grid, model and slices constant, five seeds per
#' block: an 8-cell range returned `2, 2, 2, 2, 2` cells at one slice and
#' `0, 0, 0, 0, 0` at two others -- perfect agreement across every seed, at a
#' range size where the answer is worthless. Whole slices collapsed to zero on
#' every seed for at least one block at every size up to 256 cells, and at none
#' from 512 up.
#'
#' So a narrow replicate interval is not evidence of a firm answer; below a few
#' dozen cells it is evidence that there was nothing left to resample. This is
#' why the screen is on `range_cells` and not on `cells_sd`.
#'
#' @section Why not study extent or suitable fraction:
#'
#' Both move with the fitting buffer, which is a modelling choice rather than a
#' property of the species. Forcing an exact background ring of *B* degrees on
#' *Gazella gazella*:
#'
#' | buffer | study extent | present cells | suitable fraction | traj CV |
#' |---|---|---|---|---|
#' | 1 deg | 81 | 16 | 19.8% | 71% |
#' | 2 deg | 172 | 16 | 9.3% | 125% |
#' | 4 deg | 433 | 18 | 4.2% | 64% |
#' | 8 deg | 1154 | 41 | 3.6% | 60% |
#'
#' The extent moves 14-fold and the suitable fraction 5-fold on one unchanged
#' species, while the trajectory stays equally unstable throughout. `range_cells`
#' does not move, because the buffer does not change the range.
#'
#' The same table settles a second question. Widening the background does shift
#' the level of the estimate (16 to 41 cells) and the `p10` cutoff with it
#' (0.573 to 0.217), so `p10` is not as background-insensitive as its
#' presence-referenced definition suggests: the *rule* reads off the presences,
#' but the *model* generating those predictions is fitted against the
#' background. What widening the buffer does not do is make a small range
#' resolvable.
#'
#' @section Where the default comes from:
#'
#' The eight ungulates alone leave a gap: 79 cells unstable, 747 stable,
#' nothing measured in between. A cutoff dropped into that gap would be picked,
#' not measured. So 23 further species -- rodents whose ranges fall inside the
#' same Beyer extent, modelled on the same 0.5-degree grid, five seeds each --
#' were measured to fill it, giving 31 species spanning 11 to 12549 cells.
#'
#' A trajectory is only as good as its worst slice, so the statistic below is
#' the largest per-slice coefficient of variation across seeds, not the median:
#'
#' | range cells | n | worst-slice CV: min | median | max |
#' |---|---|---|---|---|
#' | under 50 | 4 | 38.0% | 123.5% | 194.4% |
#' | 50-100 | 4 | 25.2% | 48.3% | 103.7% |
#' | 100-200 | 4 | 11.6% | 19.3% | 24.0% |
#' | 200-400 | 5 | 5.2% | 70.3% | 119.8% |
#' | 400-800 | 6 | 7.9% | 29.1% | 223.6% |
#' | over 800 | 8 | 4.4% | 17.5% | 36.5% |
#'
#' **All 8 species below 100 cells have a worst-slice CV above 25%.** That is
#' what fixes the default: 100 is the largest cutoff at which every excluded
#' species is demonstrably unstable. Raise it to 200 and the rule stops holding
#' -- only 8 of the 12 species excluded would deserve it, and *Marmota
#' marmota* (168 cells, worst-slice CV 11.6%) would be thrown out for nothing.
#'
#' The 32 Tian Shan rodents clear it by a factor of 31: the narrowest, *Marmota
#' baibacina*, covers 3161 cells at that vignette's resolution. This is why the
#' screen once looked unnecessary. The evidence it was withdrawn on came
#' entirely from wide-ranging Eurasian rodents, none of which could exhibit the
#' failure -- not because the failure was not real, but because that assemblage
#' had no species small enough to show it.
#'
#' @section What it does not fix:
#'
#' Instability does not stop at the cutoff, and range size stops explaining it
#' there. Above 100 cells, 8 of 23 species still have a worst-slice CV above
#' 25%, and the relationship is no longer monotone -- the 200-400 bin is worse
#' than the 100-200 bin. *Praomys rostratus* covers 236 cells and still varies
#' 69% across seeds.
#'
#' That residue is concentrated in the deep slices: across species above the
#' cutoff, the median across-seed CV runs 7.2% at 2 ka BP against 19.6% at 22
#' ka and 22.8% at 60 ka. The pattern is consistent with projection beyond the
#' climate the model was fitted on, which is a different problem with a
#' different remedy, and one this screen does not address. Clearing `min_cells`
#' means the grid can resolve the species. It does not mean the trajectory is
#' trustworthy; check `cells_sd` from `replicates > 1` for that, bearing in
#' mind the caveat above about ranges too small to resample.
#'
#' @section Cells, not area:
#'
#' `range_cells` counts grid cells because the grid is what limits the
#' inference. The same species is better resolved on a finer reconstruction,
#' and worse on a coarser one, with no change to its biology. A cutoff in cells
#' is therefore only meaningful alongside the resolution it was applied at,
#' which is why `$resolvability` records the `min_cells` used and
#' `$models` carries `range_cells` per species.
#'
#' @section Occurrence data measures something slightly different:
#'
#' With a polygon range, `range_cells` is a property of the species and the
#' grid. With occurrences from [gbif_occurrences()] it is a property of the
#' species, the grid **and the survey effort**: it counts the distinct cells
#' that hold a record, so a well-mapped species with few records looks small.
#'
#' Drawing occurrences inside three Levantine ungulate ranges and refitting on
#' those instead of the polygon:
#'
#' | species | polygon | 30 records | 300 | 3000 |
#' |---|---|---|---|---|
#' | *Dama mesopotamica* | 11 | 9 | 11 | 11 |
#' | *Gazella gazella* | 19 | 9 | 18 | 20 |
#' | *Capra aegagrus* | 747 | 30 | 232 | 626 |
#'
#' The narrow-ranged taxa converge on their polygon count within a few hundred
#' records, because there are only so many cells to find. *Capra aegagrus* does
#' not: at 30 records it reports 30 cells and would be screened out, though its
#' range covers 747. The screen is still telling the truth -- a model given 30
#' points really does see 30 climate vectors -- but the remedy is more records,
#' not a coarser conclusion about the species.
#'
#' So with occurrence data, read an exclusion as "not enough distinct records
#' to place this species on this grid", and check `range_cells` against the
#' record count in `$models` before concluding anything about range size.
#'
#' @section What the screen is not:
#'
#' It is not a statement that the species was absent, rare, or unimportant. It
#' says the reconstruction is too coarse to place it, which is a property of
#' the data, not the animal. Excluded species keep every row they had:
#' `$species`, `$models`, `$ranges` and `$fits` are unchanged, and only the
#' richness surfaces differ. Read `$resolvability` alongside any richness
#' figure so the assemblage is never quietly smaller than it looks.
#'
#' @seealso [run_hindcast_series()], [fit_sdm()]
#' @name richcast-resolvability
#' @keywords internal
NULL
