# richcast extdata

## zooarch_taxa.csv

The Pleistocene zooarchaeological species lists behind `zooarch_taxa()`;
sources are given per row and in `?zooarch_taxa`.

## levant_site_series.csv

Dated archaeological site-phases per richcast slice (Beyer2020, 0-120 ka BP)
in the Levantine corridor, 34.0-37.5 E by 29.5-37.0 N, used in
`vignette("middle-east")`. One row per slice; no individual dates are
included. Built outside the package from:

* ROAD, the ROCEEH Out of Africa Database (Kandel et al. 2023,
  doi:10.1371/journal.pone.0289513), retrieved with the roadDB R package
  (doi:10.32614/CRAN.package.roadDB). CC BY-SA 4.0.
* NERD, Near East Radiocarbon Dates (Palmisano et al. 2022,
  doi:10.5334/joad.90; https://github.com/apalmisano82/NERD). CC BY 4.0.
* p3k14c (Bird et al. 2022, doi:10.1038/s41597-022-01118-7). Data CC0,
  attribution requested.
* Calibrated with IntCal20 (Reimer et al. 2020, doi:10.1017/RDC.2020.41);
  taphonomic correction after Surovell et al. (2009,
  doi:10.1016/j.jas.2009.03.029).

Columns: `ka_bp`; `levant_phases`, expected site-phases in the slice's
window; `levant_per_ka`, the same per ka; `levant_per_ka_taph`, corrected
for taphonomic loss; `levant_per_ka_q05`, `levant_per_ka_q95`, the 90%
Monte Carlo envelope of `levant_per_ka`. The method is described in the
vignette.

**Licence: CC BY-SA 4.0** (https://creativecommons.org/licenses/by-sa/4.0/),
because ROAD content is published under it. This file is not covered by the
package's MIT licence. Cite the sources above when you use it.
