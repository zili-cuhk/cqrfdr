# cqrfdr 0.3.0

* `hdcqr(lambda = "zheng")`: the penalty sequence of Zheng, Peng and He (2018, Section 5.1), with a pivotal level at the first grid point and an increment chosen by cross-validation at every grid point; this is the L-HDCQR comparator of the paper. It reproduces the paper's implementation exactly on the same seed and folds.
* README and vignette describe MS-MS, the comparators and the diagnostics used in the paper (share of empty splits, agreement between runs).
* `inst/CITATION` added.

# cqrfdr 0.2.0

* New `msms()`: the multiple-splitting version of SS-MS (SS-MS on `B` random splits, aggregated by inclusion rates as MDS aggregates DS in Dai et al., 2023). The result records the share of empty splits and the split-wise selection sizes.
* New `ds()`: the single data-splitting procedure of Dai et al. (2023) (CV support without a cap, 10-fold CV, cutoff U), as a comparator of SS-MS.
* `mds()` now follows the design of Dai et al. (2023) by default: every split is a complete DS fit with its own 10-fold CV (`tune = "cv_each"`), no screening cap (`screen_cap = "none"`) and their tie rule (`ties = "dai"`). The previous behaviour is `mds(..., tune = "cv_once", nfolds = 5, screen_cap = "d2", ties = "random")`.
* `mds()` and `msms()` share one internal routine and report `split_sizes` and `share_empty`.
* With a fixed penalty multiplier (`lambda = "pivotal"`), the `"path"` and `"size"` screens now move down the multiplier grid as documented; before, they kept the support at the fixed multiplier.

# cqrfdr 0.1.0

* First version. `ssms()` (path-ordered screen, U+ threshold, several FDR levels per fit), `mds()` (CV multiplier chosen once, random tie-breaking), `hdcqr()` with `coef()`/`predict()`, the solver `hdcqr_path()`, `pivotal_lambda()`, `mirror_stat()`, `mirror_threshold()`, `sim_cqr()` and `fdp_power()`.
* Results agree exactly with the simulation code used for the paper on the same seed.
