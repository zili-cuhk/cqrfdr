# cqrfdr

FDR-controlled variable selection for high-dimensional censored quantile regression.

R package accompanying Yu, J. and Liu, Z., *FDR-Controlled Variable Selection for High-Dimensional Censored Quantile Regression* (manuscript). Everything is written in base R, so no compiler is needed.

| Function | Method |
|---|---|
| `ssms()` | **SS-MS** (proposed): one random split; penalized sequential Peng–Huang fit and path-ordered screen on the first half, unpenalized refit on the second half, mirror statistics with the cutoff U₊ |
| `msms()` | **MS-MS** (proposed): SS-MS on B random splits, aggregated by inclusion rates |
| `ds()` | DS of Dai, Lin, Xing and Liu (2023) |
| `mds()` | MDS of Dai, Lin, Xing and Liu (2023) |
| `hdcqr()` | L-HDCQR of Zheng, Peng and He (2018); `lambda = "zheng"` gives their penalty sequence |
| `sim_cqr()`, `fdp_power()` | simulation designs of the paper and evaluation |

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("zili-cuhk/cqrfdr")   # the vignette PDF is included

# or from the source tarball
install.packages("cqrfdr_0.3.0.tar.gz", repos = NULL, type = "source")
```

## Example

```r
library(cqrfdr)
set.seed(2026)
d <- sim_cqr(n = 300, p = 200, s = 10, beta = 0.8, cens = 0.3)

fit <- ssms(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2))
fit                       # selections and thresholds at both levels
selected(fit, q = 0.1)    # names of the selected covariates
fdp_power(selected(fit, 0.1), paste0("X", d$active))

fm <- msms(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2), B = 50)
selected(fm, q = 0.1)
fm$share_empty            # share of splits whose own SS-MS selection was empty
```

With real data the response can be a `survival::Surv` object, for example `msms(X, Surv(time, status), tau = 0.5, q = 0.1)`; the times are then log-transformed.

`vignette("cqrfdr")` explains the options, the comparators and the checks to run on real data before looking at a selection (screening cap and plasmode simulation).

## Defaults (the settings of the paper)

- Grid of 40 quantile levels on [0.1, τ].
- Penalty: pivotal level times a multiplier chosen by 5-fold CV of the held-out Peng–Huang objective.
- Path-ordered screen, capped at ⌊d₂/log d₂⌋ candidates (d₂ = events in the second half).
- Cutoff U₊ (`offset = 1`).
- `msms()`: B = 50 splits, multiplier chosen once, ties among inclusion rates broken at random.
- `ds()`, `mds()`: as in Dai et al. (2023): 10-fold CV within each split, candidate set = CV support without a cap, cutoff without the offset; `mds()` uses their aggregation rule.
- `hdcqr(lambda = "zheng")`: pivotal penalty at the first grid point, increment chosen by 5-fold CV at every grid point.

## Reproducing the paper

The simulation study, the real-data analyses and the SCQR comparators (He, Pan, Tan and Zhou, 2022; the authors' C++ code) are in the replication materials of the paper, https://github.com/zili-cuhk/cqrfdr-replication. They implement the same algorithms as this package.

## References

- Dai, C., Lin, B., Xing, X. and Liu, J. S. (2023). False discovery rate control via data splitting. *Journal of the American Statistical Association*, 118, 2503–2520.
- He, X., Pan, X., Tan, K. M. and Zhou, W.-X. (2022). Scalable estimation and inference for censored quantile regression process. *The Annals of Statistics*, 50, 2899–2924.
- Peng, L. and Huang, Y. (2008). Survival analysis with quantile regression models. *Journal of the American Statistical Association*, 103, 637–649.
- Zheng, Q., Peng, L. and He, X. (2018). High dimensional censored quantile regression. *The Annals of Statistics*, 46, 308–343.

## License

GPL (≥ 3).
