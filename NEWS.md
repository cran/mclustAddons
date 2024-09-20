# mclustAddons 0.9 (2024-09) NOT ON CRAN

-  Added `GMMlogreturn()` to model log-returns financial data via GMMs.
-  Added `VaR()` and `ES()` to compute risk measures from GMMs.
- `densityMclustBounded()` allows to have lambda values fixed for a 
  subset of variables while to estimate lambda parameters for the 
  remaining variables.
- `densityMclustBounded()` added noise component.
- converted all documentation to use `roxygen2` package.
- clean-up package dependencies.

# mclustAddons 0.8 (2024-02) 

- Starting with this version mclust >= 6.1 is needed because of using
  efficient `softmax()` and `logsumexp()` functions coded in Fortran.
  This replaces previously included functions of the same name (now
  removed).
- Added `prior` and `nstart` arguments to `densityMclustBounded()` 
  function.
- Added `rangepowerBackTransform()`.
- Exported both `rangepowerTransform()` and `rangepowerBackTransform()`.
- Bug fixes.

# mclustAddons 0.7.2 (2023-01)

- Bug fixes.

# mclustAddons 0.7.1 (2022-11)

- Small bug fix in help page.

# mclustAddons 0.7 (2022-11)

- Add code for computing entropy via Gaussian mixtures. See 
  `help(EntropyGMM)`.
- Add efficient Rcpp-based functions to compute log-sum-exp and softmax.

# mclustAddons 0.6 (2021-12)

- Bug fixes.
- Allowed to have both `lbound = NULL` and `ubound = NULL`, in which 
  case only power transformation is performed.

# mclustAddons 0.5 (2021-09)

- Release submitted to CRAN.

# mclustAddons 0.4 (2020-07 NOT ON CRAN)

- Add logarithm arg to `predict.densityMclustBounded()` function.

# mclustAddons 0.3 (2020-05 NOT ON CRAN)

- Release for R version 4.0 
- Inclusion of `MclustMEM()` and accompanying functions.

# mclustAddons 0.2 (2017-10 NOT ON CRAN)

- Inclusion of `densityMclustBounded()` and accompanying functions.

# mclustAddons 0.1 (2014-10 NOT ON CRAN)

- Initial developing release.
