# simproxyage: Simulate Age Uncertainty in Climate Proxy Records.

------------------------------

## Introduction

**simproxyage** implements an stochastic approach to simulate the age
uncertainty in layer-counted climate proxy records by perturbing a reference
chronology by randomly omitting or double-counting the reference layers.

The implementation is an adapted and extended version of the MATLAB code
developed by Comboul et al. (2014) who applied it on annually banded coral
archives. The current R package has been used by Münch and Laepple (2018) to
model the time uncertainty of annually dated stable isotope firn-core records
based on seasonal layer counting constrained by volcanic tie points.

The original MATLAB code is available at the [Common Climate Google Code
repository](https://code.google.com/p/common-climate/). The R package version
hosted here has been implemented by Dr. Thomas Münch. Please contact me
<<thomas.muench@awi.de>> at the Alfred Wegener Institute, Helmholtz Centre for
Polar and Marine Research, Germany, for further information.
 
## Installation

**simproxyage** can be installed directly from GitHub:

```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("EarthSystemDiagnostics/simproxyage")
```

## Literature cited

Comboul, M. , Emile-Geay, J., Evans, M. N., Mirnateghi, N., Cobb, K. M. and
Thompson, D. M.: A probabilistic model of chronological errors in layer-counted
climate proxies: applications to annually banded coral archives, Clim. Past,
10(2), 825-841, doi:
[10.5194/cp-10-825-2014](https://doi.org/10.5194/cp-10-825-2014), 2014.

Münch, T. and Laepple, T.: What climate signal is contained in decadal to
centennial scale isotope variations from Antarctic ice cores?, Clim. Past,
accepted for publication, doi:
[10.5194/cp-2018-112](https://doi.org/10.5194/cp-2018-112), 2018.


