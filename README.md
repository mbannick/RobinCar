[![R-CMD-check](https://github.com/mbannick/RobinCar/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/mbannick/RobinCar/actions/workflows/R-CMD-check.yaml) ![CRAN](https://www.r-pkg.org/badges/version/RobinCar) ![downloads](https://cranlogs.r-pkg.org/badges/grand-total/RobinCar) ![downloads](https://cranlogs.r-pkg.org/badges/RobinCar) [![DOI](https://zenodo.org/badge/506080289.svg)](https://zenodo.org/badge/latestdoi/506080289)


# RobinCar: ROBust estimation and INference for Covariate Adjustment in Randomized clinical trials

RobinCar is a package that allows for robust estimation and inference for treatment effects in randomized clinical trials when covariates are used at the design and/or analysis stages of the trial. Supported covariate-adaptive randomization schemes at the design phase are simple randomization, stratified permuted block randomization, biased coin randomization, and Pocock and Simon's minimization. Statistical methods at the analysis stage are model-assisted and assumption-lean, in accordance with [FDA guidance on covariate adjustment](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/adjusting-covariates-randomized-clinical-trials-drugs-and-biological-products). Publications describing the methods are listed [here](#publications).

See also [RobinCar2](https://openpharma.github.io/RobinCar2/main/), which is a lite version of RobinCar and is supported by the [ASA Biopharmaceutical Section Covariate Adjustment Scientific Working Group Software Subteam](https://carswg.github.io/subteam_software.html).

## Authors
[Ting Ye](https://sites.google.com/view/tingye), Yanyao Yi, [Marlena Bannick](https://marlenabannick.com/) (maintainer), Yuhan Qian, and Faith Bian

## Documentation

To view documentation about the functions, see the RobinCar website here: https://marlenabannick.com/RobinCar/. You will also find vignettes about how to use the functions.

## Installation

RobinCar is now available on CRAN!

### 1. Install via CRAN
![CRAN](https://www.r-pkg.org/badges/version-last-release/RobinCar)

```{}
install.packages("RobinCar")
```

### 2. Install with `devtools`

To get the most recent version in development, you can install the package with `devtools`:
```{bash}
devtools::install_github("mbannick/RobinCar")
```

### 3. Clone repository

Or to download the package, you may clone the repository:
```{bash}
git clone https://github.com/mbannick/RobinCar.git
```

## <a id="publications"></a>Publications

Here are publications and preprints that explain the methods in RobinCar:

* [Ting Ye, Jun Shao, Yanyao Yi (2023)](https://doi.org/10.1093/biomet/asad045)
* [Marlena Bannick, Jun Shao, Jingyi Liu, Yu Du, Yanyao Yi, Ting Ye (2025)](https://doi.org/10.1093/biomet/asaf029)
* [Ting Ye, Marlena Bannick, Yanyao Yi, Jun Shao (2023)](https://doi.org/10.1080/24754269.2023.2205802)
* [Ting Ye, Jun Shao, Yanyao Yi, Qinyuan Zhao (2023)](https://doi.org/10.1080/01621459.2022.2049278)
* [Xiaoyu Qiu, Yuhan Qian, Jaehwan Yi, Jinqiu Wang, Yu Du, Yanyao Yi, Ting Ye (2025)](https://arxiv.org/pdf/2408.12541)
