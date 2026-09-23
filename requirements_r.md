# R package requirements

Tested with R 4.3 on Linux (RHEL-like cluster) and macOS.

All packages are on CRAN. Install everything with:

```r
install.packages(c(
  # Data handling and plotting
  "dplyr", "readr", "tidyr", "purrr", "stringr", "tibble", "forcats", "ggplot2",
  # Cohort + matching
  "MatchIt", "mice", "cobalt",
  # Statistical models
  "car", "emmeans", "sandwich", "lmtest",
  # Figure composition and export
  "patchwork", "svglite", "circlize"
))
```

## Package roles

| Package | Used in | Purpose |
|---------|--------|---------|
| `dplyr`, `readr`, `tidyr`, `purrr`, `stringr`, `tibble`, `forcats` | all R scripts | data manipulation |
| `MatchIt` | `01_cohort_assembly.R` | 1:2 nearest-neighbour propensity score matching |
| `mice` | `01_cohort_assembly.R` | chained-equations imputation of clinical covariates |
| `cobalt` | `01_cohort_assembly.R` | covariate balance after matching |
| `car` | `07a`, `07b`, `07f` | Type III sums of squares |
| `emmeans` | `07a`, `07b`, `07f` | estimated marginal means for adjusted Cohen's d |
| `sandwich`, `lmtest` | `07c`, `07h`, `08a` | HC3 heteroscedasticity-consistent standard errors |
| `ggplot2` | `01b`, `07a`, `07c`, `07e`, `08a` | figures |
| `patchwork` | `08a` | multi-panel composition of Figure 3 |
| `svglite` | figure scripts | SVG export |
| `circlize` | `08b` | chord diagram of Figure 1 |

## R-on-cluster install tips

If installing on a shared cluster, you may need to set a per-user library:

```r
dir.create("~/Rlibs", showWarnings = FALSE)
.libPaths("~/Rlibs")
install.packages(c(...), lib = "~/Rlibs")
```

`dplyr`, `readr` and `svglite` have system-library dependencies (libcurl,
libxml2, libssl, fontconfig, freetype). Most managed cluster modules provide
these; on a vanilla Ubuntu/Debian box you may need:

```bash
sudo apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev \
                        libfontconfig1-dev libfreetype6-dev libpng-dev
```
