# AB-PGYI-R-Model-Compilation

This repository contains R code originally used to compile output tables from the Alberta Provincial Growth and Yield Initiative (PGYI) database.  A small package structure has been added to expose helper functions for calculating individual tree metrics such as basal area, volume, biomass and missing diameter/height estimates.

## Package Installation

```R
# install from local clone
remotes::install_local(".")
```

## Functions

The `ABTreeMetrics` package provides the following utilities:

- `calc_basal_area(dbh)` – basal area of a tree (m²).
- `calc_qmd(dbh)` – quadratic mean diameter from a vector of DBH values.
- `calc_volume(table, merch = "total")` – stem volume using the PGYI taper equations. Merchantability can be `"total"`, `"1307"`, `"1510"` or `"blank"`. Column names for species, height, dbh and natural subregion can be overridden.
- `calc_biomass(species, dbh, height, vol_1307)` – estimate tree biomass.
- `predict_height(dbh, species, natsub)` – estimate tree height from DBH.
- `predict_dbh(height, species, natsub)` – estimate DBH from height.

The lookup tables used in these functions are found under `GYPSY/GYPSY data/lookup` and are read automatically at runtime.

```R
library(ABTreeMetrics)

# example
calc_basal_area(25)                    # single tree BA
calc_qmd(c(20, 25, 30))                # QMD
```

These functions are distilled from the original scripts located in `GYPSY/GYPSY subscripts`.
