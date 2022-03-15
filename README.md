# clustack: A Stacking Approach to Spatial and Spatio-Temporal Cluster Detection

`clustack` is an R package based on *Spatial and Spatio-Temporal Cluster Detection Using Stacking* (Kamenetsky, Zhu, Gangnon). `clustack` implements a Poisson to stack cluster estimates. Adjustment for overdispersion is available. The number of clusters is selected using (quasi-) information criteria. We provide methods for calculating cell-wise confidence bounds.


Code developed by and repository maintained by M.Kamenetsky.


The main function in this package is `detectclusters()`, which identifies spatial and/or spatio-temporal clusters using stacking. The function `calcbounds()` calculates cell-wise confidence bounds using the Buckland, Burnham and Anderson, and MATA (model-averaged tail area intervals) approaches.


Please report any bugs or constructive tips [here](https://github.com/mkamenet3/clustack/issues).


## Installation

`clustack` was built on `R` version 4.1.2 "Bird Hippie"

Package dependencies:

- Matrix (>= 1.2-17)

Package imports:
-  clusso (>= 0.1)
- ggplot2 (>= 3.2.1)
- dplyr (>= 0.8.3)
- tidyr (>= 1.0.0)
- magrittr (>= 1.5)
- data.table (>= 1.12.6)
- geosphere (>= 1.5-10)
- MASS (>= 7.3-51.4)
- RColorBrewer (>= 1.1-2)
- Rcpp (>= 1.0.3)
- SOAR (>= 0.99-11)
- stringi (>= 1.4.3)
- rmarkdown (>= 2.11)
- kableExtra (>= 1.3.4)
- tibble (>= 3.1.6)
- stats
- grDevices
- methods
- graphics
- utils
- knitr


To download the latest version of `clustack`:

```
library("devtools")
devtools::install_github("mkamenet3/clustack")
```

## Website

For vignettes and more information, please visit the [clustack website](https://mkamenet3.github.io/clustack/)
