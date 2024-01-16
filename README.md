<!-- badges: start -->
[![R-CMD-check](https://github.com/raeslab/RLdbRDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/raeslab/RLdbRDA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# RLdbRDA: RaesLab distance-based Redundancy Analysis

The RLdbRDA package provides a streamlined workflow for performing distance-based redundancy analysis (db-RDA) as is done in the RaesLab on microbial datasets, enabling researchers to identify and visualize the influence of metadata variables on community composition.

## Installation


LRdbRDA needs to be installed directly from GitHub using devtools. From an R console enter the commands below.


```commandline
library(devtools)
install_github("raeslab/LRdbRDA")
```

If you are using ```renv```, instead use the commands below to install this package.

```commandline
renv::install("raeslab/LRdbRDA")
```


## Usage

### Example

```R
library(RLdbRDA)
library(vegan)

data(varespec)
data(varechem)

out <- rldbrda(varespec, varechem)
out

plot_data <- prepare_plot_data(out)
plot_data

g <- plot_dbrda(plot_data)
g
```

![bar plot showing the single and cumulative effect of various features on the varespec dataset included in vegan](./docs/img/rldbrda_example_output.png)


