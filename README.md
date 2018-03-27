
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

The `archs4` package provides utility functions to query and explore the raw archs4 data files. These files need to have been downloaded locally on the machine using this package and must all belong in a `getOption("archs4.datadir")` directory.

Installation
============

You first need to download the `human_matrix.h5` and `mouse_matrix.h5` and put them into a single directory. The download directory needs to be set in R such that `getOption("archs4.datadir")` returns the path to this directory.

The default value for `getOption("archs4.datadir")` is `~/.archs4data`.

If you are developing this package and builidng the documentation, the build happens in a vanilla R workspace, which won't set your R's `options` if they are in your `~/.Rprofile`. In this case, symlink your archs4 data directory such that `~/.archs4data` points to the directory you picked on your machine.

Usage
=====

After downloading the necessary files, you can get sample information for GEO series or sample identifiers like so:

``` r
library(archs4)

ids <- c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130")
sample.info <- archs4_sample_info(ids, source = "human")
head(sample.info)
#> # A tibble: 6 x 6
#>   series_id sample_id  sample_h5idx sample_title sample_name   query_type
#>   <chr>     <chr>             <int> <chr>        <chr>         <chr>     
#> 1 GSE89189  GSM2360252        69074 10318X2      iPS microglia series    
#> 2 GSE89189  GSM2360253        69075 7028X2       iPS microglia series    
#> 3 GSE89189  GSM2360254        69076 x2-1         iPS microglia series    
#> 4 GSE89189  GSM2360255        69077 x2-2         iPS microglia series    
#> 5 GSE89189  GSM2360256        69078 x2-3         iPS microglia series    
#> 6 GSE89189  GSM2360257        69079 x2-4         iPS microglia series
```

or create a `DGEList` from a GEO series like so:

``` r
y <- archs4::as.DGEList("GSE89189", feature_type = "gene", source = "human")
```
