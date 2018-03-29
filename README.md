
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
sample.info <- archs4_sample_info(ids)
head(sample.info)
#> # A tibble: 6 x 7
#>   series_id sample_id  query_type sample_h5idx organism sample_title
#>   <chr>     <chr>      <chr>             <int> <chr>    <chr>       
#> 1 GSE89189  GSM2360252 series            69074 human    10318X2     
#> 2 GSE89189  GSM2360253 series            69075 human    7028X2      
#> 3 GSE89189  GSM2360254 series            69076 human    x2-1        
#> 4 GSE89189  GSM2360255 series            69077 human    x2-2        
#> 5 GSE89189  GSM2360256 series            69078 human    x2-3        
#> 6 GSE89189  GSM2360257 series            69079 human    x2-4        
#> # ... with 1 more variable: sample_name <chr>
```

You can materialize a `DGEList` from a series of series and sample-level identifiers, as well. For instance, to create a DGEList from the [Zhang et al. (Barres) transcriptome database](http://www.jneurosci.org/content/34/36/11929.long). You would identify its GEO series identifier ("[GSE52564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52564)") and do this:

``` r
y <- archs4::as.DGEList("GSE52564", feature_type = "gene")
```

That command takes about a second to run on recent SSD-driven laptop.
