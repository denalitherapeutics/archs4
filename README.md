
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

This package provides utility functions to query and explore the data made available through the [ARCHS4 project](https://amp.pharm.mssm.edu/archs4/).

Installation
============

The user will first have to download the four gene and transcript level hdf5 data files made available on the [ARCHS4 downloads](https://amp.pharm.mssm.edu/archs4/download.html) page for the mouse and human data. These data files need to all be collated into the directory specified by `getOption("archs4.datadir")` (`~/.archs4data`, by default).

**NOTE:** If you are developing this package and builidng the documentation, the build happens in a vanilla R workspace, which won't set your R's `options` if they are in your `~/.Rprofile`. In this case, symlink your archs4 data directory such that `~/.archs4data` points to the directory you picked on your machine.

Usage
=====

Please note again that successful usage of this package requires a valid entry for `getOption("archs4.datadir")`, otherwise many functions will require you to explicitly pass in a `datadir` parameter.

``` r
library(archs4)

a4 <- Archs4Repository()
ids <- c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130")
sample.info <- sample_info(a4, ids)
head(sample.info)
#> # A tibble: 6 x 8
#>   series_id sample_id  query_type sample_h5idx_gene sample_h5idx_transcri…
#>   <chr>     <chr>      <chr>                  <int>                  <int>
#> 1 GSE89189  GSM2360252 series                 69074                  60169
#> 2 GSE89189  GSM2360253 series                 69075                  60170
#> 3 GSE89189  GSM2360254 series                 69076                  60173
#> 4 GSE89189  GSM2360255 series                 69077                  60182
#> 5 GSE89189  GSM2360256 series                 69078                  60188
#> 6 GSE89189  GSM2360257 series                 69079                  60163
#> # ... with 3 more variables: organism <chr>, Sample_title <chr>,
#> #   Sample_source_name_ch1 <chr>
```

You can use the `as.DGEList` function to materialize a `edgeR::DGEList` from a vector of series or sample-level identifiers, as long as they all reference the same species.

The most often use-case will likely be to create a `DGEList` for a given study. For instance, the GEO series identifier [`"GSE89189"`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89189) refers to the expression data generated to support the [Abud et al. iPSC-Derived Human Microglia-like Cells ...](https://www.ncbi.nlm.nih.gov/pubmed/28426964) paper.

You can get a `DGEList` of gene-level counts for that study using the command below. This will create a `DGEList` of 27,024 genes across 37 samples in about 1.5 seconds:

``` r
yg <- as.DGEList(a4, "GSE89189", feature_type = "gene")
```

The following command retrieives the 178,135 transcript level counts for this experiment in about 1.5 seconds, as well:

``` r
yt <- as.DGEList(a4, "GSE89189", feature_type = "transcript")
```
