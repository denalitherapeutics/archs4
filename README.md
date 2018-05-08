
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

The `archs4` package provides utility functions to query and explore the expression profiling data made available through the [ARCHS4 project](https://amp.pharm.mssm.edu/archs4/), which is described in the following publication:

[Massive mining of publicly available RNA-seq data from human and mouse](https://www.nature.com/articles/s41467-018-03751-6).

Because this package requires the user to download a number of data files that are external to the package, the [installation instructions](#installation) are *a bit* more involved than other R packages, and we leave them for [the end of this document](#installation).

Usage
=====

After [successful installation](#installation) of this package, you can query the series and samples included in the ARCHS4 repository, as well as materialize the expresion data into well-known bioconductor assay containers for downstream analysis.

To query GEO series and samples, you can use the `sample_info` function:

``` r
library(archs4)

a4 <- Archs4Repository()
ids <- c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130")
sample.info <- sample_info(a4, ids)
head(sample.info)
#> # A tibble: 6 x 8
#>   series_id sample_id  Sample_title Sample_source_name_ch1 query_type
#>   <chr>     <chr>      <chr>        <chr>                  <chr>     
#> 1 GSE89189  GSM2360252 10318X2      iPS microglia          series    
#> 2 GSE89189  GSM2360253 7028X2       iPS microglia          series    
#> 3 GSE89189  GSM2360254 x2-1         iPS microglia          series    
#> 4 GSE89189  GSM2360255 x2-2         iPS microglia          series    
#> 5 GSE89189  GSM2360256 x2-3         iPS microglia          series    
#> 6 GSE89189  GSM2360257 x2-4         iPS microglia          series    
#> # ... with 3 more variables: sample_h5idx_gene <int>,
#> #   sample_h5idx_transcript <int>, organism <chr>
```

You can use the `as.DGEList` function to materialize an `edgeR::DGEList` from a an arbitrary number of GEO sample and series identifiers. The only restriction is that the data from the series/samples must all be from the same species.

The most often use-case will likely be to create a `DGEList` for a given study. For instance, the GEO series identifier [`"GSE89189"`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89189) refers to the expression data generated to support the [Abud et al. iPSC-Derived Human Microglia-like Cells ...](https://www.ncbi.nlm.nih.gov/pubmed/28426964) paper.

Creating a `DGEList` from this study will create an object with 27,024 genes across 37 samples in about 1.5 seconds:

``` r
yg <- as.DGEList(a4, "GSE89189", feature_type = "gene")
```

The following command retrieves the 178,135 transcript level counts for this experiment in about 1.5 seconds, as well:

``` r
yt <- as.DGEList(a4, "GSE89189", feature_type = "transcript")
```

Installation
============

The installation of the `archs4` package is a bit more involved than a standard package installation and can be roughly broken down into three steps.

1.  Install the R package along with its dependencies.
2.  Download a number of (large) data files into a specific folder.
3.  Generate metadata from the files downloaded in (2) for downstream use.

We will walk you through each step in this section.

R Package Installation
----------------------

The `arcsh4` package depends on other packages that are available through both [CRAN](https://cran.r-project.org/) and [Bioconductor](http://bioconductor.org/). For that reason, we will use the [`BiocInstaller::biocLite()`](https://www.bioconductor.org/install/#why-biocLite) function to install this package, which can seamlessly install packages from github, CRAN, and Bioconductor.

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("denalitherapeutics/archs4")
library("archs4")
```

When you first load the `archs4` library, you will notice a startup message telling you that something isn't quite right with your `archs4` installation. The message will look something like this:

    Note that your default archs4 data directory is NOT setup correctly

      * Run `archs4_local_data_dir_validate()` to diagnose
      * Refer to the ARCHS4 Data Download section of the archs4 vignette for more information

    Your default archs4 data directory (`getOption("archs4.datadir")`) is:

      ~/.archs4data

In order for the package to work correctly, you must download a number of files which are enumerated in the [Data File Download](#data-file-download) section below into a single directory. You will then instruct the `archs4` package the path to the directory that holds all of these files by setting the value of R's global `"archs4.datadir`" option to be the path to that directory.

Data File Download
------------------

You will have to create a directory on your filesystem which will hold a number of data files that the `archs4` package depends on. Let's call this directory `$ARCHS4DIR`, which we will define here to be `~/archs4v2data`.

The `archs4` package provides the `archs4_local_data_dir_create()` convenience function which creates this directory and copies over a `meta.yaml` file into that directory. The purpose of this file is to specify the names of the downloaded files that correspond to the human and mouse-level gene and transcript-level data.

``` r
library(archs4)
archs4dir <- "~/archs4v2data"
archs4_local_data_dir_create(archs4dir)
```

Once this directory is created successfully, you will then have to download the following files into it:

-   archs4
    -   [`human_matrix.h5`](https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5): human gene-level counts
    -   [`human_hiseq_transcript_v2.h5`](https://s3.amazonaws.com/mssm-seq-matrix/human_hiseq_transcript_v2.h5): human transcript-level counts
    -   [`mouse_matrix.h5`](https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5): mouse gene-level counts
    -   [`mouse_hiseq_transcript_v2.h5`](https://s3.amazonaws.com/mssm-seq-matrix/mouse_hiseq_transcript_v2.h5): mouse transcript-level counts
-   ensembl
    -   `Homo_sapiens.GRCh38.90.gtf.gz`: gtf used for human transcript annotations <ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz>
    -   `Mus_musculus.GRCm38.90.gtf.gz`: gtf used for mouse transcript annotations <ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz>

The enumerated items above contain links to the files that need to be downloaded. You can right-click on them and select `Save As ...` and instruct your web-browser to save them to your local `$ARCHS4DIR`.

**NOTE**: Most all of the `archs4` functions accept a `datadir` parameter, which should be the path to `$ARCHS4DIR`. For convenience, the default value of this parameter is always set to `getOption("archs4.datadir")`. This means that you can modify your `~/.Rprofile` file to set the value of this option to `"~/archs4v2data"` (for instance), so that the package will always look there by default. If this option is not set in your `~/.Rprofile`, the default value for this option is "~/.archs4data".

Feature-Level Metadata Generation
---------------------------------

The datasets currently made available by the [ARCHS4 Project](https://amp.pharm.mssm.edu/archs4/) only provide minimal feature-level metadata:

-   the features in the gene-level datasets are identified only by their symbol; and
-   only the ensembl transcript id's are provided for the features in the transcript-level datasets

We want to augment these features with richer annotations, such as the ensembl gene identifiers or gene biotypes, for instance.

To make such data generation automatic and easy for the user, once you have downloaded the Ensembl GTF files listed above into the `$ARCHS4DIR`, you can run the `create_augmented_feature_info()` to extract these extra feature-level metadata from the GTF files and store them as tables inside `$ARCHS4DIR` for later use.

``` r
create_augmented_feature_info(archs4dir)
```

This function will load and parse the GTF files from human and mouse, and create gene- and transcript-level `*.csv.gz` files in the `$ARCHS4DIR` which the `archs4` package will then later use downstream.

Once your `$ARCHS4DIR` is setup, you may find it convenient to set the default value for R's global `"archs4.datadir"` option to the `$ARCHS4DIR` directory you just setup. To do so, you can put the following line in your `~/.Rprofile` file:

``` r
options(archs4.datadir = "~/archs4v2data")
```

ARCHS4 Installation Heatlh
--------------------------

Because the installation of this package is a bit more involved than most, we have also provided an `archs4_local_data_dir_validate()` function, which you can run over your `$ARCHS4DIR` in order to check on "the health" of your install.

This function will simply look at your `$ARCHS4DIR` to ensure that the required files are there, and tries to give you helpful error messages if not.

For instance, if the first two files enumerated in the [Data File Download](#data-file-download) section were missing from your `$ARCHS4DIR` (ie. `human_matrix.h5` and `human_hiseq_transcript_v2.h5`), you would be warned that "something isn't right" when you first load the `archs4` package. You could then run the `archs4_local_data_dir_validate()` to see what is wrong:

``` r
archs4_local_data_dir_validate(archs4dir)
#> The following ARCHS4 files are missing, please download them:
#>   * human_matrix.h5: https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5
#>   * human_hiseq_transcript_v2.h5: #> https://s3.amazonaws.com/mssm-seq-matrix/human_hiseq_transcript_v2.h5
```

**NOTE:** If all installation and data download/processing steps have been completed successfully, a call to `archs4_local_data_dir_validate()` will simply return `TRUE`.

Package Development
-------------------

If you are developing this package, you will find that it will be convenient to symlink the package's default `archs4.datadir` path (`~/.arcsh4data`) to the `$ARCHS4DIR` you just setup. This is because often times things like roxygen2 document compilation, unit testing, etc. happen in a vanilla R workspace, which won't run the configuration that is prescribed in your `~/.Rprofile` file.
