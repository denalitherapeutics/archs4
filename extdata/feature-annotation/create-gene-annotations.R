# Defines functions to create gene annotations from featurecount results.
#
# Currently we just create feature metadata once and store it here so we can
# easily append it to outgoing `as.DGEList`s.

library(biomaRt)
library(GenomicsStudyDb)
library(checkmate)
library(rprojroot)
library(readr)

root <- find_root(is_r_package)

kBucket <- "com-dnli-ngs"
kFeatureCountsMouse <- "NGS000003/309O/featurecounts/gene_counts.txt"
kFeatureCountsHuman <- "NGS000001/09_49_TG/featurecounts/gene_counts.txt"

# NOTE: When updating the annotation files, make sure you:
# 1. correct their paths (if changed) where they are laoded. This is currently
#    (Feb 14, 2018) only in the .load_results.featurecount funtion. Look for
#    "aug.fn" in the R/aws.s3.R to see what the names are.
# 2. Update the col_types encoding in the read_tsv call. Again this will be
#    near the "aug.fn" parts of aws.s3.R
#  We should think of a more robust way to encode these paths and their columns
#  in the package itself, probably.

#' Create a gene annotation file from a featurecount result
#'
#' We are assuming a lot ... featureCounts was run over an
#' ENSEMBL / GENCODE BASIC gtf and that the ensembl features are versioned.
#' Essentially assuming we are only using gencode-basic annotations for now.
#'
#' @param x a path to a full featureCounts output, or the parsed version of it
#'   that comes from GenomicsStudyDb:::read_featurecount
#' @return a tibble of gene annotations for the features.
create_gene_annotations <- function(x, type = c("gene", "transcript")) {
  type <- match.arg(type)
  if (type == "transcript") {
    stop("Not dealing with transcripts right now.")
  }

  if (is.character(x)) {
    assert_file_exists(x)
    x <- GenomicsStudyDb:::read_featurecount(fn = x)
  }
  fc.info <- GenomicsStudyDb:::featurecount_info(x) # ensures `x` is legit

  # setup the parameters to interact with biomaRt for the given species and
  # feature type
  if ("human" == fc.info$species) {
    ds <- "hsapiens_gene_ensembl"
    symbol <- "hgnc_symbol"
  } else if ("mouse" == fc.info$species) {
    ds <- "mmusculus_gene_ensembl"
    symbol <- "mgi_symbol"
  } else {
    stop("Unknown species: ", fc.info$species)
  }
  fltr <- switch(type,
                 gene = "ensembl_gene_id_version",
                 transcript = "ensembl_transcript_id_version")
  atrs <- c("ensembl_gene_id", "ensembl_gene_id_version",
            "entrezgene", "gene_biotype", symbol)

  # do the query
  ensembl <- useMart("ensembl", dataset = ds)
  info <- getBM(attributes = atrs, filters = fltr, values = x[["feature_id"]],
                mart = ensembl)
  # runaround NSE/tidyeval
  info[["symbol"]] <- info[[symbol]]
  info[["query"]] <- info[[fltr]]

  # resort data by query, symbol, entrezgene so that NAs in symbol and
  # entrez go to bottom of each query
  out <- info %>%
    dplyr::arrange(query, symbol, entrezgene) %>%
    dplyr::distinct(query, .keep_all = TRUE) %>%
    dplyr::select(query, ensembl_gene_id, entrez_id = entrezgene, symbol,
                  gene_biotype)
    mutate(feature_type = "ensembl")
  out
}

# Let's me easily source this script then run code in if (FALSE) block
if (FALSE) {
  # Create mouse annotations
  fcs <- GenomicsStudyDb:::read_featurecount(key = kFeatureCountsMouse,
                                             bucket = kBucket)
  info <- GenomicsStudyDb:::featurecount_info(fcs)
  annos <- create_gene_annotations(fcs, "gene")
  fn.out <- file.path(root, "inst", "extdata", "feature-annotation",
                      sub(".gtf", ".csv.gz", info$annotation))
  readr::write_csv(annos, fn.out)

  # Create human annotations
  fcs <- GenomicsStudyDb:::read_featurecount(key = kFeatureCountsHuman,
                                             bucket = kBucket)
  info <- GenomicsStudyDb:::featurecount_info(fcs)
  annos <- create_gene_annotations(fcs, "gene")
  fn.out <- file.path(root, "inst", "extdata", "feature-annotation",
                      sub(".gtf", ".csv.gz", info$annotation))
  readr::write_csv(annos, fn.out)
}
