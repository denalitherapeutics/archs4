#' Retrieve expression data for genes/transcripts across samples
#'
#' Returns the estimated counts for the genes/transcripts enumerated in the
#' `features` table for the samples enumerated in the `samples` table.
#'
#' Note that the values returned are simply estimated counts. They are not
#' normalized by sequencing depth. For now, the only use for this function is
#' to compare how ratios of genes compare across samples.
#'
#' @export
#' @importFrom edgeR cpm
#' @importFrom reshape2 melt
#'
#' @param a4 An `Archs4Repository`
#' @param features a `tibble` of feature descriptors (as returned from
#'   [feature_info()]). Really this table simply needs the following columns:
#'   * `"a4name"` or `"gene_id" or `;
#'   * `"feature_type"`: (gene or transcript)
#'   * `"source"` (human or mouse),;
#' @param samples a samples identifier: ie. a tibble with series_id and
#'   sample_id columns.
#' @param type "gene" or "transcript"-level data?
#' @param source "mouse" or "human" pass this in explicitly if there is no
#'   "source" column in your `features` or `samples` data.frame
#' @param feature_meta additional columns from
#'   `feature_info(a4, feature_type = feature_type)`
#' @return a data.frame of expression data. Each row is an observation
#'   (one gene one sample)
#' @examples
#' a4 <- Archs4Repository()
#' gnz <- feature_lookup(a4, c("CFAP65", "PECR"), "gene")
#' gexprs <- fetch_expression(a4, gnz)
fetch_expression <- function(a4, features, samples = NULL,
                             feature_type = "gene", source = NULL,
                             feature_meta = "symbol",
                             sample_meta = c("Sample_title", "Sample_source_name_ch1"),
                             ...) {
  assert_class(a4, "Archs4Repository")
  assert_class(features, "data.frame")
  source <- .extract_source(features, samples, source)

  # Extract Features
  if (is.character(features)) {
    features <- feature_lookup(a4, features, feature_type = feature_type)
  }
  assert_data_frame(features)
  assert_subset(c("ensembl_id", "feature_type"), colnames(features))
  features <- distinct(features, h5idx, .keep_all = TRUE)

  # Extract samples
  if (is.null(samples)) {
    samples <- a4 %>%
      sample_table(feature_type = type) %>%
      filter(organism == source) %>%
      select(series_id, sample_id)
  }
  assert_data_frame(samples)
  assert_subset(c("series_id", "sample_id"), colnames(samples))
  sids <- unique(samples$sample_id)

  y <- as.DGEList(a4, sids, features = features, feature_type = feature_type)

  counts <- reshape2::melt(y$counts)
  counts$Var1 <- as.character(counts$Var1)
  counts$Var2 <- as.character(counts$Var2)
  colnames(counts) <- c("ensembl_id", "sample_id", "count")

  cpms <- reshape2::melt(edgeR::cpm(y, prior.count = 5, log = TRUE))
  counts$cpm <- cpms[[3]]

  out <- left_join(samples, counts, by = c("sample_id"))

  if (is.character(feature_meta)) {
    fi.meta <- feature_info(a4, feature_type = feature_type, source = source)
    feature_meta <- unique(c("ensembl_id", feature_meta))
    fi.meta <- fi.meta[, colnames(fi.meta) %in% feature_meta, drop = FALSE]
    if (ncol(fi.meta) > 1L) {
      out <- left_join(out, fi.meta, by = "ensembl_id")
    }
  }

  if (is.character(sample_meta)) {
    si <- tryCatch({
      sample_info(a4, sids, sample_meta, check_missing_samples = FALSE)
    }, error = function(e) NULL)
    if (is.data.frame(si)) {
      scols <- unique(c("sample_id", sample_meta))
      si <- si[, intersect(scols, colnames(si)), drop = FALSE]
      out <- left_join(out, si, by = "sample_id")
    } else {
      warning("There was an error fetching the sample metadata requested")
    }
  }

  out
}

.extract_source <- function(features, samples, source) {
  if (is.null(source)) {
    source <- unique(c(features$source))
    if (is.data.frame(samples)) {
      source <- unique(c(source, samples$source))
    }
    source <- source[!is.null(source)]
    if (length(source) == 0) {
      stop("Don't know where to fetch sources from")
    }
  }
  if (is.character(source)) {
    if (length(source) != 1L) {
      stop("Can't specificy multiple organisms for source -- one at a time!")
    }
  } else {
    stop("source expected to be a character by know")
  }
  if (!source %in% archs4_sources()) {
    stop("Unknown source:", source)
  }
  source
}

#' Retrieves (organism) source of gene/transcript identifiers
#'
#' Parses ensembl identifiers and determines if they are for genes, or
#' transcripts, as well as the organism they should belong to (human, mouse)
#'
#' @param x A vector of identifiers
#' @
feature_source <- function(x, ...) {

}

# library size and normalization factor estimation =============================
#
# The code here was adapted from the code in edgeR.
# We want to calculate and store library sizes and normalization factors
# for an Archs4Repository, so we can calculate things like "cpm" on the fly ...
# the super-fly.

# Default top-level TMM normalization parameters:
#
#   logratioTrim = 0.3
#   sumTrim = 0.05
#   doWeighting = TRUE
#   Acutoff = -1e10
#   p = 0.75
#
# Main TMM Calculation
#   f75 <- .calcFactorQuantile(data=x, lib.size=lib.size, p=0.75)
#   refColumn <- which.min(abs(f75-mean(f75)))
#   f <- rep(NA,ncol(x)) # holds norm factors
#   for(i in 1:ncol(x))
#     f[i] <- .calcFactorWeighted(obs=x[,i], ref=x[,refColumn],
#                                 libsize.obs=lib.size[i],
#                                 libsize.ref=lib.size[refColumn],
#                                 logratioTrim=logratioTrim,
#                                 sumTrim=sumTrim,
#                                 doWeighting=doWeighting,
#                                 Acutoff=Acutoff)
#
# .calcFactorQuantile
# 	y <- t(t(data)/lib.size)           # divides each column by its library size
#   f <- apply(y,2,function(x) quantile(x,p=p)) # find value at upper quartile

# .expr.datasets <- c("human_gene",       "mouse_gene",
#                     "human_transcript", "mouse_transcript")
.expr.datasets <- c("human_gene", "mouse_gene")

#' Estimate normalization factors for datasets in the Archs4Repository
#'
#' This function will serialize the results of the library size and
#' normalization factors into files inside `datadir(a4)`
#'
#' @export
#' @param a4 The `Arcsh4Repository`
#' @return Invisibly returns a tibble of the the library size and normalization
#'   factors for the expression data in `a4` (invisibly).
estimate_repository_norm_factors <- function(a4) {
  info <- lapply(.expr.datasets, function(key) {
    estimate_norm_factors(a4, key)
  })
  invisible(info)
}

#' Calculates library size and norm factors for a specific dataset
#'
#' This function will serialize the results of the library size and
#' normalization factors into files inside `datadir(a4)`
#'
#'
estimate_norm_factors <- function(a4, key, n = 500, logratioTrim = 0.3,
                                  sumTrim = 0.05, doWeighting = TRUE,
                                  Acutoff = -1e10, p = 0.75, ...) {
  if (FALSE) {
    key = "mouse_gene"
    n = 500; logratioTrim = 0.3; sumTrim = 0.05; doWeighting = TRUE;
    Acutoff = -1e10; p = 0.75
  }

  key <- match.arg(key, .expr.datasets)
  h5.fn <- file_path(a4, key)
  fn.out <- paste0(key, "-normfactors.csv")
  fn.out <- file.path(datadir(a4), fn.out)

  # Almost could have saved a lot of time using the v4++ matrices since they
  # include a "meta/reads_aligned" vector, however to use TMM, we still need
  # to calculate the percentiles ...
  # sinfo <- tibble(
  #   series_id = trimws(h5read(h5.fn, "meta/Sample_series_id")),
  #   sample_id = trimws(h5read(h5.fn, "meta/Sample_geo_accession")),
  #   libsize = h5read(h5.fn, "meta/reads_aligned"))
  sinfo <- .crankLibSize(h5.fn, n = n, p = p)
  fn.tmp <- sub("\\.csv$", "-tmp.csv", fn.out)
  write.csv(sinfo, fn.tmp, row.names = FALSE)

  # Identify the sample to use as the reference profile
  refColumn <- which.min(abs(sinfo$quantile - mean(sinfo$quantile, na.rm = TRUE)))
  ref.vals <- h5read(h5.fn, "data/expression", list(NULL, refColumn))

  # Run a TMM for all samples against the reference profile
  nf <- .crankNormFactors(h5.fn, sinfo, ref.vals, n = n,
                          logratioTrim = logratioTrim, sumTrim = sumTrim,
                          doWeighting = doWeighting, Acutoff = Acutoff)
  sinfo$normfactor <- nf

  write.csv(sinfo, fn.out, row.names = FALSE)
  invisible(sinfo)
}

# before the v4 gene-level count matrices, the read counts per sample were not
# included, so we were calculating them manually.
.crankLibSize <- function(h5.fn, n = 500, p = 0.75) {
  message("Calculating library sizes ...")
  sinfo <- tibble(
    series_id = trimws(h5read(h5.fn, "meta/Sample_series_id")),
    sample_id = trimws(h5read(h5.fn, "meta/Sample_geo_accession")))

  n.batches <- nrow(sinfo) / n
  remainder <- nrow(sinfo) %% n
  if (remainder > 0) n.batches <- n.batches + 1
  n.batches <- floor(n.batches)

  offset <- 0

  res <- vector("list", n.batches)
  for (i in 1:n.batches) {
    idx.start <- offset + 1
    idx.end <- idx.start + (if (i != n.batches) n else remainder) - 1L
    idxs <- idx.start:idx.end
    offset <- idx.end

    # lib size and percentile
    dat <- h5read(h5.fn, "data/expression", list(NULL, idxs))
    lib.size <- colSums(dat)
    dat <- t(t(dat)/lib.size)
    # Somehow some samples had NA counts, so we need to protect from that
    fq <- apply(dat, 2, function(x) quantile(x, p = p, na.rm = TRUE))
    res[[i]] <- data.frame(index = idxs, libsize = lib.size, quantile = fq)
  }
  res <- bind_rows(res)
  bind_cols(sinfo, res)
}

.crankNormFactors <- function(h5.fn, sinfo, ref.vals, n = n,
                              logratioTrim = logratioTrim, sumTrim = sumTrim,
                              doWeighting = doWeighting, Acutoff = Acutoff) {
  message("Calculating normalization factors ...")
  ref.size <- sum(ref.vals)
  # Figure out how to batch iterate over our data, we now have to reload the
  # data to calculate the normalization factors.
  n.batches <- nrow(sinfo) / n
  remainder <- nrow(sinfo) %% n
  if (remainder > 0) n.batches <- n.batches + 1
  n.batches <- floor(n.batches)

  # Load batches again and calc norm factors
  nf <- numeric(nrow(sinfo))
  offset <- 0

  for (i in 1:n.batches) {
    idx.start <- offset + 1
    idx.end <- idx.start + (if (i != n.batches) n else remainder) - 1L
    idxs <- idx.start:idx.end
    offset <- idx.end

    # lib size and percentile
    dat <- h5read(h5.fn, "data/expression", list(NULL, idxs))
    for (j in 1:ncol(dat)) {
      idx <- idxs[j]
      vals <- dat[,j]
      if (any(is.na(vals))) {
        nf[idx] <- NA_real_
      } else {
        nf[idx] <- .calcFactorWeighted(obs=vals, ref=ref.vals,
                                       libsize.obs=sum(vals),
                                       libsize.ref=ref.size,
                                       logratioTrim=logratioTrim,
                                       sumTrim=sumTrim,
                                       doWeighting=doWeighting,
                                       Acutoff=Acutoff)
      }
    }
  }

  nf
}

# Calculates the TMM between two libraries.
#
# This is lifted straight from edgeR_3.22.2. The only reason it is included here
# is because it is not exported (rightfully) from the package.
.calcFactorWeighted <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL,
                                logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,
                                Acutoff=-1e10) {
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
  if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

  logR <- log2((obs/nO)/(ref/nR))			# log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2	# absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref	 # estimated asymptotic variance

  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if(max(abs(logR)) < 1e-6) return(1)

  #	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS

  #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #	a fix from leonardo ivan almonacid cardenas, since rank() can return
  #	non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

  if(doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  else
    f <- mean(logR[keep], na.rm=TRUE)

  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}
