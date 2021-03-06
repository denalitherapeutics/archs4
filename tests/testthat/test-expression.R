context("Expression")

if (!exists("a4")) {
  a4 <- Archs4Repository()
}

test_that("fetch_expression grabs the right count data", {
  gnz <- c("IGFL3", "GSTA4", "MRPS21")
  h5.fn <- file_path(a4, "human_gene")
  h5.idx <- match(gnz, rhdf5::h5read(h5.fn, "meta/genes"))

  for (i in seq(h5.idx)) {
    gene <- gnz[i]
    h5idx <- h5.idx[i]

    counts <- rhdf5::h5read(h5.fn, "data/expression", list(h5idx, NULL))
    counts <- as.vector(counts)
    names(counts) <- rhdf5::h5read(h5.fn, "meta/Sample_geo_accession")
    counts <- counts[!is.na(counts)]

    samples <- sample_table(a4) %>%
      filter(!is.na(libsize), !is.na(normfactor)) %>%
      select(series_id, sample_id) %>%
      semi_join(tibble(sample_id = names(counts)), by = "sample_id") %>%
      distinct(sample_id, .keep_all = TRUE)

    counts <- counts[samples$sample_id]

    res <- fetch_expression(a4, gene, samples = samples,
                            feature_type = "gene", source = "human")
    sample_key <- paste(res$series_id, res$sample_id, sep = "_")
    expect_true(sum(duplicated(sample_key)) == 0, info = gene)
    expect_setequal(names(counts), res$sample_id)
    xref <- match(names(counts), res$sample_id)
    expect_equal(unname(counts), res$count[xref])
  }
})
