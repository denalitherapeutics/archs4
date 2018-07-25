context("Expression Concordance")

# There have been a few occassions where the ARCHS4 gene expression quantitation
# has lined up far from what is expected. For example:
#
# 1. the values in the human v4 gene-level expression matrices lined up poorly
#    with internal estimates of the same, and this concordance was resolved when
#    the v5 gene-level matrices were published.
# 2. The v5 mouse gene-level matrices look a bit off compared to internally
#    processed data, but the v4 versions of the seemed more inline with
#    expectation.

# As a "validation set", a subset of gene level
# `edgeR::cpm(., prior.count = 3, log = TRUE)` values are provided for studies
# processed by:
# 1. A STAR + salmon + GENOCDE ref piline; and
# 2. recount2

test_that("human expression is concordant with expectation", {

})

test_that("mouse expression is concoordant with expectation", {

})