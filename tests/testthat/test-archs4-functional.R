context("ARCHS4 Functional Interface")

if (!exists("a4")) {
  # This is loaded by the testthat/helper-all.R script when testthat is running
  # the unit tests, but included here for convenience when doing interactive
  # test development
  a4 <- Archs4Repository()
}

