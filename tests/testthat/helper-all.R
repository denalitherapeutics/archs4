# These commands are run before the tests and are made available for everyone
# to use.

# Instantiate a universal Archs4Repository because it takes a little while
# (~4s) to construct one.
if (!interactive() && identical(Sys.getenv("NOT_CRAN"), "true")) {
  # Putting it in this block, because helper-* function are run for convenience
  # in several scenarios, including when we:
  #
  #   1. Run devtools::load_all(); and
  #   2. Run devtools::document()
  #
  # In these scenarios, I don't want to instantiated the Archs4Repository, so
  # the if () shoudl only evaluate to true when we are actually *testing*.
  a4 <- Archs4Repository()
}
