# These commands are run before the tests and are made available for everyone
# to use.

# Instantiate a universal Archs4Repository because it takes a little while
# (~4s) to construct one.
if (!interactive()) {
  # Putting it in this block, because a devtools::load_all() loads the helper
  # functions, and I don't want that to happen.
  a4 <- Archs4Repository()
}
