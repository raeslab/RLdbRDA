library(vegan)
library(RLdbRDA)

test_that("RLdbRDA works", {
  data(varespec)
  data(varechem)

  expect_no_error(rldbrda(varespec, varechem))
})
