library(vegan)
library(RLdbRDA)

test_that("RLdbRDA works", {
  data(varespec)
  data(varechem)

  expect_no_error(rldbrda(varespec, varechem))
})

test_that("Prepare plot data works", {
  data(varespec)
  data(varechem)

  out <- rldbrda(varespec, varechem)

  expect_no_error(prepare_plot_data(out))
})

test_that("Plot dbRDA works", {
  data(varespec)
  data(varechem)

  out <- rldbrda(varespec, varechem)
  plot_data <- prepare_plot_data(out)

  expect_no_error(plot_dbrda(plot_data))
})

test_that("Validation colnames metadata", {
  data(varespec)
  data(varechem)

  colnames(varechem)[1] <- "N test"

  expect_error(rldbrda(varespec, varechem))
})
