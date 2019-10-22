context("Test Plotting Functions")
library("BAGEL")

test_that(desc = "Testing VCF Input", {
  result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
  p <- plot_samples(result)
  expect_true(ggplot2::is.ggplot(p))
})
