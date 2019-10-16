context("Test Plotting Functions")
library("MotifSig")

test_that(desc = "Testing VCF Input", {
  sample_df <- readRDS(system.file("testdata", "sample_df.rds", package =
                           "MotifSig"))
  p <- MotifSig::plot_full(sample_df)
  expect_true(ggplot2::is.ggplot(p))
})
