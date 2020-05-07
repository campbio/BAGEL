context("Test Find Signature Functions")
library("BAGEL")

test_that(desc = "Inputs are correct", {
  incomplete_bagel <- readRDS(system.file("testdata", "incomplete-bagel.rds",
                          package = "BAGEL"))

  expect_error(discover_signatures(incomplete_bagel, "SNV96"),
               regexp = "count_tables")
  expect_error(predict_exposure(incomplete_bagel, "SNV96"),
               regexp = "count_tables")

  expect_error(discover_signatures(data.frame(0)), regexp = "bagel object")
  expect_error(predict_exposure(data.frame(0)), regexp = "bagel object")
})
