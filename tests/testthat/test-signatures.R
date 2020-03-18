context("Test Find Signature Functions")
library("BAGEL")

test_that(desc = "Inputs are correct", {
  incomplete_bagel <- readRDS(system.file("testdata", "incomplete-bagel.rds",
                          package = "BAGEL"))

  expect_error(find_signatures(incomplete_bagel, "SNV96"),
               regexp = "count_tables")
  expect_error(infer_signatures(incomplete_bagel, "SNV96"),
               regexp = "count_tables")

  expect_error(find_signatures(data.frame(0)), regexp = "bagel object")
  expect_error(infer_signatures(data.frame(0)), regexp = "bagel object")
})
