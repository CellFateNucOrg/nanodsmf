library(testthat)
library(nanodsmf)

context("The DNAstringset2gr function")
test_that("The conversion works",{
  expect_equal(DNAstringset2gr(genomeFile),genomeGR)
})
