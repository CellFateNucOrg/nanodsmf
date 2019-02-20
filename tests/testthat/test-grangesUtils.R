library(testthat)
library(nanodsmf)

context("The DNAstringset2gr function")
test_that("The conversion works",{
  expect_equal(DNAstringset2gr(genomeFile),genomeGR)
})



context("The BSgenome2gr function")
test_that("The conversion works",{
  expect_equal(BSgenome2gr(BSgenome.Celegans.UCSC.ce11::Celegans),genomeGR)
})
