context("Testing helper functions showTumorType()")

################################################################################
# START TESTS
################################################################################
#showTumorType is a simple function that takes no input and return a list of 
# description and id. There is little test needed here
#test_that("showTumorType()", {
  #Integrity check
  #expect_equal(length(showTumorType()), 58, label="testing length of shoTumorType output")
#})

################################################################################
#showCancerStydy() 
#-------------------------------------------------------------------------------
test_that("showCancerStudy()", {
  #General tests on I/O
  expect_equal(dim(showCancerStudy("")), c(0,2),label="Expecting an empty output due to empty string")
  expect_equal(dim(showCancerStudy("askdjakfiwrgtbaf")), c(0,2),label="Expecting an empty output due to no match")
  expect_equal(dim(showCancerStudy(c("askdjakfiwrgtbaf", ""))), c(0,2),label="Expecting an empty output due to no match")
  #check that output if correct 
  # These tests are not performed anymore since cBioPortal changes frequently
  #expect_equal(dim(showCancerStudy("brca")), c(7,2),label="Expecting 7 rows and 2 cols output")
  #expect_equal(dim(showCancerStudy("luad")), c(5,2),label="Expecting 5 rows and 2 cols output")
  #expect_equal(dim(showCancerStudy(c("brca", "luad"))), c(12,2),label="Expecting 12 rows and 2 cols output")
})