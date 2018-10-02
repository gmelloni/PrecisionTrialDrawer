context("Testing helper functions")

#create df for testing
test_df = data.frame(col1=factor(c("a", "a","a", "b","b", 1))
                     ,col2= factor(c(1, 3, 2, 5, 3, 3))
)

#create fake test_panel for simple initial testing
test_panel <- data.frame(drug=""
                         , gene_symbol="MAP2K1"
                         , alteration="SNV"
                         , exact_alteration=""
                         , mutation_specification=""
                         , group=""
)


######################################################################
# START TESTS
######################################################################
test_that(".annotateGeneLength()", {
  #General tests on I/O
  expect_error(test_out <- .annotateGeneLength(c("KRAS", "BRC1", "ERBB2")), label="Expected warning due to added typo")
  expect_equal(as.character(.annotateGeneLength(c("KRAS"))), c("KRAS", "570", "1119"))
  expect_error(.annotateGeneLength(), label="Expect argument missing with no default")
  expect_error(.annotateGeneLength(c()), "*Submitted gene list is empty*")
  expect_error(.annotateGeneLength(c("")), "*Submitted gene list contains empty gene names*" , label="Submitted gene list contains empty gene names")
  expect_error(.annotateGeneLength(""),  "*Submitted gene list contains empty gene names*",  label="Submitted gene list contains empty gene names")
  expect_error(.annotateGeneLength(c("KRAS", "BRC1", "ERBB2", "")),  "*Submitted gene list contains empty gene names*",  label="Submitted gene list contains empty gene names")
  expect_error(.annotateGeneLength(c("KRAS", "BRC1", "ERBB2")),  "*The following gene symbols were not found:*")
  expect_error(.annotateGeneLength("BRRR1"),  "*The following gene symbols were not found:*")
  #Test on the gene length
  expect_equal(.annotateGeneLength(c("KRAS", "KRAS", "KRAS", "KRAS")), .annotateGeneLength(c("KRAS")), label="results shuold generate the same length")
})


######################################################################
# .rsToGenomicPosition()
# --------------------------------------------------------------------
test_that(".rsToGenomicPosition", {
  #I/O tests
  expect_equal(.rsToGenomicPosition("rs2396789")$genomic_range, "17:8547-8547", label="Test RS location")
  expect_equal(.rsToGenomicPosition("rs2396789")$rs, "rs2396789", label="Test RS name imported correctly")
  expect_error(.rsToGenomicPosition(c("","rs2396789")), "*Not all the RS ids where found in Biomart*", label="Submitted empty rs id")
  expect_error(.rsToGenomicPosition(""), "*Input rs list was empty*", label="Test RS location")
  expect_error(.rsToGenomicPosition("asdaklsskghiruhgszfga"), "*Not all the RS ids where found in Biomart*", label="Submitted wrong rs id")
})
