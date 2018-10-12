context("Testing class Method - getAlterations()")
######################################################################

#create fake test_panel for simple initial testing
# test_design <- data.frame(drug=""
#                          , gene_symbol="MAP2K1"
#                          , alteration="SNV"
#                          , exact_alteration=""
#                          , mutation_specification=""
#                          , group=""
# )
# 
# test_panel <- newCancerPanel(test_design)
data(cpObj2)
test_panel <- cpObj2

######################################################################
# START TESTS
######################################################################
test_that("getAlterations", {
  
  #------------------------------------------
  # General tests on I/O
  #------------------------------------------
  #expect_error(getAlterations(test_panel, tumor_type="brca"), "*Error: could not find function*")
  
  #------------------------------------------
  # Test tumor type parameter
  #------------------------------------------
  expect_error(getAlterations(test_panel, tumor_type=""), "*You must specify one or more tumor types or tumor studies*")
  #generate mismatches
  expect_error(getAlterations(test_panel, tumor_type="asdasdad"), "*The following tumor types or cancer studies are not available*")
  expect_error(getAlterations(test_panel, tumor_type=c("asdasdad", "brca")), "*The following tumor types or cancer studies are not available*")
  expect_error(getAlterations(test_panel, tumor_type=c("asdasdad", "")), "*The following tumor types or cancer studies are not available*")
  expect_error(getAlterations(test_panel, tumor_type=c("brca", "")), "*The following tumor types or cancer studies are not available*")
  #mixing tumortype with cancer studies
  expect_error(getAlterations(test_panel, tumor_type=c("brca", "brca_tcga_pub")), "*You probably mixed cancer studies and tumor types*")
  
  #this tests takes too long to run
  #getAlterations(test_panel, tumor_type="all_tumors")
})
  

#test repo

#test expression z-score

#test fusions

#test SNV

#test CNA

#test expression




