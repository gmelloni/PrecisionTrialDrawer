context("Testing helper function - .checkDataAvailability()")
######################################################################


tumor_type=""
genProfile="Mutations"

######################################################################
# START TESTS
######################################################################
test_that("Testing .checkDataAvailability", {

  #Testing I/O
  expect_equal(.checkDataAvailability(tumor_type = tumor_type, genProfile = genProfile)
               , list(), label="Exection an empty list"
  )
 
  # expect_equal(.checkDataAvailability(tumor_type="lgggbm_tcga_pub" , genProfile="Mutations")
  #              , list(), label="Exection an empty list"
  # )
  # 
   
})

#.checkDataAvailability(tumor_type="lgggbm_tcga_pub" , genProfile="Mutations")