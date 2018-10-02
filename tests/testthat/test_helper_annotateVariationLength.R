context("Testing .annotateVariationLength()")
######################################################################

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
test_that(".annotateVariationLength()", {
  #General tests on I/O
  .annotateVariationLength(test_panel, .annotateGeneLength(test_panel$gene_symbol), utr = FALSE, padding_length=0 ) 
  expect_error(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="error"
                                                   , mutation_specification=""
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=0)
               , "*Mutation Specification is not correct at panel line 1*", 
               label="simulate error in the exact_alteration field")
  
  expect_equal(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="amino_acid_position"
                                                   , mutation_specification="300-300"
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=0)$variation_len
               , 3
               , label="Simulate length from amino_acid_position 300-300, to be equal to 1")
  
  expect_equal(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="amino_acid_variant"
                                                   , mutation_specification="V600E"
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=0)$variation_len
               , 3
               , label="simulate length from amino_acid_variant V600E, to be equal to 3")
  
  expect_equal(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="genomic_position"
                                                   , mutation_specification="12:2000-3000"
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=0)$variation_len
               , 1001
               , label="simulate length from RS id, to be equalt to 1001")
  
  expect_equal(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="genomic_variant"
                                                   , mutation_specification=""
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=1)$variation_len
               , 3
               , label="simulate length from genomic_variant with no specifications, to be equal to 3")
  
  expect_equal(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="genomic_variant"
                                                   , mutation_specification=""
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=1)$variation_len
               , 3
               , label="simulate length from genomic_variant with no specifications, to be equal to 3")
  
  expect_equal(.annotateVariationLength(data.frame(drug=""
                                                   , gene_symbol="MAP2K1"
                                                   , alteration="SNV"
                                                   , exact_alteration="dbSNP_rs"
                                                   , mutation_specification=""
                                                   , group="")
                                        , .annotateGeneLength(test_panel$gene_symbol)
                                        , utr = FALSE
                                        , padding_length=0)$variation_len
               , 1
               , label="simulate length from RS id, to be equalt to 1")
})
