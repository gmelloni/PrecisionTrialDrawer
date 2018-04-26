context("Testing CancerPanel Class")

#Testing if the library "testthat" works correctly
test_that("testthat works", {
  expect_equal(1+1, 2)
  expect_equal(1+2, 3)
  })

#create fake test_panel for simple initial testing
test_panel <- data.frame(drug=""
                         , gene_symbol="MAP2K1"
                         , alteration="SNV"
                         , exact_alteration=""
                         , mutation_specification=""
                         , group=""
                         )

#newCancerPanel(test_panel)

test_that("Test class errors", {
  #expect_error(newCancerPanel(test_panel), info="wrong colname for the panel data.frame")
})


#Example dataset
# drug	gene_symbol	alteration	exact_alteration	mutation_specification	group
# Selumetinib	MAP2K1	SNV	amino_acid_position	68-361	group1
# Selumetinib	MAP2K1	CNA	amplification	NA	group1
# Selumetinib	MAP2K2	SNV	NA	NA	group1
# Selumetinib	MAP2K2	CNA	amplification	NA	group1
# Idelalisib	PIK3CA	SNV	NA	NA	group1
# Idelalisib	PIK3CA	CNA	amplification	NA	group1
# Trastuzumab	ERBB2	CNA	amplification	NA	group1
# Trastuzumab	ERBB2	expression	up	NA	group1
# Vemurafenib	BRAF	SNV	amino_acid_variant	V600E	group1
# Vemurafenib	BRAF	SNV	amino_acid_variant	V600K	group1
# Olaparib	BRCA1	CNA	deletion	NA	group1
# Olaparib	BRCA1	SNV	NA	NA	group1
# Olaparib	BRCA2	CNA	deletion	NA	group1
# Olaparib	BRCA2	SNV	NA	NA	group1
# Vandetanib	ERBB3	SNV	NA	NA	group2
# Vandetanib	ERBB2	SNV	NA	NA	group2
# Vandetanib	EGFR	SNV	NA	NA	group2
# Vandetanib	ABL2	SNV	NA	NA	group2
# Vandetanib	ABL1	SNV	NA	NA	group2
# Vandetanib	KDR	SNV	NA	NA	group2
# Vandetanib	RIT1	SNV	NA	NA	group2
# Vandetanib	KIT	SNV	NA	NA	group2
# ""	KRAS	SNV	genomic_position	12:25358500-25359500	group2
# ""	KRAS	SNV	genomic_position	12:25359000-25359000	group2
# ""	KRAS	SNV	dbSNP_rs	rs121913535	group2
# ""	KRAS	SNV	genomic_position	12:25398225-25398300	group2
# ""	HRAS	SNV	amino_acid_position	5-165	group2
# ""	NRAS	SNV	amino_acid_position	5-165	group2
# ""	NRAS	SNV	amino_acid_position	150-189	group2
# Erlotinib	EGFR	SNV	mutation_type	missense	group2
# ""	EML4-ALK	fusion	""	""	group2
# ""	PML	fusion	""	""	group2


#test_panel = read.csv("data/test.csv", sep=";")
#newCancerPanel(test_panel)
#panel@arguments$panel

################################################################################
## THINGS TO BE tested: 
##
## - Submit empty panel
## - TEST FOR WRONG COLUMN NAME is detected
## - TEST the 6th colum is added correctly
## - TEST if exact_alteration and mutation_specificaiton are both specified or not. 


