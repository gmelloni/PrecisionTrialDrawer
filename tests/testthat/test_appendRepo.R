context("Testing appendRepo function")


###############################
# Testing datasets
# ----------------------------
data(panelexample)
mypanel <- newCancerPanel(panelexample)
#data(cpObj)
#mypanel_full <- cpObj
repos_empty <- lapply(cpData(mypanel) , '[' , c('data' , 'Samples'))
mypanel_full <- getAlterations(mypanel, tumor_type = "brca")
#mypanel_full2 <- getAlterations(mypanel, tumor_type = "luad")
repos_full <- lapply(cpData(mypanel_full) , '[' , c('data' , 'Samples'))

################################################################################
# Tests
# ------------------------------------------------------------------------------
# Test for wrong panel input
# fake_panel = c("one","two")
# testthat::expect_error(appendRepo(fake_panel , repos_full)
#                        , "first argument of the function should be an instance of Class CancerPanel")

# Test for empty panel input
testthat::expect_error(appendRepo(mypanel , repos_full), "*Nothing to append to: The CancerPanel you have entered has no data.*")

# Test for wrong repo - not a list
testthat::expect_error(appendRepo(mypanel_full , repos = c("a", "b")),
                       "*repos should be a list*")

# Test for wrong repo - wrong length
testthat::expect_error(appendRepo(mypanel_full , repos = repos_empty),
                       "*repos should be a list of length 4*")

# Test for wrong repo - wrong names()
testthat::expect_error(appendRepo(mypanel_full, 
                                  repos =list("a"="a", "b"="a", "d"="a", "e"="a")),
                       "repos is not correctly formatted. There are alteration types missing in the list")

#test malformed panel: NO ("data", "Samples")
#test malformed repo: NO ("data", "Samples")
wrong_mypanel <- mypanel_full
#remove the 
wrong_mypanel@dataFull$expression <- list()
wrong_mypanel@dataFull$expression<- list(data=mypanel_full@dataFull$expression$data)
# run check, should generate an error
testthat::expect_error(suppressWarnings(suppressMessages(
  appendRepo(wrong_mypanel , repos = wrong_mypanel@dataFull)))
  , "*Each alteration type should contain two elements named data and Samples*"
)

# TEST results from the function
# ------------------------------------------------------------------------------
#warning dupicated Samples
testthat::expect_warning(suppressMessages(appendRepo(mypanel_full , repos = mypanel_full@dataFull)),
               "*Duplicated sample names in*")

#check that the fusions data table is merged
testthat::expect_equal(suppressMessages(suppressWarnings(
  nrow(appendRepo(mypanel_full , repos = mypanel_full@dataFull)@dataFull$fusion$data)))
  , 1
  , label = "Check that the fusions are merged")


#test data merging in copynumber
mypanel_test <- mypanel_full
mypanel_test@dataFull$copynumber$data <- mypanel_full@dataFull$copynumber$data[1:10,]

testthat::expect_equal(suppressMessages(suppressWarnings(
  nrow(appendRepo(mypanel_test , repos = mypanel_test@dataFull)@dataFull$copynumber$data)))
  , 10
  , label = "Check that the copynumbers are merged")
