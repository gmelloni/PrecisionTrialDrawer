################################################################################
#Helper Function
################################################################################
# The two functions in the present file allow the user to display the codes that
# needs to be implement for calling getAlteration() from the main package. More 
# specifically, using the "codes" here shown, we can select the population of 
# reference that will be used for running the simulations.
#
# Two functions are available:
#
#     1) showTumorType()
#           Aggregate all the studies by tumor type
#     2) showCancerStudy()
#           what studies to include in the analysis
#

################################################################################
# Show named vector of all the available tumor type as TCGA code + description
# Names are fetched from cbioportal, formated and output
# ------------------------------------------------------------------------------
showTumorType <- function() {
  # create CGDS object
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
  # Fetch cancer study from cbioportal
  all_cancer_studies <- cgdsr::getCancerStudies(mycgds)[,c(1,2)]
  
  # From "cancer_study_id" field: 
  #      1) spit by "_" and 
  #      2) select the first element of the split. 
  #          This will correspond to the TCGA disease code
  # From the "name" field: 
  #      1) Split string to get the disease name
  #      2) trim unwanted spaces
  all_cancer_studies2 <- unique(
    data.frame(
      Code=vapply(all_cancer_studies$cancer_study_id
                  , function(x) {
                    strsplit(x , "_")[[1]][1]
                  } , character(1))
      , Full_Name=vapply(all_cancer_studies$name
                         , function(x) {
                           .myTrimmer(strsplit(x , "\\(")[[1]][1])
                        } , character(1))
    )
  )
  
  # Since one TCGA symbol can correspond to 1 or more description, reduce the 
  # dataframe with the aggregate function, putting together description from the
  # same TCGA code separated by "|"
  all_cancer_studies3 <- aggregate(Full_Name~Code, all_cancer_studies2
                                   , FUN=function(x) {paste(x , collapse="|")})
  # reformat output
  out <- data.frame(tumor_type = as.character(all_cancer_studies3$Code)
                  , name = all_cancer_studies3$Full_Name
                  , stringsAsFactors = FALSE)
  # out <- as.character(all_cancer_studies3$Code)
  # names(out) <- all_cancer_studies3$Full_Name
  return(out)
}

################################################################################
# With this function we can see all the studies and choose the ones we like
# ------------------------------------------------------------------------------
showCancerStudy <- function(tumor_type=NULL) {
  # create CGDS object
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
  # Fetch cancer study from cbioportal
  all_cancer_studies <- cgdsr::getCancerStudies(mycgds)[,c(1,2)]
  if(is.null(tumor_type))
    return(all_cancer_studies)
  else
    # transform the tumor_type entreis in a grep pattern
    # run a grep on all the entries using sapply
    # unlist and get the matches
    all_cancer_studies[unlist(lapply(paste0("^" , tumor_type , "_")
                            , function(x) grep(x, all_cancer_studies[,1]))),]
}