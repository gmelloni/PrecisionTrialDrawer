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
  mycgds <- cBioPortalData::cBioPortal(
    hostname = "www.cbioportal.org",
    protocol = "https",
    api. = "/api/api-docs")
  all_cancer_studies <- cBioPortalData::getStudies(mycgds)
  all_cancer_studies2 <- unique(
    data.frame(
      # Code=sapply(all_cancer_studies$cancer_study_id
      #     , function(x) strsplit(x , "_")[[1]][1])
      Code=all_cancer_studies$cancerTypeId
      , Full_Name=sapply(all_cancer_studies$name
                         , function(x) .myTrimmer(strsplit(x , "\\(")[[1]][1]))
      , studyId = all_cancer_studies$studyId
    ))
  # all_cancer_studies3 <- aggregate(Full_Name~Code, all_cancer_studies2[ , setdiff(colnames(all_cancer_studies2),"studyId")]
  #                                  , FUN=function(x) {paste(x , collapse="|")})
  # all_cancer_studies4 <- aggregate(studyId~Code, all_cancer_studies2[ , setdiff(colnames(all_cancer_studies2),"Full_Name")]
  #                                  , FUN=function(x) {paste(x , collapse="|")})
  # all_cancer_studies5 <- merge(all_cancer_studies3 , all_cancer_studies4 , by = "Code")
  # # reformat output
  # out <- data.frame(tumor_type = as.character(all_cancer_studies3$Code)
  #                 , name = all_cancer_studies3$Full_Name
  #                 , stringsAsFactors = FALSE)
  out <- all_cancer_studies2
  return(out)
}

################################################################################
# With this function we can see all the studies and choose the ones we like
# ------------------------------------------------------------------------------
showCancerStudy <- function(tumor_type=NULL) {
  # create CGDS object
  mycgds <- cBioPortalData::cBioPortal(
    hostname = "www.cbioportal.org",
    protocol = "https",
    api. = "/api/api-docs")
  all_cancer_studies <- cBioPortalData::getStudies(mycgds)
  if(is.null(tumor_type))
    return(all_cancer_studies)
  else
    # transform the tumor_type entreis in a grep pattern
    # run a grep on all the entries using sapply
    # unlist and get the matches
    all_cancer_studies[unlist(lapply(paste0("^" , tumor_type , "_")
                            , function(x) grep(x, all_cancer_studies$cancerTypeId))),]
}