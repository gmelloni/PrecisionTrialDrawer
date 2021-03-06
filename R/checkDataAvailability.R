################################################################################
# Helper Function
################################################################################
#
# Check if data is available from mutations/CNA/expression
#
# INPUT:
#   tumor_type: my list of cancer studies or tumor types 
# OUTPUT:
#   a list with the availabel studies.
#
# By defaults it will download all the possible tumor_types.
.checkDataAvailability <- function(
                            tumor_type="all_tumors"
                            # ,tumor_study=NULL
                            , genProfile=c("mutations$" , "gistic|cna$" 
                              , "_rna_seq_v2_mrna_median_Zscores$")
                            ){
    #select only first element of the list
    genProfile <- genProfile[1]
  
    #get cbioportal datasets
    mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
    all_cancer_studies <- cgdsr::getCancerStudies(mycgds)[,c(1,2)]
    all_cancer_studies$tumor_type <- vapply(
      strsplit(all_cancer_studies[,1] , "_") , '[' , character(1) , 1)
    
    #---------------------------------------------------------------------------
    # FIND WHAT TUMOR INFO REQUIRED
    #---------------------------------------------------------------------------
    # All tumors
    if(tumor_type[1]=="all_tumors") {
        chosenTumors <- all_cancer_studies[,1]
    } else {
        #IF tumor_type contains "tumor_types"
        chosenTumors <- all_cancer_studies[ 
            all_cancer_studies[,'tumor_type'] %in% tumor_type , 1]
        #IF tumor_tupe contains Cancer Studies 
        if(length(chosenTumors)==0){
            chosenTumors <- all_cancer_studies[ 
                all_cancer_studies[,'cancer_study_id'] %in% tumor_type, 1]
        }
    }
    # TEMPORARY!!!
    pancan <- all_cancer_studies[,1][ grepl("_pan_can_atlas_2018" 
      , all_cancer_studies[,1])]
    chosenTumors <- setdiff(chosenTumors , pancan)
    #---------------------------------------------------------------------------
    # Verify Availability
    #---------------------------------------------------------------------------
    out_double <- lapply(chosenTumors , function(i){
        geneticProfile <- cgdsr::getGeneticProfiles(mycgds, i)[ ,c(1,2)]
        geneticProfile <- grep(genProfile 
                , geneticProfile$genetic_profile_id 
                , value=TRUE , ignore.case=TRUE)
        })
    return(out_double)
}