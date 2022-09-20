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
                            , genProfile=c("MUTATION" , "COPY_NUMBER_ALTERATION" 
                              , "MRNA_EXPRESSION")
                            ){
    #select only first element of the list
    genProfile <- genProfile[1]
  
    #get cbioportal datasets
    mycgds <- cBioPortalData::cBioPortal(
      hostname = "www.cbioportal.org",
      protocol = "https",
      api. = "/api/api-docs")
    allCanStudy <- cBioPortalData::getStudies(mycgds)[ , c("cancerTypeId" , "studyId")]
    allCanStudy <- unique(allCanStudy)
    allCanStudy$tumor_type <- allCanStudy$cancerTypeId
    #---------------------------------------------------------------------------
    # FIND WHAT TUMOR INFO REQUIRED
    #---------------------------------------------------------------------------
    # All tumors
    if(tumor_type[1]=="all_tumors") {
        chosenTumors <- allCanStudy[,1,drop=TRUE]
    } else {
        #IF tumor_type contains "tumor_types"
        chosenTumors <- allCanStudy[ 
            allCanStudy$tumor_type %in% tumor_type , 1 , drop = TRUE]
        #IF tumor_tupe contains Cancer Studies 
        if(length(chosenTumors)==0){
            chosenTumors <- allCanStudy[ 
                allCanStudy$studyId %in% tumor_type, 1 , drop = TRUE]
            names(chosenTumors) <- allCanStudy[ 
              allCanStudy$studyId %in% tumor_type, 2 , drop = TRUE]
        } else {
          names(chosenTumors) <- allCanStudy[ 
            allCanStudy$tumor_type %in% tumor_type , 2 , drop = TRUE]
        }
    }
    # TEMPORARY!!!
    # pancan <- allCanStudy[,1][ grepl("_pan_can_atlas_2018" 
    #   , allCanStudy[,1])]
    # chosenTumors <- setdiff(chosenTumors , pancan)
    #---------------------------------------------------------------------------
    # Verify Availability
    #---------------------------------------------------------------------------
    out_double <- lapply(names(chosenTumors) , function(i){
        geneticProfile <- cBioPortalData::molecularProfiles(mycgds, i)
        if(genProfile=="MUTATION"){
          geneticProfile2 <- geneticProfile[ grepl(genProfile 
                                 , geneticProfile$molecularAlterationType , ignore.case=TRUE) & geneticProfile$datatype == "MAF"
                                , "molecularProfileId" , drop = TRUE]
        } else if(genProfile=="COPY_NUMBER_ALTERATION"){
          geneticProfile2 <- geneticProfile[ grepl(genProfile 
                                                  , geneticProfile$molecularAlterationType , ignore.case=TRUE) & geneticProfile$datatype == "DISCRETE"
                                            , "molecularProfileId" , drop = TRUE]
        } else if(genProfile=="MRNA_EXPRESSION"){
          geneticProfile2 <- geneticProfile[ grepl(genProfile , geneticProfile$molecularAlterationType , ignore.case=TRUE) & 
                                                    geneticProfile$datatype == "Z-SCORE" &
                                                    grepl("mRNA" , geneticProfile$name) &
                                                    grepl("RNA Seq V2 RSEM" , geneticProfile$name)
                                                    # ( grepl("RNA Seq V2 RSEM" , geneticProfile$name) | grepl("\\(RPKM, log2 transformed\\)" , geneticProfile$name))
                                            , "molecularProfileId" , drop = TRUE]
          
        }
    })
    return(out_double)
}