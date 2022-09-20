# Given a vector of genes and a vector of tumor type as resulted 
# from showTumorType() or showCancerStudy()
# It returns a nested list composed by n elements, one for each tumor type
# Each element is composed by a CNA df and a vector with all patient 
# barcodes for that tumor type

.getCNA <- function(myGenes=myGenes
                    ,tumor_type="all_tumors"
                    ,block=NULL) {
  mycgds <- cBioPortalData::cBioPortal(
    hostname = "www.cbioportal.org",
    protocol = "https",
    api. = "/api/api-docs")
  allCanStudy <- cBioPortalData::getStudies(mycgds)[ , c("cancerTypeId" , "studyId")]
  allCanStudy <- unique(allCanStudy)
  # allCanStudy$tumor_type <- vapply(strsplit(allCanStudy$cancerTypeId
  #                                       , "_") , '[' , character(1) , 1)
  allCanStudy$tumor_type <- allCanStudy$cancerTypeId
  if(tumor_type[1]=="all_tumors") {
    chosenTumors <- allCanStudy[,1,drop=TRUE]
  } else {
    #FIND All cancer studies associated with the tumor_type specified
    chosenTumors <- allCanStudy[ allCanStudy$tumor_type %in% 
                                   tumor_type, 1,drop=TRUE]
    names(chosenTumors) <- allCanStudy[ allCanStudy$tumor_type %in% 
                                          tumor_type, 2,drop=TRUE]
    #In case we are not looking for tumor ID but cancer studies ID
    if(length(chosenTumors)==0){
      chosenTumors <- allCanStudy[ allCanStudy$studyId %in% 
                                     tumor_type, 1,drop=TRUE]
      names(chosenTumors) <- allCanStudy[ allCanStudy$studyId %in% 
                                            tumor_type, 2,drop=TRUE]
    }
    #if still we have not found anything
    if (length(chosenTumors)==0){
      stop(paste("Could not find the tumor_type nor"
                 ,"the cancer_study_id specified"))
    }
  }
    
  out_double <- lapply(names(chosenTumors) , function(i){
    geneticProfile <- cBioPortalData::molecularProfiles(mycgds, i)
    geneticProfile <- geneticProfile[ grepl("COPY_NUMBER_ALTERATION" , geneticProfile$molecularAlterationType) &
                                        grepl("DISCRETE" , geneticProfile$datatype)
                                        # grepl("LOG2-VALUE" , geneticProfile$datatype) 
                                      , ]
    if(nrow(geneticProfile)==0){
      message(paste("geneticProfile" , i , "has no copynumber data"))
      return( list( out=NULL , patients=NULL) )
    }
    # Fetch the patient list for the specified molecular profile and data type (mutations)
    caseList <- cBioPortalData::sampleLists(mycgds, i)
    # caseList <- caseList[ caseList$category == "all_cases_with_log2_cna_data" , ]
    caseList <- caseList[ caseList$category == "all_cases_with_cna_data" , ]
    if(nrow(caseList)==0){
      message(paste("ProfileId" , i , "has no samples with copynumber"))
      return( list( out=NULL , patients=NULL) )
    }
    # Get the actual case names (e.g. TCGA-AR-A1AR)
    caseListId <- cBioPortalData::getSampleInfo(api = mycgds 
                                                , studyId = geneticProfile$studyId 
                                                , sampleListIds = caseList$sampleListId)$sampleId
    tryCatch(
      cna <- cBioPortalData::molecularData(
        api = mycgds
        , molecularProfileIds = geneticProfile$molecularProfileId
        , entrezGeneIds = as.numeric(names(myGenes))
        , sampleIds = caseListId
      )[[1]]
      , error=function(e) message(paste("Impossible to retrive copynumber from" , i , "study or no copynumber are present on the selected genes"))
    )
    if(!exists("cna")) {
      cna <- NULL
      patients <- NULL
    } else if(nrow(cna) == 0 ){
      cna <- NULL
      patients <- NULL
    } else {
      patients <- unique(ifelse(grepl("^TCGA" , caseListId) 
                         , unlist(lapply(strsplit(caseListId , "-") , function(x) paste(x[seq_len(3)] , collapse="-")))
                         , caseListId))
      cna$genetic_profile_id <- i
      cna$tumor_type <- chosenTumors[i]
      cna$case_id <- ifelse(grepl("^TCGA" , cna$sampleId) 
                             , unlist(lapply(strsplit(cna$sampleId , "-") , function(x) paste(x[seq_len(3)] , collapse="-")))
                             , cna$sampleId)
      cna$gene_symbol <- myGenes[ as.character(cna$entrezGeneId) ]
      cna <- cna[ , c("case_id" , "entrezGeneId" , "gene_symbol" , "value" , "genetic_profile_id" , "tumor_type")]
    }
    # message(paste("Copynumber retrieved from" , i , "study"))
    return( list( out=cna , patients=patients) )
  })
  names(out_double) <- names(chosenTumors)
  return(out_double)
}

.subsetCNA <- function(panel , cna_full)
{
    panel_cna <- panel[ panel$alteration=="CNA" , ]
    cna_subset <- merge(cna_full , panel_cna , all.x=TRUE , by="gene_symbol")
    cna_subset <- cna_subset[ !is.na(cna_subset$group) , ]
    idx_ok <- ifelse( cna_subset$CNA==cna_subset$exact_alteration | 
                     (cna_subset$exact_alteration=="" & cna_subset$CNA %in% 
                      c("amplification" , "deletion")), TRUE , FALSE)
    cna_subset <- cna_subset[ idx_ok , ] %>% 
                    unique %>% 
                    .[ , c("drug" , "group" 
                      , "gene_symbol" , "tumor_type" , "case_id")]
    if(nrow(cna_subset)>0){
      cna_subset$alteration_id <- paste0("cna_" , seq_len(nrow(cna_subset)))
    } else {
      cna_subset$alteration_id <- character(0)
    }
    return(cna_subset)
}