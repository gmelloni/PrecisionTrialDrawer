# Given a vector of genes and a vector of tumor type
# as resulted from showTumorType() or showCancerStudy()
# It returns a nested list composed by n elements, one for each tumor type
# Each element is composed by a CNA df and a vector 
# with all patient barcodes for that tumor type
.getExpression <- function(myGenes=myGenes
                            ,tumor_type="all_tumors"
                            ,block=NULL)
{
  mycgds <- cBioPortalData::cBioPortal(
    hostname = "www.cbioportal.org",
    protocol = "https",
    api. = "/api/api-docs")
  allCanStudy <- cBioPortalData::getStudies(mycgds)[ , c("cancerTypeId" , "studyId")]
  allCanStudy <- unique(allCanStudy)
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
    geneticProfile <- geneticProfile[ grepl("MRNA_EXPRESSION" , geneticProfile$molecularAlterationType , ignore.case=TRUE) & 
                                         geneticProfile$datatype == "Z-SCORE" &
                                         grepl("mRNA" , geneticProfile$name) & grepl("diploid" , geneticProfile$name)
                                      # &
                                      #    grepl("RNA Seq V2 RSEM" , geneticProfile$name)
                                       # ( grepl("RNA Seq V2 RSEM" , geneticProfile$name) | grepl("\\(RPKM, log2 transformed\\)" , geneticProfile$name))
                                       , ]

    if(nrow(geneticProfile)>1){
      geneticProfile <- geneticProfile[ grepl("RNA Seq V2 RSEM" , geneticProfile$name) , ]
    }
    if(nrow(geneticProfile)==0){
      message(paste("geneticProfile" , i , "has no expression data"))
      return( list( out=NULL , patients=NULL , type=NULL) )
    }
    type <- if(any(grepl("array" , geneticProfile$molecularProfileId))) "array" else "rnaseq"
    # Fetch the patient list for the specified molecular profile and data type (mutations)
    caseList <- cBioPortalData::sampleLists(mycgds, i)
    # caseList <- caseList[ caseList$category == "all_cases_with_log2_cna_data" , ]
    caseList <- caseList[ caseList$category %in% c("all_cases_with_mrna_rnaseq_data","all_cases_with_mrna_array_data") , ]
    if(nrow(caseList)==0){
      message(paste("ProfileId" , i , "has no samples with expression"))
      return( list( out=NULL , patients=NULL , type=NULL) )
    }
    # Get the actual case names (e.g. TCGA-AR-A1AR)
    caseListId <- cBioPortalData::getSampleInfo(api = mycgds 
                                                , studyId = geneticProfile$studyId 
                                                , sampleListIds = caseList$sampleListId)$sampleId
    tryCatch(
      express <- cBioPortalData::molecularData(
        api = mycgds
        , molecularProfileIds = geneticProfile$molecularProfileId
        , entrezGeneIds = as.numeric(names(myGenes))
        , sampleIds = caseListId
      )[[1]]
      , error=function(e) message(paste("Impossible to retrive copynumber from" , i , "study or no copynumber are present on the selected genes"))
    )
    if(!exists("express")) {
      express <- NULL
      patients <- NULL
      type <- NULL
    } else if(nrow(express) == 0 ){
      express <- NULL
      patients <- NULL
      type <- NULL
    } else {
      patients <- unique(ifelse(grepl("^TCGA" , caseListId) 
                         , unlist(lapply(strsplit(caseListId , "-") , function(x) paste(x[seq_len(3)] , collapse="-")))
                         , caseListId))
      express$genetic_profile_id <- i
      express$tumor_type <- chosenTumors[i]
      express$case_id <- ifelse(grepl("^TCGA" , express$sampleId) 
                           , unlist(lapply(strsplit(express$sampleId , "-") , function(x) paste(x[seq_len(3)] , collapse="-")))
                           , express$sampleId)
      express$gene_symbol <- myGenes[ as.character(express$entrezGeneId) ]
      express <- express[ , c("case_id" , "entrezGeneId" , "gene_symbol" , "value" , "genetic_profile_id" , "tumor_type")]
    }
    return( list( out=express , patients=patients , type=type) )
  # } , .inform = TRUE)
  })
  names(out_double) <- names(chosenTumors)
  return(out_double)
}

.subsetExpression <- function(panel , expr_full)
{
    panel_expr <- panel[ panel$alteration=="expression" , ]
    expr_subset <- merge(expr_full , panel_expr 
                         , all.x=TRUE , by="gene_symbol")
    expr_subset <- expr_subset[ !is.na(expr_subset$group) , ]
    idx_ok <- ifelse( expr_subset$expression==expr_subset$exact_alteration | 
                     (expr_subset$exact_alteration=="" & 
                        expr_subset$expression %in% c("up" , "down"))
                     , TRUE , FALSE)
    expr_subset <- expr_subset[ idx_ok , ] %>% 
                    unique %>% 
                    .[ , c("drug" , "group" 
                           , "gene_symbol" , "tumor_type" , "case_id")]
    if(nrow(expr_subset)>0){
      expr_subset$alteration_id <- paste0("expr_" , seq_len(nrow(expr_subset)))
    } else {
      expr_subset$alteration_id <- character(0)
    }
    return(expr_subset)
}
