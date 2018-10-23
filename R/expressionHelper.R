# Given a vector of genes and a vector of tumor type as resulted from showTumorType() or showCancerStudy()
# It returns a nested list composed by n elements, one for each tumor type
# Each element is composed by a CNA df and a vector with all patient barcodes for that tumor type
.getExpression <- function(myGenes=myGenes
                            ,tumor_type="all_tumors"
                            ,block=NULL)
{
    mycgds <- cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
    all_cancer_studies <- cgdsr::getCancerStudies(mycgds)[,c(1,2)]
    all_cancer_studies$tumor_type <- vapply(strsplit(all_cancer_studies[,1] , "_") , '[' , character(1) , 1)
    if(tumor_type[1]=="all_tumors") {
        chosenTumors <- all_cancer_studies[,1]
    } else {
        # chosenTumors <- all_cancer_studies[grepl( paste(tumor_type , collapse="|") 
        #                                     , all_cancer_studies[,1] 
        #                                     , ignore.case=TRUE) , 1]
        chosenTumors <- all_cancer_studies[ all_cancer_studies[,'tumor_type'] %in% tumor_type, 1]
        if(length(chosenTumors)==0){
            chosenTumors <- all_cancer_studies[ all_cancer_studies[,'cancer_study_id'] %in% tumor_type, 1]
        }
    }
    out_double <- lapply(chosenTumors , function(i)
        {
            geneticProfile <- tryCatch({
                cgdsr::getGeneticProfiles(mycgds, i)
            } , error = function(e) {
                 url <- paste(mycgds$.url, "webservice.do?cmd=getGeneticProfiles&&cancer_study_id=", i, sep = "")
                 res <- httr::GET(url)
                 if(res$status_code!=200){
                     stop(paste("Problems with cBioPortal Connection at" , url))
                 }
                 df <- strsplit(httr::content(res) , "\n") %>% unlist %>% strsplit("\t")
                 df <- do.call("rbind" , df[-1]) %>% as.data.frame(stringsAsFactors=FALSE) %>% setNames(df[[1]])
                 message("\n")
                 return(df)
                #df <- read.table(url, skip = 0, header = TRUE, as.is = TRUE, sep = "\t", quote = "")
             })
            geneticProfile <- geneticProfile[ ,c(1,2)]
            # First choice RNAseq data, if not available use microarray
            if(any(grepl("_rna_seq_v2_mrna_median_Zscores$" , geneticProfile$genetic_profile_id))) {
                geneticProfile <- grep("_rna_seq_v2_mrna_median_Zscores$" , geneticProfile$genetic_profile_id 
                                    , value=TRUE , ignore.case=TRUE)
                caseList <- tryCatch({
                    cgdsr::getCaseLists(mycgds, i)
                    } , error = function(e) {
                        url <- paste(mycgds$.url, "webservice.do?cmd=getCaseLists&cancer_study_id=", i, sep = "")
                        res <- httr::GET(url)
                        if(res$status_code!=200){
                         stop(paste("Problems with cBioPortal Connection at" , url))
                        }
                        df <- strsplit(httr::content(res) , "\n") %>% unlist %>% strsplit("\t")
                        df <- do.call("rbind" , df[-1]) %>% as.data.frame(stringsAsFactors=FALSE) %>% setNames(df[[1]])
                        return(df)
                        #df <- read.table(url, skip = 0, header = TRUE, as.is = TRUE, sep = "\t", quote = "")
                    }
                    )
                sel <- caseList$case_list_name=="Tumor Samples with mRNA data (RNA Seq V2)"
                type <- "RNAseq"
            # TO DO:
                # Microarray seemed legit at the beginning, but the output is not actually the same as RNAseq
                # For now, we only rely on RNAseq
            # } else if(any(grepl("_mrna_median_Zscores$" , geneticProfile$genetic_profile_id))) {
            #     geneticProfile <- grep("_mrna_median_Zscores$" , geneticProfile$genetic_profile_id 
            #                         , value=TRUE , ignore.case=TRUE)
            #     caseList <- getCaseLists(mycgds, i)
            #     sel <- caseList$case_list_name=="Tumors with mRNA data (Agilent microarray)"
            #     type <- "Agilent Microarray"
            } else {
                sel <- FALSE
            }
            if(any(sel)) {
                if(is.null(block)){
                    message(paste("getting Expression from this cancer study:" , i ))
                } else {
                    message(paste("getting Expression from this cancer study:" , i, paste0("(" , block , ")")))
                }
                caseListID <- caseList[sel, 1]
                error <- tryCatch(
                    cna <- cgdsr::getProfileData( mycgds 
                        , caseList=caseListID 
                        , geneticProfile=geneticProfile 
                        , genes=myGenes)
                    , error=function(e) {
                        message(paste("Impossible to retrive Expression from" , i , "study"))
                        return(e)
                    }
                    )
                # if(error$message=="object 'myGenes' not found") stop("No gene data provided")
                if(!exists("cna")) {
                    cna <- NULL
                    patients <- NULL
                    type <- NULL
                } else {
                    patients <- strsplit(caseList[sel, 'case_ids'] , split=" ")[[1]]
                }
            } else {
                # message(paste("No CNA data from" , i , "study"))
                cna <- NULL
                patients <- NULL
                type <- NULL
            }
            return( list( out=cna , patients=patients , type=type) )
        })
        names(out_double) <- chosenTumors
        return(out_double)
}

.subsetExpression <- function(panel , expr_full)
{
    panel_expr <- panel[ panel$alteration=="expression" , ]
    expr_subset <- merge(expr_full , panel_expr , all.x=TRUE , by="gene_symbol")
    expr_subset <- expr_subset[ !is.na(expr_subset$group) , ] # there should be no NA, but who knowns
    idx_ok <- ifelse( expr_subset$expression==expr_subset$exact_alteration | 
                     (expr_subset$exact_alteration=="" & expr_subset$expression %in% c("up" , "down")), TRUE , FALSE)
    expr_subset <- expr_subset[ idx_ok , ] %>% 
                    unique %>% 
                    .[ , c("drug" , "group" , "gene_symbol" , "tumor_type" , "case_id")]
    if(nrow(expr_subset)>0){
      expr_subset$alteration_id <- paste0("expr_" , seq_len(nrow(expr_subset)))
    } else {
      expr_subset$alteration_id <- character(0)
    }
    return(expr_subset)
}