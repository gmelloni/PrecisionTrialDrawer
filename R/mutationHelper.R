# Given a vector of genes and a vector of tumor 
# type as resulted from showTumorType() or showCancerStudy()
# It returns a nested list composed by n elements, one for each tumor type
# Each element is composed by a mutation df and a 
# vector with all patient barcodes for that tumor type
.getMutations <- function(myGenes=myGenes
                            , tumor_type="all_tumors"
                            , block=NULL
                            , totalBlocks=NULL)
{
    mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
    allCanStudy <- cgdsr::getCancerStudies(mycgds)[,c(1,2)]
    allCanStudy$tumor_type <- vapply(strsplit(allCanStudy[,1] 
                                          , "_") , '[' , character(1) , 1)
    if(tumor_type[1]=="all_tumors") {
        chosenTumors <- allCanStudy[,1]
    } else {
        #FIND All cancer studies associated with the tumor_type specified
        chosenTumors <- allCanStudy[ allCanStudy[,'tumor_type'] %in% 
                                       tumor_type, 1]
        #In case we are not looking for tumor ID but cancer studies ID
        if(length(chosenTumors)==0){
            chosenTumors <- allCanStudy[ allCanStudy[,'cancer_study_id'] %in% 
                                           tumor_type, 1]
        } 
        #if still we have not found anything
        if (length(chosenTumors)==0){
          stop(paste("Could not find the tumor_type nor"
                     ,"cancer_study_id specified"))
        }
    }
    # Remove all old TCGA instances if the new pancanatlas is present
    newTCGA <- grep("_tcga_pan_can_atlas_2018$" , chosenTumors , value=TRUE)
    if(length(newTCGA)>0){
      newTCGAtum <- lapply( strsplit(newTCGA , "_") , '[' , 1) %>% 
        unlist %>% unique
      oldTCGA <- c( paste0(newTCGAtum , "_tcga") 
                    , paste0(newTCGAtum , "_tcga_pub") )
      chosenTumors <- setdiff(chosenTumors 
                            , oldTCGA)
    }
    out_double <- lapply(chosenTumors , function(i)
        {
            #for each Cancer Study, fetch the type 
            # of alteration (genetic profile) 
            #to be considered
            geneticProfile <- tryCatch({
                cgdsr::getGeneticProfiles(mycgds, i)
            } , error = function(e) {
                 url <- paste0(mycgds$.url
                      , "webservice.do?cmd=getGeneticProfiles&&cancer_study_id="
                      , i)
                 res <- httr::GET(url)
                 if(res$status_code!=200){
                     stop(paste("Problems with cBioPortal Connection at" , url))
                 }
                 df <- strsplit(httr::content(res) , "\n") %>% 
                   unlist %>% strsplit("\t")
                 df <- do.call("rbind" , df[-1]) %>% 
                   as.data.frame(stringsAsFactors=FALSE) %>% setNames(df[[1]])
                 return(df)
             })
            geneticProfile <- geneticProfile[ ,c(1,2)]
            # select mutations
            # NEW APPROACH, we retrieve more
            sel <- grepl("mutations$" , geneticProfile$genetic_profile_id)
            # in 99% of the cases mutations are present 
            # but in case, just return NULL
            if(!any(sel)){
              return( list( out=NULL , patients=NULL) )
            }
            geneticProfile <- geneticProfile[sel, 1]
            caseList <- tryCatch({
                cgdsr::getCaseLists(mycgds, i)
                } , error = function(e) {
                    url <- paste(mycgds$.url
                      , "webservice.do?cmd=getCaseLists&cancer_study_id="
                      , i, sep = "")
                    res <- httr::GET(url)
                    if(res$status_code!=200){
                     stop(paste("Problems with cBioPortal Connection at" 
                                , url))
                    }
                    df <- strsplit(httr::content(res) , "\n") %>% 
                      unlist %>% strsplit("\t")
                    df <- do.call("rbind" , df[-1]) %>% 
                      as.data.frame(stringsAsFactors=FALSE) %>% 
                      setNames(df[[1]])
                    return(df)
            })
            toLowCL <- tolower(caseList$case_list_name)
            #find if we have any sequenced tumor
            sel <- toLowCL=="sequenced tumors"
            # sometime 'Sequenced Tumors' is not the name of the case_list_name for mutations
            # If 'Sequenced Tumors' does not exist, try 'Samples with mutation data' first and 'All Tumors' next
            if(!any(sel)){
              sel <- toLowCL=="samples with mutation data" | toLowCL=="samples with mutation data."
            }
            if(!any(sel)){
              sel <- toLowCL=="all samples"
            }
            if(any(sel)) {
                if(is.null(block)){
                    message(paste("getting mutations from this cancer study:" 
                          , i ))
                } else {
                    message(paste("getting mutations from this cancer study:" 
                          , i , paste0("(" , block , "/" , totalBlocks , ")")))
                }
                caseListID <- caseList[sel, 1]
                error <- tryCatch(
                    muts <- cgdsr::getMutationData( mycgds 
                        , caseList=caseListID 
                        , geneticProfile=geneticProfile 
                        , genes=myGenes)
                    , error=function(e) {
                        message(paste("Impossible to retrive mutations from" 
                                      , i , "study"))
                        }
                    )
                if(!exists("muts")) {
                    muts <- NULL
                    patients <- NULL
                } else if(nrow(muts) == 0 ){
                  muts <- NULL
                  patients <- NULL
                } else {
                    if(ncol(muts)!=22) {
                        muts <- NULL
                        patients <- NULL
                    } else {
                        colsIlike <- c(
                                "entrez_gene_id"
                                ,"gene_symbol"
                                ,"case_id"
                                # ,"sequencing_center"
                                # ,"mutation_status"
                                ,"mutation_type"
                                # ,"validation_status"
                                ,"amino_acid_change"
                                # ,"functional_impact_score"
                                # ,"xvar_link"
                                # ,"xvar_link_pdb"
                                # ,"xvar_link_msa"
                                ,"chr"
                                ,"start_position"
                                # ,"end_position"
                                ,"reference_allele"
                                ,"variant_allele"
                                # ,"reference_read_count_tumor"
                                # ,"variant_read_count_tumor"
                                # ,"reference_read_count_normal"
                                # ,"variant_read_count_normal"
                                # ,"genetic_profile_id"
                                )
                        muts <- muts[ , colsIlike]
                        # Fix chromosome number
                        muts$chr <- muts$chr %>% as.character %>%
                          sub("chr" , "" , .) %>%
                          sub("23" , "X" , .) %>%
                          sub("24" , "Y" , .) %>%
                          sub("M" , "MT" , .)
                        muts <- muts[ muts$chr %in% c(1:22, "X", "Y", "MT") , ]
                        muts$genetic_profile_id <- i
                        patients <- strsplit(caseList[sel, 'case_ids'] 
                                             , split=" ")[[1]]
                    }
                }
            } else {
                muts <- NULL
                patients <- NULL
            }
            return( list( out=muts , patients=patients) )
        })
        names(out_double) <- chosenTumors
        return(out_double)
}

.subsetMutations <- function(panel , muts_full , rs_df)
{
    emptyDF <- data.frame(drug = character(0)
                          ,group = character(0)
                          ,gene_symbol = character(0)
                          ,tumor_type = character(0)
                          ,case_id = character(0)
                          ,genomic_position=character(0)
                          , stringsAsFactors = FALSE)
    if(any(panel$exact_alteration=="amino_acid_position")){
        panel_mut_1 <- panel[ panel$exact_alteration=="amino_acid_position" 
                            , c("drug", "gene_symbol" 
                                , "mutation_specification" ,"group")]
        mut_specs <- strsplit(panel_mut_1$mutation_specification , "-")
        panel_mut_1$start <- vapply(mut_specs , function(x) {
          as.numeric(x[1])} , numeric(1))
        panel_mut_1$end <- vapply(mut_specs , function(x) {
          as.numeric(x[2])} , numeric(1))
        panel_mut_1$mutation_specification <- NULL
        mut_data1 <- merge(muts_full , panel_mut_1 
                           , all.y=TRUE , by="gene_symbol")
        if(any(!is.na(mut_data1$case_id))){
            mut_data1 <- mut_data1[ !is.na(mut_data1$case_id) , ]
            # watch out from mutations without amino acid position
            mut_data1 <- mut_data1[!is.na(mut_data1$amino_position) 
                                   , , drop=FALSE]
            mut_data1 <- mut_data1[mut_data1$amino_position>=mut_data1$start & 
                mut_data1$amino_position<=mut_data1$end , , drop=FALSE] %>% 
                unique
        } else {
            mut_data1 <- emptyDF
        }        
        mut_data1 <- mut_data1[ , colnames(emptyDF)]
    } else {
        mut_data1 <- emptyDF
    }

    if(any(panel$exact_alteration=="amino_acid_variant")){
        panel_mut_2 <- panel[ panel$exact_alteration=="amino_acid_variant" 
                            , c("drug", "gene_symbol" 
                                , "mutation_specification" ,"group")]
        mut_data2 <- merge(muts_full , panel_mut_2 
                           , all.y=TRUE , by="gene_symbol")
        if(any(!is.na(mut_data2$case_id))){
            mut_data2 <- mut_data2[ !is.na(mut_data2$case_id) , ]
            mut_data2 <- mut_data2[mut_data2$amino_acid_change==
                                     mut_data2$mutation_specification 
                                   , , drop=FALSE] %>% unique
        } else {
            mut_data2 <- emptyDF
        }
        mut_data2 <- mut_data2[ , colnames(emptyDF)]
    } else {
        mut_data2 <- emptyDF
    }

    if(any(panel$exact_alteration=="mutation_type")){
        ktcga_types <- c(
        "In_Frame_Del", "In_Frame_Ins","Missense_Mutation"
        ,"Frame_Shift_Del", "Frame_Shift_Ins"
        , "Nonsense_Mutation",  "Splice_Site"
        , "Translation_Start_Site", "Nonstop_Mutation"
        , "Silent","3'UTR", "3'Flank"
        , "5'UTR", "5'Flank", "IGR1 ", "Intron", "RNA", "Targeted_Region")
        knotTransc <- c("3'UTR", "3'Flank", "5'UTR", "5'Flank"
                ,"IGR1", "IGR", "Intron", "RNA", "Targeted_Region")
        kmiss_type <- ktcga_types[c(1,2,3)]
        ktrunc_type <- ktcga_types[c(4,5,6)]
        panel_mut_3 <- panel[ panel$exact_alteration=="mutation_type" 
                            , c("drug", "gene_symbol" 
                                , "mutation_specification" ,"group")]
        mut_data3 <- merge(muts_full , panel_mut_3 , all.y=TRUE 
                           , by="gene_symbol")
        if(any(!is.na(mut_data3$case_id))){
            mut_data3 <- mut_data3[ !is.na(mut_data3$case_id) , ]
            mut_data3 <- mut_data3[(mut_data3$mutation_type %in% 
                    kmiss_type & mut_data3$mutation_specification=="missense") |
                    (mut_data3$mutation_type %in% ktrunc_type & 
                       mut_data3$mutation_specification=="truncating")
                    , , drop=FALSE] %>% unique
        } else {
            mut_data3 <- emptyDF    
        }
        mut_data3 <- mut_data3[ , colnames(emptyDF)]
    } else {
        mut_data3 <- emptyDF
    }

    if(any(panel$exact_alteration=="genomic_position" | 
           panel$exact_alteration=="dbSNP_rs")){
        panel_mut_4 <- panel[ panel$exact_alteration %in% 
                                c("genomic_position","dbSNP_rs")
                , c("drug", "gene_symbol" , "mutation_specification" ,"group")]
        if(!is.null(rs_df)){
            panel_mut_4$mutation_specification <- 
              .mapvalues(panel_mut_4$mutation_specification 
                                                    , from=rs_df$rs
                                                    , to=rs_df$genomic_range
                                                    ,warn_missing=FALSE)
        }
        mut_specs2 <- strsplit(panel_mut_4$mutation_specification , "-|:")
        panel_mut_4$start <- vapply(mut_specs2 , function(x) {
          as.numeric(x[2])} , numeric(1))
        panel_mut_4$end <- vapply(mut_specs2 , function(x) {
          as.numeric(x[3])} , numeric(1))
        panel_mut_4$mutation_specification <- NULL
        mut_data4 <- merge(muts_full , panel_mut_4 
                           , all.y=TRUE , by="gene_symbol")
        if(any(!is.na(mut_data4$case_id))){
            mut_data4 <- mut_data4[!is.na(mut_data4$case_id) , ]
            mut_data4$genomic_position_num <- 
              strsplit(mut_data4$genomic_position , ":") %>%
                                        vapply(. , '[' , character(1) , 2) %>%
                                        as.numeric
            mut_data4 <- mut_data4[ 
              mut_data4$genomic_position_num>=mut_data4$start & 
                mut_data4$genomic_position_num<=mut_data4$end , ] %>% unique
        } else {
            mut_data4 <- emptyDF
        }
        mut_data4 <- mut_data4[ , colnames(emptyDF)]
    } else {
        mut_data4 <- emptyDF
    }

    if(any(panel$exact_alteration=="genomic_variant")){
        panel_mut_5 <- panel[ panel$exact_alteration=="genomic_variant" 
                            , c("drug", "gene_symbol" 
                                , "mutation_specification" ,"group")]
        mut_data5 <- merge(muts_full , panel_mut_5 
                           , all.y=TRUE , by="gene_symbol")
        if(any(!is.na(mut_data5$case_id))){
            mut_data5 <- mut_data5[!is.na(mut_data5$case_id) , ]
            mut_data5 <- mut_data5[ mut_data5$genomic_variant==
                                mut_data5$mutation_specification, ] %>% unique
        } else {
            mut_data5 <- emptyDF
        }
        mut_data5 <- mut_data5[ , colnames(emptyDF)]
    } else {
        mut_data5 <- emptyDF
    }

    if(any(panel$exact_alteration=="")){
        panel_mut_6 <- panel[ panel$exact_alteration=="" 
                            , c("drug", "gene_symbol" 
                            , "mutation_specification" ,"group")]
        mut_data6 <- merge(muts_full , panel_mut_6 , all.y=TRUE 
                           , by="gene_symbol")
        if(any(!is.na(mut_data6$case_id))){
            mut_data6 <- mut_data6[!is.na(mut_data6$case_id) , ]
        } else {
            mut_data6 <- emptyDF
        }
        mut_data6 <- mut_data6[ , colnames(emptyDF)]
    } else {
        mut_data6 <- emptyDF
    }
    muts_subset <- do.call("rbind" 
          , list(mut_data1,mut_data2,mut_data3,mut_data4,mut_data5,mut_data6))
    muts_subset <- muts_subset[ , c("drug", "group", "gene_symbol", "tumor_type"
                                    , "case_id","genomic_position")] %>% unique
    if(nrow(muts_subset)>0){
      muts_subset$alteration_id <- paste0("mut_" , seq_len(nrow(muts_subset)))
    } else {
      muts_subset$alteration_id <- character(0)
    }
    muts_subset$genomic_position <- NULL
    return(muts_subset)
}
