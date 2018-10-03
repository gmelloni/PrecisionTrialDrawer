###############################################################################
# This function gets the length of all the transcripts
# INPUT:
#   gene names list
# Option:
#   - canonicalTranscript
#    If TRUE, than it goes through a selection to find the canonical transcript
#    If FALSE, than it merges the possible isoform of the transcripts
# OUTPUT:
# a series of genes and transcripts length 
.annotateGeneLength <- function(genes 
                                , canonicalTranscript=TRUE 
                                , myhost="www.ensembl.org") {
  
  ########################################
  #Initial safe check
  if(length(genes)==0){
  stop("Submitted gene list is empty")
  }
  if(any(genes=="")){
  stop("Submitted gene list contains empty gene names")
  }    
  if(any(genes==" ")){
  stop("Gene list contains empty strings instead of gene names: \" \" ")
  }
  
  
  #########################################
  # Initial set up for biomart queries
  # ---------------------------------------
  # Why 2 marts?
  # measures are taken from archived version because it is still hg19
  # attributes like transcript tsl and transcript source 
  # can only be found in the most recent versions
  message("Connecting to ensembl biomart...")
  ensembl=biomaRt::useMart(host=myhost
                           , biomart="ENSEMBL_MART_ENSEMBL" 
                           , dataset="hsapiens_gene_ensembl")
  ensemblold=biomaRt::useMart(host="feb2014.archive.ensembl.org" 
                              , biomart="ENSEMBL_MART_ENSEMBL" 
                              , dataset="hsapiens_gene_ensembl")
  
  # ---------------------------------------
  # Run first Biomart Query on hgnc symbol to fetch info
  # ---------------------------------------
  # Due to the Biomart design, attributes from multiple attribute pages are not 
  # allow. In order to walkaround the problem, we have therefore performed 
  # separete queries, and marged them back together. This will be the main 
  # biomart dataset of refence upon which to run further queries. 
  dframe <- biomaRt::getBM(attributes=c("hgnc_symbol"
                                        , "ensembl_transcript_id"
                                        ,"ensembl_gene_id")
                           , filters=c("hgnc_symbol")
                           , values=genes
                           , mart=ensemblold)
  #Check for missing genes symbols in the query results
  if(any(genes %notin% dframe$hgnc_symbol)){
    genesNotFound <- setdiff(unique(genes) , unique(dframe$hgnc_symbol))
    errorMex <- "The following gene symbols were not found:" %++% 
      paste(genesNotFound , collapse=", ")
    stop(errorMex)
  }
  
  # ---------------------------------------
  # Get transcript associated attributes from previous query using transcriptID
  # ---------------------------------------
  dframe_att <- biomaRt::getBM(attributes=c(
    "ensembl_gene_id" 
    ,"ensembl_transcript_id"
    ,"transcript_tsl"
    ,"transcript_source"
    #,"transcript_status"
    #,"transcript_version"
    #,"hgnc_transcript_name"
    #,"uniprot_genename_transcript_name"
  )
  , filters=c("ensembl_transcript_id")
  , values=dframe[["ensembl_transcript_id"]]
  , mart=ensembl)
  # ---------------------------------------
  # Get transcript lengths running a new Biomart query (on the archive mart) 
  # starting from the very first BM query
  # ---------------------------------------
  dframe_len <- biomaRt::getBM(attributes=c(
    "cds_length"
    ,"chromosome_name"
    ,"genomic_coding_start"
    ,"genomic_coding_end"
    ,"exon_chrom_start"
    ,"exon_chrom_end"
    ,"ensembl_gene_id" 
    ,"ensembl_transcript_id"
  )
  , filters=c("ensembl_transcript_id")
  , values=dframe[["ensembl_transcript_id"]]
  , mart=ensemblold)
  
  #########################################
  # MERGING BACK QUERIES RESULTS TO ONE SINGLE DATAFRAME
  # ---------------------------------------
  # merge the results from the queries
  # the priority is given to hg19. If there is no attribute it is not important
  dframe2 <- merge(dframe_len , dframe_att , all.x=TRUE)
  # merge the latest dataframe with the 
  # very first biomart query and discard all 
  # regions with CDS length as NA value
  dframe_merge <- merge(dframe , dframe2 , all.x=TRUE)
  # split dataframe for each gene name. This df will provide us all possible 
  # transcripts for each gene
  dframe_merge <- split(dframe_merge , dframe_merge$hgnc_symbol)
  
  #########################################
  # ADJUST TSL NUMBER TO REMOVE NAs
  #----------------------------------------  
  
  #clean up attribute fetching only the tsl number
  dframe_merge <- lapply(dframe_merge , function(x) {
    if(all( x$transcript_tsl %in% c("tslNA" , NA)) ){
      x$transcript_tsl <- 1
      return(x)
    }
    x$transcript_tsl <- strsplit(x$transcript_tsl , " ") %>% 
      vapply(. , '[' , character(1) , 1 ) %>%
      sub("^tsl" , "" , .)
    x$transcript_tsl <- suppressWarnings(as.numeric(x$transcript_tsl))
    # Substitute NA transcript tsl with the highest number of tsl
    tsl_max <- ifelse( is.na(max(x$transcript_tsl , na.rm=TRUE)) 
                       , 1 
                       , max(x$transcript_tsl , na.rm=TRUE))
    x$transcript_tsl[is.na(x$transcript_tsl)] <- tsl_max
    return(x)
  })  
  
  #########################################
  # ADJUST GENES WITH NO CDS (like pseudo genes or RNA genes)
  #----------------------------------------
  # Remove genes with no canonical chromosome (if possible)
  # Keep only the genes with a cds_length (if possible)
  # If no cds_length available, raise a warning and use the exon length instead
  chrs <- c(seq_len(22) , "X" , "Y" , "MT" , "M")
  dframe_merge <- lapply(dframe_merge , function(x) {
    # Check that the gene is mapped to a chromosome if possible
    x_chr <- x[ as.character(x$chromosome_name) %in% chrs , ]
    if(nrow(x_chr)>0){
      x <- x_chr
    }
    # Check to make sure the gene in consideration is coding. If not skip it
    if(nrow(x[ !is.na(x$cds_length) , ])>0){
      # get coding lenghts (discarding those regions non coding)
      x <- x[ !is.na(x$cds_length) , , drop=FALSE]
      return(x)
    }
    if(nrow(x[ !is.na(x$cds_length) , ])==0){
      warning(paste("The following gene has no coding regions:" 
                    , unique(x$hgnc_symbol) 
                    , ". Only full exons length will be used"))
      x$genomic_coding_start <- x$exon_chrom_start
      x$genomic_coding_end <- x$exon_chrom_end
      fakecds_length <- split(x , x$ensembl_transcript_id) %>%
        lapply(. , function(trans){
          cds_length <- trans$genomic_coding_end - trans$genomic_coding_start
          cds_length <- sum(cds_length , na.rm=TRUE)
          return(data.frame(
                  ensembl_transcript_id=unique(trans$ensembl_transcript_id)
                  ,fakecds_length=cds_length
                  ,stringsAsFactors=FALSE))
        }) %>% do.call("rbind" , .)
      x <- merge(x , fakecds_length , all.x=TRUE)
      x$cds_length <- x$fakecds_length
      x$fakecds_length <- NULL
    }
    return(x)
  })
  #########################################
  # TRANSCRIPT OPTION1: SELECT CANONICAL
  # ---------------------------------------
  # For each gene we need to choose 1 transcript (the canonical) using the 
  #following decisional schema (in the shown order).
  # 1) Gene, if possible, should be mapped to a proper chromosome 1:22,X,Y,M,MT
  # 2) Gene with the longest coding length
  # 3) Gene in ensembl_havana merged db
  # 4) Gene with the highest TSL (transcript support level)
  # 5) The ENSTid with lower number
  if(canonicalTranscript){
    chrs <- c(seq_len(22) , "X" , "Y" , "MT" , "M")
    dframe_merge <- lapply(dframe_merge , function(x) {
      # Check that the gene is mapped to a chromosome if possible
      # for(x in dframe_merge){
      x_chr <- x[ as.character(x$chromosome_name) %in% chrs , ]
      if(nrow(x_chr)>0){
        x <- x_chr
      }
      # Check to make sure the gene in consideration is coding. If not skip it
      if(nrow(x[ !is.na(x$cds_length) , ])>0){
        # get coding lenghts (discarding those regions non coding)
        x <- x[ !is.na(x$cds_length) , , drop=FALSE]
      }
      # Select transcript with longest cds 
      chosenTransc <- unique(x[ x$cds_length==max(x$cds_length , na.rm=TRUE) 
                                , "ensembl_transcript_id"])
      # return it if the selection was successful and return ONLY 1 transcript
      if(length(chosenTransc)==1 & !is.na(chosenTransc[1])){
        return(x[x$ensembl_transcript_id==chosenTransc , ,drop=FALSE])
      } else {
        # select all the remaining transcripts
        longest <- x[x$ensembl_transcript_id %in% chosenTransc , ]
        # select transcripts from havana and in case is only one, return it
        havana <- longest[ longest$transcript_source %in% 
                             c("ensembl_havana" , "havana") , , drop=FALSE]
        if(length(unique(havana$ensembl_transcript_id))==1){
          return(havana)
        } else {
          #there are either more than 1 transcript from havana or 0.
          if(nrow(havana)!=0){
            # 2 or more havana transcripts:
            # choose the best tsl support among the havana transcripts
            tsl <- havana[ havana$transcript_tsl %in% 
                             min(havana$transcript_tsl , na.rm=TRUE) 
                           ,  , drop=FALSE]
            lastChance <- tsl[ 
              tsl$ensembl_transcript_id==
                sort(unique(tsl$ensembl_transcript_id))[1] 
              ,  , drop=FALSE]
            return(lastChance)
          } else {
            # 0 havana transcripts: 
            # choose the best tsl support among the longest transcripts
            tsl <- longest[ longest$transcript_tsl %in% 
                              min(longest$transcript_tsl , na.rm=TRUE) 
                            ,  , drop=FALSE]
            lastChance <- tsl[ tsl$ensembl_transcript_id==
                                 sort(unique(tsl$ensembl_transcript_id))[1] 
                               ,  , drop=FALSE]
            return(lastChance)
          }
        }
      }
    })
  }
  #########################################
  # TRANSCRIPT OPTION2: MERGE all TRANSCRIPTs
  # ---------------------------------------
  # In case we don't select one transcript
  # we have to collapse all exons of all transcripts
  codingIR <- lapply(dframe_merge 
                     , function(df) {
                df <- df[ !is.na(df$genomic_coding_start) , ] #rm NAs
                ir <- IRanges(start=df$genomic_coding_start 
                              , end=df$genomic_coding_end)
                ir <- reduce(ir , min.gapwidth=1L) #merge transcripts
                out <- data.frame(gene_symbol=unique(df$hgnc_symbol)
                    , cds_len=sum(ir@width)
                    , stringsAsFactors=FALSE)
  }) %>% do.call("rbind" , .)
  # for the UTRs
  codingUTRIR <- lapply(dframe_merge 
                        , function(df) {
                          df <- df[ !is.na(df$exon_chrom_start) , ]
                          ir <- IRanges(start=df$exon_chrom_start 
                                  , end=df$exon_chrom_end) #Select UTS too here
                          ir <- reduce(ir , min.gapwidth=1L)
                          out <- data.frame(gene_symbol=unique(df$hgnc_symbol)
                                            , cds_and_utr_len=sum(ir@width)
                                            , stringsAsFactors=FALSE)
                        }
  ) %>% do.call("rbind" , .)
  #merge the two dafatrame into one table
  final_df <- merge(codingIR , codingUTRIR , all.x=TRUE) %>% unique
  return(final_df)
}
