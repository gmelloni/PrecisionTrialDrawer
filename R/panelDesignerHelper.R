################################################################################
# This function takes a list of RS id and fetched the genomic coordinates.
# 
# INPUT: 
#   List of RS ids
# OUTPUT:
#   df with rs id and genomic coordinates on hg19
.rsToGenomicPosition <- function(rs)
{
  #if rs == "" it would download all the RS in the database. 
  if(length(rs)==1 && rs==""){
    stop("Input rs list was empty - Function: .rsToGenomicPosition()")
  }
  # Do not look at the same rs twice
  rs <- unique(rs)
  # Query Biomart HG19
  snp_mart = biomaRt::useMart(biomart="ENSEMBL_MART_SNP"
                    , host="feb2014.archive.ensembl.org"
                    , dataset="hsapiens_snp")
  #run the query
  BM_rs = biomaRt::getBM(attributes=c( "chr_name"
              , "chrom_start"
              , "refsnp_id"
          )
          , filters="snp_filter"
          , values=rs
          , mart=snp_mart
      )
  # If we cannot find all the snps, this function should break
  if(!all(rs %in% BM_rs$refsnp_id)){
    stop(".rsToGenomicPosition() - Not all the RS ids where found in Biomart")
  }
  # If some SNPs are mapped on multiple positions or they are mapped on non canonical chr
  # We first delete all non canonical chromosomes, then check if they are uniquely mapped
  # If they are not, we choose the first result
  BM_rs_split <- split(BM_rs , BM_rs$refsnp_id)
  BM_rs_split_adj <- lapply(BM_rs_split , function(df) {
                if(nrow(df)==1){
                    return(df)
                } else {
                    df2 <- df[ as.character(df$chr_name) %in% c(1:22 , "X" , "Y" , "M" , "MT") , , drop=FALSE]
                    if(nrow(df2)==1){
                        return(df2)
                    } else if(nrow(df2)==0){
                        warning("The following rs ID is not mapped on canonical chromosome:\n"
                            , immediate.=TRUE)
                        print(df)
                        return(df[1 , , drop=FALSE])
                    } else {
                        warning("The following rs ID is mapped on multiple positions, the first one was chosen:\n"
                            , immediate.=TRUE)
                        print(df2)
                        return(df2[1 , , drop=FALSE])
                    }
                }
            })
  BM_rs <- do.call("rbind" , BM_rs_split_adj)
  BM_rs_out <- data.frame(
      rs=BM_rs$refsnp_id
      ,genomic_range=paste( 
        BM_rs$chr_name
        , paste(BM_rs$chrom_start , BM_rs$chrom_start , sep="-")
        , sep=":")
      ,stringsAsFactors=FALSE)
  # Order by the original order
  BM_rs_out <- BM_rs_out[ order(match(BM_rs_out$rs , rs)) , ]
  return(BM_rs_out)
}

################################################################################
# This function takes a list of gene-aminoacid positions and fetched the genomic coordinates.
# 
# INPUT: 
#   panel_aa
# OUTPUT:
#   df with genomic coordinates for each AA change or AA 
.fromAAtoGenomic <- function(panel_aa , hgnc_query , BPPARAM=bpparam("SerialParam") , myhost="www.ensembl.org") {
    if(any(panel_aa$exact_alteration=="amino_acid_variant")){
        aa <- panel_aa[panel_aa$exact_alteration=="amino_acid_variant" , "mutation_specification"]
        aa <- data.frame(aa=aa , aa_num=stringr::str_extract(aa,"\\d+") , stringsAsFactors=FALSE)
        aa$aa_num <- paste(aa$aa_num , aa$aa_num , sep="-")
        panel_aa$mutation_specification <- .mapvalues(panel_aa$mutation_specification , aa$aa , aa$aa_num , warn_missing=FALSE)
    }
    genes <- unique(panel_aa$gene_symbol)
    # uniprot_can <- readRDS(file.path(system.file(package="PrecisionTrialDesigner")
    #                                 , "extdata"
    #                                 , "uniProtCanonical.rds"))
    # uniprot_can <- LowMACAAnnotation::getMyUni()
    # uniprot_can$prot_len <- nchar(uniprot_can$AMINO_SEQ)
    # Here we use just the latest version of ensembl because the REST API for proteins works in gchr38
    # We have to manually revert the positions to hg19 with annotationhub
    ensembl=biomaRt::useMart(host=myhost , biomart="ENSEMBL_MART_ENSEMBL" , dataset="hsapiens_gene_ensembl")
    # Retrieve proteins in biomart
    # We only want genes which are listed in HGNC
    # We only want proteins that are in uniprot
    bm <- biomaRt::getBM(mart=ensembl, values = genes
             ,filters = "hgnc_symbol"
             ,attributes = c(
                 "peptide"
                 # ,"peptide_location"
                 #,"uniprot_swissprot"
                 ,"uniprotswissprot"
                 ,"ensembl_gene_id"
                 ,"ensembl_peptide_id"
                 ,"hgnc_symbol"
                 ,"ensembl_transcript_id"
                 # ,"cds_length"
                 # ,"chromosome_name"
                 )
            ) %>% 
            split(. , .$hgnc_symbol) %>%
            lapply(. , function(x) {
                    x$prot_len <- nchar(x$peptide) - 1
                    out <- x[ x$ensembl_gene_id %in% hgnc_query[ , 'ensembl_gene_id'] , ]
                    out <- out[out$uniprotswissprot!="" , ]
                    if(nrow(out)==0)
                        stop("you look for a gene with no uniprot entry")
                    else
                        return(out)
                    }) %>%
            do.call("rbind" , .)
    # DEPRECATED 04/20/2018
    # Removed dependency from LowMACAAnnotation
    # bm_uniprot <- merge(bm , uniprot_can 
    #                     , by.x=c("uniprotswissprot") 
    #                     , by.y=c("Entry") 
    #                     , all.x=TRUE) %>%
    #                 .[ .$prot_len.x==.$prot_len.y , ]
    # bm_uniprot$Edit <- sapply(1:nrow(bm_uniprot) 
    #                     , function(i) adist(bm_uniprot[i,"AMINO_SEQ"] 
    #                                         , bm_uniprot[i,"peptide"]) 
    #                     )
    bm_uniprot_split <- split(bm , bm$hgnc_symbol) %>%
                        lapply(. , function(x) {
                            if(nrow(x)==1){
                                return(x)
                            } else {
                                # return the longest protein
                                # TODO: adopt a thing similar to canonical transcript
                                out <- out[ out$prot_len == max(out$prot_len , na.rm=TRUE) , ]
                                # out <- x[ x$Edit==min(x$Edit) , ]
                                if(nrow(out)>1)
                                    return(out[1 , , drop=FALSE])
                                else
                                    return(out)
                            }
                            }) %>% sapply(. , function(x) x[ , "ensembl_peptide_id"])
    #hub <- AnnotationHub::AnnotationHub()
    #chain <- AnnotationHub::query(hub, 'hg38ToHg19')[[1]]
    hg38Rest <- "http://rest.ensembl.org/map/translation/"
    hg19Rest <- "http://grch37.rest.ensembl.org/map/translation/"
    ranges_bedstyle <- bplapply(genes , function(gene) {
                        gene_df <- panel_aa[ panel_aa$gene_symbol==gene, ]
                        df <- lapply(unique(gene_df$mutation_specification) , function(x) {
                            # browser()
                            query <- paste0(hg19Rest
                                           ,bm_uniprot_split[gene] 
                                           , "/"
                                           , sub("-" , ".." , x)
                                           ,"?content-type=application/json")
                            ########## httr could be a valid alternative jsonlite, but requires more operations
                            ########## Consider adopting it to remove dependency to jsonlite
                            # json_stream <- httr::GET(query) %>%
                            #                 httr::content("parsed")
                            # json_stream <- json_stream[["mappings"]]#[[1]] %>% as.data.frame(stringsAsFactors=FALSE)
                            # json_stream <- lapply(json_stream , unlist) %>%
                            #                         do.call("rbind" , .) %>%
                            #                         as.data.frame(stringsAsFactors=FALSE)
                            ###############################
                            mycon <- url(query , open="r")
                            json_stream <- suppressWarnings(
                                            suppressMessages(
                                                jsonlite::stream_in(con=mycon , verbose=FALSE)$mappings[[1]]
                                            ))
                            close(mycon)
                            #Sys.sleep(1)
                            return(json_stream)
                            }) %>% do.call("rbind" , .)
                        gr <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(paste0("chr" , unique(df$seq_region_name)) , nrow(df))
                                        ,ranges=IRanges::IRanges(start=df$start , end=df$end)
                                        ,strand=df$strand)
                        # res <- suppressMessages(rtracklayer::liftOver(gr, chain)) %>% unlist
                        res <- gr
                        res <- IRanges::reduce(res , min.gapwidth=1L)
                        } , BPPARAM=BPPARAM)
    names(ranges_bedstyle) <- genes
    final_df <- lapply(genes , function(x) {
                    df <- as.data.frame(ranges_bedstyle[[x]])
                    df$gene_symbol <- x
                    df$strand <- NULL
                    return(df)
        }) %>% do.call("rbind" , .)
    final_df$target_chr <- sub("chr" , "" , final_df$seqnames)
    colnames(final_df)[2] <- "target_start"
    colnames(final_df)[3] <- "target_end"
    colnames(final_df)[4] <- "target_width"
    final_df <- final_df[ ,c("target_chr" , "target_start" 
                                        , "target_end" , "target_width" 
                                        , "gene_symbol")]
    return(final_df)
}

.fromIntervalstoGenomic <- function(panel_gn) {
    # Convert rs number into genomic locations and substitute them permanently (ex. rs1234567 becomes 1:1000-1000)
    if(any(panel_gn$exact_alteration=="dbSNP_rs")){
        rs <- panel_gn[panel_gn$exact_alteration=="dbSNP_rs" , "mutation_specification"]
        rs_df <- .rsToGenomicPosition(rs)
        panel_gn$mutation_specification <- .mapvalues(panel_gn$mutation_specification , rs_df$rs , rs_df$genomic_range)
    } else {
        rs_df <- NULL
    }
    # Convert genomic variants in genomic locations and substitute them permanently (ex. 1:1000:A,C becomes 1:1000-1000)
    if(any(panel_gn$exact_alteration=="genomic_variant")){
        gv <- panel_gn[panel_gn$exact_alteration=="genomic_variant" , "mutation_specification"]
        gv <- data.frame(gv=gv 
                        , genomic_range=strsplit(gv , ":") %>% sapply(. , function(x) paste(x[1] , paste(x[2] , x[2] , sep="-") , sep=":"))
                        , stringsAsFactors=FALSE
                        )
        panel_gn$mutation_specification <- .mapvalues(panel_gn$mutation_specification , gv$gv , gv$genomic_range)
    }
    mut_spec_split <- strsplit(panel_gn$mutation_specification , ":|-") %>% do.call("rbind" , .)
    # class(mut_spec_split) <- "numeric"
    colnames(mut_spec_split) <- c("target_chr" , "target_start" , "target_end")
    # mut_spec_split <- cbind(mut_spec_split 
    #             , target_width=(as.numeric(mut_spec_split[,"target_end"]) - as.numeric(mut_spec_split[,"target_start"]) + 1))
    
    panel_gn <- cbind(panel_gn[ , "gene_symbol" , drop=FALSE] , mut_spec_split)
    panel_gn <- .changeFactor(panel_gn)
    # Collapse redundant regions
    panel_gn_split <- split(panel_gn , panel_gn$gene_symbol)
    panel_gn_split <- lapply(names(panel_gn_split) , function(x) {
                        # browser()
                        ir <- IRanges::IRanges(start=as.numeric(panel_gn_split[[x]]$target_start)
                                    , end=as.numeric(panel_gn_split[[x]]$target_end)
                                            )
                        red <- IRanges::reduce(ir , min.gapwidth=1L)
                        red <- as.data.frame(red)
                        red$gene_symbol <- x
                        red$target_chr <- panel_gn_split[[x]]$target_chr %>% unique
                        colnames(red) <- c("target_start" , "target_end" 
                                        , "target_width" , "gene_symbol" 
                                        , "target_chr")
                        return(red)
        }) %>% do.call("rbind" , .)
    panel_gn_split <- panel_gn_split[ , c("target_chr" , "target_start" 
                                        , "target_end" , "target_width" 
                                        , "gene_symbol")]
    return(list(panel_gn_split , rs_df))
}


################################################################################
# This is the core function of the panelDesigner Method
# 
# INPUT: 
#   panel
# OUTPUT:
#   panel updated with gene length
.calculateGenomicSpace <- function(panel , canonicalTranscript , BPPARAM=bpparam("SerialParam") , myhost="www.ensembl.org"){
  # genes <- unique(panel$gene_symbol)
  genes <- panel$gene_symbol %>% strsplit(. , "__") %>% unlist %>% unique
  # message("Retrieving gene transcripts ranges...")
  # Retrieve full exon length from ENSEMBL
  #########################################
  # Initial set up for biomart queries
  # ---------------------------------------
  # Why 2 marts?
  # measures are taken from archived version because it is still hg19
  # attributes like transcript tsl and transcript source can only be found in the most recent versions
  message("Connecting to ensembl biomart...")
  ensembl=biomaRt::useMart(host=myhost 
                  , biomart="ENSEMBL_MART_ENSEMBL" 
                  , dataset="hsapiens_gene_ensembl")
  ensemblold=biomaRt::useMart(host="feb2014.archive.ensembl.org" 
                     , biomart="ENSEMBL_MART_ENSEMBL" 
                     , dataset="hsapiens_gene_ensembl")
  dframe <- biomaRt::getBM(attributes=c("hgnc_symbol", "ensembl_transcript_id","ensembl_gene_id")
                  , filters=c("hgnc_symbol"), values=genes, mart=ensemblold)# %>%
              # .[.$ensembl_gene_id %in% hgnc_query[ , 'ensembl_gene_id'] , ]
  # Get transcript associated attributes
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
  #clean up attribute fetching only the tsl number
  dframe_att$transcript_tsl <- strsplit(dframe_att$transcript_tsl , " ") %>% 
    sapply(. , '[' , 1) %>%
    sub("^tsl" , "" , .)
  dframe_att$transcript_tsl <- suppressWarnings(as.numeric(dframe_att$transcript_tsl))
  # Substitute NA transcript tsl with the highest number of tsl
  tsl_max <- ifelse( is.na(max(dframe_att$transcript_tsl , na.rm=TRUE)) 
                     , 1 
                     , max(dframe_att$transcript_tsl , na.rm=TRUE))
  dframe_att$transcript_tsl[is.na(dframe_att$transcript_tsl)] <- tsl_max
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
  # merge the latest dataframe with the very first biomart query and discard all 
  # regions with CDS length as NA value
  dframe_merge <- merge(dframe , dframe2 , all.x=TRUE) # %>% .[!is.na(.$cds_length) , ]
  # split dataframe for each gene name. This df will provide us all possible 
  # transcripts for each gene
  dframe_merge <- split(dframe_merge , dframe_merge$hgnc_symbol)

  #########################################
  # ADJUST GENES WITH NO CDS (like pseudo genes or RNA genes)
  #----------------------------------------
  # Remove genes with no canonical chromosome (if possible)
  # Keep only the genes with a cds_length (if possible)
  # If no cds_length available, raise a warning and use the exon length instead
  chrs <- c(1:22 , "X" , "Y" , "MT" , "M")
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
                          lapply(. , function(transcript){
                            cds_length <- transcript$genomic_coding_end - transcript$genomic_coding_start
                            cds_length <- sum(cds_length , na.rm=TRUE)
                            return(data.frame(ensembl_transcript_id=unique(transcript$ensembl_transcript_id)
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
    chrs <- c(1:22 , "X" , "Y" , "MT" , "M")
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
      }
      # Select transcript with longest cds 
      chosenTransc <- unique(x[ x$cds_length==max(x$cds_length , na.rm=TRUE) , "ensembl_transcript_id"])
      # return it if the selection was successful and return ONLY 1 transcript
      if(length(chosenTransc)==1 & !is.na(chosenTransc[1])){
        return(x[x$ensembl_transcript_id==chosenTransc , ,drop=FALSE])
      } else {
        # select all the remaining transcripts
        longest <- x[x$ensembl_transcript_id %in% chosenTransc , ]
        # select transcripts from havana and in case is only one, return it
        havana <- longest[ longest$transcript_source %in% "ensembl_havana" , , drop=FALSE]
        if(length(unique(havana$ensembl_transcript_id))==1){
          return(havana)
        } else {
          #there are either more than 1 transcript from havana or 0.
          if(nrow(havana)!=0){
            # 2 or more havana transcripts:
            # choose the best tsl support among the havana transcripts
            tsl <- havana[ havana$transcript_tsl %in% min(havana$transcript_tsl , na.rm=TRUE) ,  , drop=FALSE]
            lastChance <- tsl[ tsl$ensembl_transcript_id==sort(unique(tsl$ensembl_transcript_id))[1] ,  , drop=FALSE]
            return(lastChance)
          } else {
            # 0 havana transcripts: 
            # choose the best tsl support among the longest transcripts
            tsl <- longest[ longest$transcript_tsl %in% min(longest$transcript_tsl , na.rm=TRUE) ,  , drop=FALSE]
            lastChance <- tsl[ tsl$ensembl_transcript_id==sort(unique(tsl$ensembl_transcript_id))[1] ,  , drop=FALSE]
            return(lastChance)
          }
        }
      }
    })
  }
  
  dframe_merge <- lapply(dframe_merge , function(x) 
          x[ , colnames(x) %notin% c("transcript_tsl" , "transcript_source")]
          ) %>% do.call("rbind" , .)
  # These genes will be considered in their full sequence
  full_sequence_genes <- panel[ panel$alteration %in% c("CNA" , "expression") | 
                              panel$exact_alteration %in% c("","mutation_type") 
                              , "gene_symbol"] %>% 
                          strsplit(. , "__") %>%
                          unlist %>%
                          unique
  missingDframe <- setdiff(full_sequence_genes , unique(dframe_merge$hgnc_symbol) )
  if(length(missingDframe)>0){
    warning("We could not find these genes in biomaRt:" %++% paste(missingDframe , collapse=", "))
  }
  # Genes annotated as aminoacid variants/amino acid intervals
  panel_aa <- panel[ panel$exact_alteration %in% c("amino_acid_position","amino_acid_variant") , , drop=FALSE]
  #panel_aa <- panel_aa[ panel_aa$gene_symbol %notin% full_sequence_genes , , drop=FALSE]
  if(nrow(panel_aa)!=0){
      message("Reverse mapping of amino acids to genomic...")
      # Make a REST query to HGNC to retrieve a unique couple gene_symbol:ENSEMBLid
      hgnc_query <- BiocParallel::bplapply(genes , function(gene) {
       query <- "http://rest.genenames.org/fetch/symbol/" %+% gene
       # tmp <- tempfile()
       # download(query , tmp , mode="wb" , quiet=TRUE)
       tmp <- httr::GET(query)
       tmp <- XML::xmlParse(tmp)
       # check if the query was good
       found <- XML::xmlToList(tmp)$result$.attr['numFound'] %>% as.numeric
       if(found>1)
           stop("You didn't provide a unique valid official gene symbol for: " %+% gene)
       if(found==0){
           # If I can't find the official gene symbol, I look at previous symbols
           # The reason for this is that Ensembl is not always up-to-date
           query <- "http://rest.genenames.org/fetch/prev_symbol/" %+% gene
           tmp <- httr::GET(query)
           tmp <- XML::xmlParse(tmp)
           # tmp <- tempfile()
           # download(query , tmp , mode="wb" , quiet=TRUE)
           found <- XML::xmlToList(tmp)$result$.attr['numFound'] %>% as.numeric
           if(found>1)
             stop("You didn't provide a unique valid official gene symbol for: " %+% gene)
           if(found==0)
             stop("This gene symbol is not an official hgnc symbol: " %+% gene)
           hgnc_df <- XML::xmlToList(tmp) %>% 
             .[['result']] %>% 
             .[['doc']] %>% 
             lapply(. , function(x) {
                   if("str" %in% names(x)) 
                     unname(c(x$str , x$.attrs))
                   else
                     unname(c(x$text , x$.attrs))
               }) %>%
             do.call("rbind" , .) %>%
             .[ .[,2] %in% c("hgnc_id" , "prev_symbol" , "name" , "entrez_id" , "ensembl_gene_id"), ] %>%
             t(.)
           hgnc_df_names <- as.vector(hgnc_df[2 , , drop=TRUE])
           hgnc_df <- as.vector(hgnc_df[1 , , drop=TRUE])
           names(hgnc_df) <- hgnc_df_names
           names(hgnc_df)[names(hgnc_df)=="prev_symbol"] <- "symbol"
           hgnc_df <- hgnc_df[c("hgnc_id" , "symbol" , "name" , "entrez_id" , "ensembl_gene_id")]
           return(hgnc_df)
       }
       hgnc_df <- XML::xmlToList(tmp) %>% 
                   .[['result']] %>% 
                   .[['doc']] %>% 
                   lapply(. , function(x) unname(c(x$text , x$.attrs))) %>%
                   do.call("rbind" , .) %>%
                   .[ .[,2] %in% c("hgnc_id" , "symbol" , "name" , "entrez_id" , "ensembl_gene_id"), ] %>%
                   t(.)
       hgnc_df_names <- as.vector(hgnc_df[2 , , drop=TRUE])
       hgnc_df <- as.vector(hgnc_df[1 , , drop=TRUE])
       names(hgnc_df) <- hgnc_df_names
       hgnc_df <- hgnc_df[c("hgnc_id" , "symbol" , "name" , "entrez_id" , "ensembl_gene_id")]
       return(hgnc_df)
      } , BPPARAM=BPPARAM) %>% do.call("rbind" , .)
      aa_intervals <- .fromAAtoGenomic(panel_aa, hgnc_query=hgnc_query , BPPARAM=BPPARAM , myhost=myhost)
  } else {
      aa_intervals <- NULL
  }
  # Genes annotated as genomic ranges , genomic alterations and rs numbers
  panel_gn <- panel[ panel$exact_alteration %in% c("genomic_position","genomic_variant","dbSNP_rs") , , drop=FALSE]
  #panel_gn <- panel_gn[ panel_gn$gene_symbol %notin% full_sequence_genes , , drop=FALSE]
  if(nrow(panel_gn)!=0){
      message("Unifying genomic annotations...")
      gn_intervals <- .fromIntervalstoGenomic( panel_gn=panel_gn )
  } else {
      gn_intervals <- list(NULL , NULL)
  }
  intervals <- rbind(aa_intervals , gn_intervals[[1]])
  # rs_df <- gn_intervals[[2]]
  return(list(GeneIntervals=dframe_merge 
              , TargetIntervals=intervals 
              , FullGenes=full_sequence_genes 
              # , dbSNP_rs=rs_df
              ))
}