# Return a list with:
    # Gene length divided by exons and cds regions
    # All single alteration in genomic ranges style (both from aminoacid ranges/genomic ranges)
    # Vector of genes in the panel taken as full sequence
    # bed style dataframe ready for submission to your preferred sequencing company 
setGeneric('panelDesigner', function(object 
                            , alterationType = c("copynumber", "expression", "mutations", "fusions") 
                            , padding_length=100 , merge_window=50 , utr=FALSE 
                            , canonicalTranscript=TRUE , BPPARAM=bpparam("SerialParam")
                            , myhost="www.ensembl.org")  {
    standardGeneric('panelDesigner')
    })
setMethod('panelDesigner', 'CancerPanel', function(object 
                            , alterationType = c("copynumber", "expression", "mutations", "fusions") 
                            , padding_length=100 , merge_window=50 , utr=FALSE 
                            , canonicalTranscript=TRUE , BPPARAM=bpparam("SerialParam")
                            , myhost="www.ensembl.org") {
    panel <- object@arguments$panel
    alterationType <- .mapvalues(from=c("copynumber", "expression", "mutations", "fusions") 
                                , to=c("CNA" , "expression" , "SNV" , "fusion")
                                , alterationType
                                , warn_missing=FALSE)
    panel <- panel[ panel$alteration %in% alterationType , ]
    panel_design <- .calculateGenomicSpace(panel , canonicalTranscript=canonicalTranscript , BPPARAM = BPPARAM , myhost=myhost)
    geneIntervals <- panel_design$GeneIntervals
    padding_range <- padding_length*2 + 1
    if(utr){
        full_gene_ranges <- geneIntervals[geneIntervals$hgnc_symbol %in% panel_design$FullGenes 
                                    , c("chromosome_name" , "exon_chrom_start" , "exon_chrom_end")
                                    ]
    } else {
        full_gene_ranges <- geneIntervals[geneIntervals$hgnc_symbol %in% panel_design$FullGenes 
                                    , c("chromosome_name" , "genomic_coding_start" , "genomic_coding_end")
                                    ]
    }
    full_gene_ranges <- full_gene_ranges[ !is.na(full_gene_ranges[ , 2]) , ]
    colnames(full_gene_ranges) <- c("chr" , "start" , "end")
    if(!is.null(panel_design$TargetIntervals)){
        target_intervals <- panel_design$TargetIntervals[ , c(1,2,3)]
        colnames(target_intervals) <- c("chr" , "start" , "end")
        rownames(target_intervals) <- paste(panel_design$TargetIntervals$gene_symbol  
                                            , "target" 
                                            , rownames(target_intervals) , sep=".")
    } else {
        target_intervals <- NULL
    }
    toBeRanged <- rbind(full_gene_ranges , target_intervals)
    toBeRanged$chr <- paste0("chr" , toBeRanged$chr)
    width <- toBeRanged$end - toBeRanged$start
    toBeRanged$start <- vapply(seq_len(nrow(toBeRanged)) , function(i) {
                                    if(width[i]>=padding_range) 
                                        toBeRanged$start[i]
                                    else
                                        toBeRanged$start[i] - round((padding_range - width[i])/2)
                                    } , numeric(1))
    toBeRanged$end <- vapply(seq_len(nrow(toBeRanged)) , function(i) {
                                    if(width[i]>=padding_range)
                                        toBeRanged$end[i]
                                    else
                                        toBeRanged$end[i] + round((padding_range - width[i])/2)
                                    } , numeric(1))
    toBeRanged$Annotation <- rownames(toBeRanged) %>% 
                                strsplit(. , "\\.") %>%
                                vapply(. , '[' , character(1) , 1) %>%
                                unlist
    toBeRanged2 <- GenomicRanges::makeGRangesFromDataFrame(toBeRanged , keep.extra.columns = TRUE)
    toBeRanged3 <- IRanges::reduce(toBeRanged2 , min.gapwidth=merge_window)
    # Now we want to keep the annotation of toBeRanged2 in toBeRanged3
        # Find the overlaps between the original dataset and the reduce/bedmerged one
    IDX <- GenomicRanges::findOverlaps(toBeRanged2, toBeRanged3) %>% S4Vectors::subjectHits(.)
        # Aggregate all the rows of the bedmerged dataset with their gene symbol
    agg_annotation <- data.frame(rowToAgg=IDX , annotation=toBeRanged$Annotation , stringsAsFactors=FALSE) %>%
                        aggregate(annotation ~ rowToAgg , data=. , FUN=function(x) paste(unique(x) , collapse="|")) %>%
                        .[ order(.$rowToAgg) , ]
    agg_annotation$rowToAgg <- as.character(agg_annotation$rowToAgg)
        # Merge the annotation with the reduced/bedmerged dataframe
    toBeRanged_final <- as.data.frame(toBeRanged3)[ , c(1,2,3)]
    toBeRanged_final <- .mergeOrder(toBeRanged_final , agg_annotation , by.x="row.names" , by.y="rowToAgg" , all.x=TRUE)
    toBeRanged_final$Row.names <- NULL
        # Put the annotation
    colnames(toBeRanged_final) <- c("chr" , "start" , "end" , "annotation")
    toBeRanged_final$chr <- as.character(toBeRanged_final$chr)
        # Bed files are 0 based!
    toBeRanged_final$start <- toBeRanged_final$start - 1L
        # Add rs ID annotations
    dfRS <- object@arguments$dbSNP_rs
    if(nrow(dfRS)==1 & dfRS$rs[1]==""){
      toBeRanged_final_rs <- toBeRanged_final
    } else {
      splitdbSNP_rs <- strsplit(dfRS$genomic_range , ":|-")
      dfRS$chr <- paste("chr" , vapply(splitdbSNP_rs , '[' , character(1) , 1) , sep="")
      dfRS$startRS <- as.integer(vapply(splitdbSNP_rs , '[' , character(1) , 2))
      rs_region <- vapply(seq_len(nrow(dfRS)) , function(x) {
                        dfRS_row <- dfRS[x , , drop=FALSE]
                        toBeRanged_final_red <- toBeRanged_final[ toBeRanged_final$chr==dfRS_row$chr , , drop=FALSE]
                        out <- "none"
                        for(row in rownames(toBeRanged_final_red)){
                            if(dfRS_row$startRS>toBeRanged_final_red[row , "start"]){
                                if(dfRS_row$startRS<=toBeRanged_final_red[row , "end"]){
                                    out <- row
                                    break
                                }
                            }
                        }
                        return(out)
                    } , character(1))
      dfRS$bedRow <- rs_region
      dfRS_agg <- aggregate(rs ~ bedRow , dfRS , FUN=function(x) paste(x , collapse=";"))
      toBeRanged_final_rs <- merge(toBeRanged_final , dfRS_agg , by.x="row.names" , by.y="bedRow" , all.x=TRUE)
      toBeRanged_final_rs$Row.names <- NULL
      toBeRanged_final_rs$annotation <- ifelse(is.na(toBeRanged_final_rs$rs) 
                                            , toBeRanged_final_rs$annotation 
                                            , paste(toBeRanged_final_rs$annotation , toBeRanged_final_rs$rs , sep="|"))
      toBeRanged_final_rs$start <- as.integer(toBeRanged_final_rs$start)
      toBeRanged_final_rs$rs <- NULL
    }
    panel_design <- c(panel_design , list(BedStylePanel=toBeRanged_final_rs))
    return(panel_design)
})
