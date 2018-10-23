#----------------------------------
# CHECK BED FORMAT SPECIFICATIONS
#

.bedCheck <- function(bed){
    # bed must be a dataframe, char num num
    if( !is.data.frame(bed) ){
        stop("bed must be a data.frame")
    }
    # bed must have 3 columns
    if(ncol(bed)!=3)
        stop("bed data.frame must have 3 columns")
    # chr num num format
    # if(!is.character(bed[,1]))
    #     stop("chr column must be character")
    if(!is.numeric(bed[,2]))
    if(!is.numeric(bed[,3]))
        stop("end column must be numeric")
    if(! all(grepl("^chr" , bed[ , 1])))
        stop("first column must be in the format chr1-chr22,chrX,chrY,chrM or chrMT")
    # chr must be chrX format
    if( ! all(bed[,1] %in% paste0("chr" , c(seq_len(22) , "X" , "M" , "MT" , "Y"))) )
        stop("chr name not in chr1-chr22,chrX,chrY,chrM or chrMT")
    # check start end
    if( !all( bed[,2] < bed[,3]) )
        stop("some end are greater than start")
    # Give name to the column and change a few things
    colnames(bed) <- c("chr" , "start" , "end")
    if(any(bed$chr == "chrMT")){
        bed[ bed$chr == "chrMT" , "chr"] <- "chrM"
    }
    bed$chr <- gsub("chr" , "" , bed$chr)
    return(bed)
}

setGeneric('filterMutations'
           , function(object , filtered = NULL 
                  , bed = NULL , mode = c("exclude" , "keep") , tumor_type=NULL) {
    standardGeneric('filterMutations')
    })
setMethod('filterMutations', 'CancerPanel'
          , function(object , filtered 
                  , bed=NULL , mode = c("exclude" , "keep") , tumor_type=NULL)
{
    # browser()
    if( (is.null(filtered) && is.null(bed)) || (!is.null(filtered) && !is.null(bed)) ) {
        stop("Set filtered or bed parameters, not both")
    }
    if(identical(object@dataFull , list())){
        stop("dataFull slot is empty, no data to filter. run getAlterations")
    }
    if(is.null(object@dataFull$mutations$data)){
        stop("No mutation data to filter")
    }
    mode <- mode[1]
    if(mode %notin% c("exclude" , "keep")){
        stop("mode can only be exclude or keep")
    }
    # Take data
    muts <- object@dataFull$mutations$data
    if( !is.null(tumor_type) ){
      mutsToSpare <- muts[ muts$tumor_type %notin% tumor_type , , drop=FALSE]
      muts <- muts[ muts$tumor_type %in% tumor_type , , drop=FALSE]
    } else {
      mutsToSpare <- NULL
    }
    if(nrow(muts)==0){
      object@dataFull$mutations$data <- rbind( muts , mutsToSpare)
      object <- subsetAlterations(object)
      return(object)
    }
    # If filtered is set
    if(!is.null(filtered)){
        if(!is.data.frame(filtered)){
            stop("filtered must be a data.frame")
        }
        filtered <- filtered[ , sort(colnames(filtered))]
        possibleColNames <- list(
            sort(c("gene_symbol" , "amino_acid_change"))
            , sort(c("gene_symbol" , "amino_position"))
            , c("genomic_position")
            )
        #---------------------------------
        # Check DF colnames
        #
        if(any(colnames(filtered) %notin% unlist(possibleColNames))){
            noGood <- setdiff(colnames(filtered) , unlist(possibleColNames))
            stop(paste("The following colnames are not recognized:" , paste(noGood , collapse=", ")))
        }
        rightColNames <- lapply(possibleColNames , function(x) identical(sort(x) , colnames(filtered))) %>% unlist
        if(sum(rightColNames)!=1){
            stop("filtered data.frame do not contain a valid mutation identifier. Check the manual")
        }
        #-----------------------------------
        # NOW LET'S FILTER
        #
        colToBeFiltered <- muts[ , colnames(filtered) , drop=FALSE] %>% apply(. , 1 , paste , collapse="")
        colToFilter <- filtered[ , colnames(filtered) , drop=FALSE] %>% apply(. , 1 , paste , collapse="")
        if(mode == "exclude"){
            myNewRows <- colToBeFiltered %notin% colToFilter
        } else {
            myNewRows <- colToBeFiltered %in% colToFilter
        }
        muts <- muts[ myNewRows , , drop=FALSE]
    }
    # If bed is set
    if(!is.null(bed)){
        # bed <- data.frame(chr = paste0("chr" , c(1:23 , "X" , "M")) , start = sample(10000:100000 , 25) , end = sample(100000:1000000 , 25) , stringsAsFactors=FALSE)
        bed <- .bedCheck(bed)
        bed$start <- bed$start + 1L
        bedGR <- GenomicRanges::makeGRangesFromDataFrame(bed)
        mutsGR <- strsplit(muts$genomic_position , ":") %>% 
                    lapply(. , '[' , c(1,2)) %>% 
                    do.call("rbind" , .)
        mutsGR <- data.frame(chr = mutsGR[ , 1] , start = as.integer(mutsGR[ , 2])-1L , end = as.integer(mutsGR[ , 2]) )
        mutsGR <- GenomicRanges::makeGRangesFromDataFrame(mutsGR)
        IDX <- GenomicRanges::findOverlaps(bedGR , mutsGR) %>% S4Vectors::subjectHits(.)
        IDXtruefalse <- seq_len(nrow(muts)) %in% IDX
        if(mode == "exclude"){
            muts <- muts[ !IDXtruefalse ,  ]
        } else {
            muts <- muts[ IDXtruefalse , ]
        }
    }
    object@dataFull$mutations$data <- rbind( muts , mutsToSpare)
    object <- subsetAlterations(object)
    return(object)
})


