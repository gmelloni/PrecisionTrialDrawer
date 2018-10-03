#-------------------------------------------
# SET OF CLASS AND METHODS ADAPTED FROM LOWMACA FOR panelOptimizer
#-------------------------------------------

.createEmptyLowObj <- function()
    list(
        arguments=list(
            genes=NULL
            , pfam=NULL
            , pfamAllGenes=''
            , input=''
            , mode=''
            , params=list(
                mutation_type=c('missense', 'all', 'truncating', 'silent')[1]
                , tumor_type='all_tumors'
                , min_mutation_number=-1L
                , density_bw=0
                , clustal_cmd='clustalo'
                , use_hmm=FALSE
                , datum=FALSE
                )
            , parallelize=list(
                getMutations=FALSE
                , makeAlignment=TRUE
                )
            )
        )

.newLowMACAsimple <- function(genes=NULL, pfam=NULL)
{
    # check arguments
    if( is.null(genes) & is.null(pfam) )
        stop('Either a gene or a pfam should be specified.')
    if( length(pfam)>1 )
        stop('Only one Pfam ID can be evaluated')
    # initialize the LowMaca object
    # object <- new('LowMACAsimple')
    object <- .createEmptyLowObj()
    # class(object) <- "LowMACAsimple"
    if( !is.null(pfam) ) {
        object$arguments$mode <- 'pfam'
        object$arguments$pfam <- pfam
    } else {
        object$arguments$mode <- 'genes'
        object$arguments$genes <- genes
    }
    mode <- object$arguments$mode
    if( mode == 'genes' )
    {
        # load annotation files
        myAlias <- LowMACAAnnotation::getMyAlias()
        myUni <- LowMACAAnnotation::getMyUni()
        #
        selectedColumns <- c('Gene_Symbol', 'Entrez', 'UNIPROT', 'AMINO_SEQ')
        geneID <- .checkGene_to_geneID(genes, myUni, myAlias)$Gene_Symbol
        selectedRows <- myUni$Gene_Symbol %in% geneID
        genesData <- myUni[selectedRows, selectedColumns]
        genesData$Pfam_ID <- '-'
        genesData$Envelope_Start <- 1
        genesData$Envelope_End <- nchar(genesData$AMINO_SEQ)
        seq_names <- paste(genesData[, 'Gene_Symbol']
            , genesData[, 'Pfam_ID']
            , genesData[, 'Entrez'] 
            , paste(genesData[, 'Envelope_Start'], 
                genesData[, 'Envelope_End'], sep=";") 
            , sep="|")
        rownames(genesData) <- seq_names
    } else {
    ## mode == 'pfam'
        # load annotation files
        myPfam <- LowMACAAnnotation::getMyPfam()
        # select rows and colums to match the input
        selectedColumns <- c('Gene_Symbol', 'Pfam_ID', 'Entrez'
            , 'Envelope_Start', 'Envelope_End', 'UNIPROT', 'Pfam_Fasta')
        selectedRows <- myPfam$Pfam_ID %in% pfam
        if( !any(selectedRows) ) stop('Pfam name is not correct or is not mapped by LowMACAsimple')
        genesData <- myPfam[selectedRows, selectedColumns]
        seq_names <- paste(genesData[, 'Gene_Symbol']
            , genesData[, 'Pfam_ID']
            , genesData[, 'Entrez'] 
            , paste(genesData[, 'Envelope_Start'], 
                genesData[, 'Envelope_End'], sep=";") 
            , sep="|")
        colnames(genesData)[colnames(genesData) == 'Pfam_Fasta'] <- 'AMINO_SEQ'
        rownames(genesData) <- seq_names
        if( !is.null(genes) ) {
            object$arguments$genes <- genes
            object$arguments$pfamAllGenes <- genesData
            # load annotation files
            myAlias <- LowMACAAnnotation::getMyAlias()
            myUni <- LowMACAAnnotation::getMyUni()
            geneID <- .checkGene_to_geneID(genes, myUni, myAlias)$Gene_Symbol
            selectedRows <- genesData$Gene_Symbol %in% geneID
            genesData <- genesData[selectedRows, ]
        }
    }
    genesData[, 'Gene_Symbol'] <- factor(genesData[, 'Gene_Symbol'])
    genesData[, 'Pfam_ID'] <- factor(genesData[, 'Pfam_ID'])
    genesData[, 'Entrez'] <- factor(genesData[, 'Entrez'])
    genesData[, 'UNIPROT'] <- factor(genesData[, 'UNIPROT'])
    ## convert non-canonical amnio acids into their
    ## most similar and canonical one
    ##        U (selenocysteine) -> A (alanine)
    genesData$AMINO_SEQ <- gsub('U','A', genesData$AMINO_SEQ)
    ## return the updated object
    object$arguments$input <- genesData
    return(object)

}

# setGeneric('lmParamsSimple', function(object) standardGeneric('lmParamsSimple'))
.lmParamsSimple <- function(object) {

    return(object$arguments$params)
    }

# setGeneric('lmParamsSimple<-', function(object, value) standardGeneric('lmParamsSimple<-'))
# setReplaceMethod('lmParamsSimple', 'LowMACAsimple', function(object, value) {
#     object$arguments$params <- value
#     # validObject(object)
#     object
# })


# setGeneric('setupSimple', function(object, repos=NULL) standardGeneric('setupSimple'))
.setupSimple <- function(object, repos=NULL) {

    object <- .alignSequencesSimple(object)
    object <- .getMutationsSimple(object, repos=repos)
    object <- .mapMutationsSimple(object)
    return(object)
    }

# setGeneric('alignSequencesSimple', function(object) standardGeneric('alignSequencesSimple'))
.alignSequencesSimple <- function(object) {

    genesData <- object$arguments$input
    object$alignment <- .clustalOAlign(genesData)
    object$alignment$df <- data.frame(
                consensus=strsplit(genesData$AMINO_SEQ,'')[[1]]
                , conservation=rep(1, length(genesData$AMINO_SEQ))
                )
    return(object)
    }

# setGeneric('getMutationsSimple', function(object, repos=NULL) standardGeneric('getMutationsSimple'))
.getMutationsSimple <- function(object, repos=NULL) {

    genes <- object$arguments$input
    mutation_type <- object$arguments$params$mutation_type
    tumor_type <- object$arguments$params$tumor_type
    parallelGetMut <- object$arguments$parallelize$getMutationsSimple
    # outputFolder <- object$arguments$paths$output_folder
    gmOut <- .getLocalGeneMutations(genes, 
            mutation_type=mutation_type, tumor_type=tumor_type, 
            localData=repos)
    object$mutations$data <- gmOut$Mutations
    object$mutations$freq <- gmOut$AbsFreq
    return(object)
    }

# setGeneric('mapMutationsSimple', function(object) standardGeneric('mapMutationsSimple'))
.mapMutationsSimple <- function(object) {

    mut <- object$mutations$data
    if( nrow(mut)==0 )
        return(object)
    alignment <- object$alignment$ALIGNMENT
    alignmentLength <- max(alignment$Align)
    # # elaborate the mutations and map them with the alignment
    mut_aligned <- merge( mut , alignment[!is.na(alignment$Ref) , ]
                    , by.x=c("Entrez" , "Gene_Symbol" , "Amino_Acid_Position") 
                    , by.y=c("Entrez" , "Gene_Symbol" , "Ref") 
                    , all.x=TRUE)
    mut_aligned <- mut_aligned[!is.na(mut_aligned$Align), ]
    # parameters for the analysis [can't be set by the user, so far]
    splitCriterion <- c('Gene_Symbol', 'Tumor_Type', 'domainID')[1]
    # mode <- object$arguments$mode
    # if( mode == 'genes' ) {
    #     splitCriterion <- splitCriterion[1]
    # } else splitCriterion <- splitCriterion[3]
    
    # split data
    mut_aligned.split <- split(mut_aligned, mut_aligned[[splitCriterion]])
    mut_aligned.extended <- matrix(0, nrow=length(mut_aligned.split)
        , ncol=alignmentLength)
    for(i in seq_len(length(mut_aligned.split))) {
        tmp <- lengths(
            split(mut_aligned.split[[i]]$Align, mut_aligned.split[[i]]$Align)
            )
        if( length(tmp)>0 ) 
            mut_aligned.extended[i,as.numeric(names(tmp))] <- tmp
    }
    rownames(mut_aligned.extended) <- names(mut_aligned.split)
    # filter data based on number of mutations
    minNMut <- object$arguments$params$min_mutation_number
    tokeep <- rowSums(mut_aligned.extended) >= minNMut
    if(any(!tokeep))
    {
        warning(
            paste("We excluded these genes (or domains) because they have less than"
                , minNMut, "mutations")
            , immediate.=TRUE)
        excludedGenes <- rownames(mut_aligned.extended[!tokeep, ])
        print(excludedGenes)
        object$mutations$excluded <- excludedGenes
    }
    mut_aligned.extended <- mut_aligned.extended[tokeep, , drop=FALSE]
    object$mutations$aligned <- mut_aligned.extended
    return(object)

    }


# setGeneric('entropySimple', function(object, bw=NULL ) standardGeneric('entropySimple'))
.entropySimple <- function(object, bw=NULL ) {

    # if(conservation > 1 || conservation < 0)
        # stop("Conservation threshold must be a number between 0 and 1")
    conservation <- 0
    myMut <- object$mutations
    if(identical(myMut , list()))
        stop("mutations slot is empty. Launch the method getMutationsSimple first!")
    myAln <- object$alignment
    if(identical(myAln , list()))
        stop("alignment slot is empty. Launch the method alignSequencesSimple first!")
    if( nrow(object$mutations$data)==0 ) {
        object$entropy$absval <- NA
        object$entropy$log10pval <- NA
        object$entropy$pvalue <- NA
        return(object)
    }
    mut_extended <- object$mutations$aligned
    alignment    <- object$alignment
    # Uniform variable
    message("Making uniform model...")
    if( is.null(bw) ) bw <- object$arguments$params$density_bw else 
        object$arguments$params$density_bw <- bw
    if( bw=='auto' ) 
    {
        # calculate the bw from the global profile
        bw <- .profileDensity(colSums(mut_extended))$bw
    }
    message(paste('Assigned bandwidth:', round(bw,2)))
    object$entropy$bw <- bw
    weights <- .alnWeights(alignment)
    # if( plotUniform ) pdf(file.path(outputFolder , "Uniform_Model.pdf"))
    object$entropy$uniform <- .makeUniformModel(mut_extended, bw=bw, nboot=1000, 
        weights=weights, plotOUT=FALSE)
    # if( plotUniform ) dev.off()

    # Calculate the entropy values
    globalProfile <- colSums(mut_extended)
    uniform <- object$entropy$uniform
    absval <- .profileEntropy(globalProfile, norm=FALSE, bw=bw)
    log10pval <- .profileEntropy(globalProfile, norm=TRUE, bw=bw, model=uniform)
    pvalue <- 10^log10pval

    object$entropy$absval <- absval
    object$entropy$log10pval <- log10pval
    object$entropy$pvalue <- pvalue
    object$entropy$conservation_thr <- conservation

    # null profile
    # Null profile calculated on the global profile
    nullOut <- .makeNullProfile(mut_extended, bw=bw, 
        nboot=1000, weights=.alnWeights(alignment))

    pvals <- nullOut$pvalue
    filteredPvals <- pvals
    filteredPvals[object$alignment$df$conservation < conservation] <- NA
    qvals <- p.adjust(filteredPvals, method='BH')

    object$alignment$df <- cbind(
        object$alignment$df[, c('consensus', 'conservation')]
        , nullOut, qvalue=qvals
        #, misalnFreq=object$alignment$df[, 'misalnFreq']
        )
    return(object)
    }


# setGeneric('lfmSimple', function(object, metric='qvalue', threshold=.05)
    # standardGeneric('lfmSimple'))
.lfmSimple <- function(object, metric='qvalue', threshold=.05) {

    # if(is.null(conservation))
    #     conservation <- object$entropy$conservation_thr
    #Checks on parameters
    if(! (metric %in% c("qvalue" , "pvalue") ))
        stop("metric parameter can be only 'pvalue' or 'qvalue'")
    if(!(threshold <= 1 && threshold>=0))
        stop("pvalue/qvalue threshold must be a number between 0 and 1")
    # if(!(conservation <= 1 && conservation>=0))
    #     stop("conservation score threshold must be a number between 0 and 1")
    entrop <- object$entropy
    if(identical(entrop , list()))
        stop("Entropy slot is empty. Launch the method entropy first!")
    if( nrow(object$mutations$data)==0 ) {
        message('No mutations available for this object.')
    } else {
        if( metric == 'qvalue' ) {
            # if( conservation != object$entropy$conservation_thr ) {
            #     ## in case a new conservation threshold is given,
            #     ## recalculate the qvalue
            #     filteredPvals <- object$alignment$df$pvalue
            #     filteredPvals[object$alignment$df$conservation < conservation] <- NA
            #     qvalue <- p.adjust(filteredPvals, method='BH')
            # } else {
                ## use qvalue already stored
                qvalue <- object$alignment$df$qvalue
            # }
            signifPos <- which(qvalue < threshold)
        } else {
            ## use pvalue
            pvalue <- object$alignment$df$pvalue
            signifPos <- which(pvalue < threshold)
        }
        out <- lapply(signifPos, function(pos) {
            selectedRows <- object$alignment$ALIGNMENT$Align == pos
            alnData <- object$alignment$ALIGNMENT[selectedRows, ]
            selectedColumns <- c('Gene_Symbol', 'Ref', 'Envelope_Start'
                , 'Envelope_End')
            alnData <- alnData[!is.na(alnData$Ref), selectedColumns]
            colnames(alnData) <- c('Gene_Symbol', 'Amino_Acid_Position'
                , 'Envelope_Start', 'Envelope_End')
            selectedColumns <- c('Gene_Symbol', 'Amino_Acid_Position'
                , 'Amino_Acid_Change', 'Sample', 'Tumor_Type')
            mutData <- object$mutations$data[, selectedColumns]
            lfmSimple <- merge(mutData, alnData)
            if( nrow(lfmSimple)>0 ){
                lfmSimple$Multiple_Aln_pos <- pos
                lfmSimple$metric <- switch(metric, 'pvalue'=pvalue[pos], 'qvalue'=qvalue[pos])
            } else {
                lfmSimple$Multiple_Aln_pos <- character(0)
                lfmSimple$metric <- numeric(0)
            }
            return(lfmSimple)
            })
        out <- do.call('rbind', out)
        # load data
        myUni <- LowMACAAnnotation::getMyUni()
        selectedColumns <- c('Gene_Symbol', 'Entrez', 'Entry', 'UNIPROT'
            , 'Chromosome', 'Protein.name')
        myUniSmall <- myUni[, selectedColumns]
        out <- merge(out, myUniSmall)
        return(out)
    }
    }




# setGeneric('nullProfileSimple', function(object) standardGeneric('nullProfileSimple'))
.nullProfileSimple <- function(object) {


    if( length(object$alignment)==0 )
        stop('Perform alignment on the object before plotting.')

    if( length(object$mutations)==0 )
        stop('Retrieve mutations for the object before plotting.')

    if( length(object$entropy)==0 )
        stop('Perform entropy calculation on the object before plotting.')

    if( nrow(object$mutations$data)==0 ) {
        message('No mutations available for this object.')
    } else {

        # if( is.null(windowlimits) )
            windowlimits <- seq_len(nrow(object$alignment$df))
        # if(is.null(conservation))
            conservation <- object$entropy$conservation_thr
        mean <- object$alignment$df$mean#[windowlimits]
        lowerThreshold <- object$alignment$df$lTsh#[windowlimits]
        upperThreshold <- object$alignment$df$uTsh#[windowlimits]
        profile <- object$alignment$df$profile#[windowlimits]
        pvalue <- object$alignment$df$pvalue#[windowlimits]
        # if( conservation != object$entropy$conservation_thr ) {
        #     ## in case a new conservation threshold is given,
        #     ## recalculate the qvalue
        #     filteredPvals <- object$alignment$df$pvalue
        #     filteredPvals[object$alignment$df$conservation < conservation] <- NA
        #     qvalue <- p.adjust(filteredPvals, method='BH')
        # } else {
            ## use qvalue already stored
            qvalue <- object$alignment$df$qvalue
        # }
        qvalue <- qvalue#[windowlimits]
        qvalSignif_x <- which(qvalue < 5e-2)
        qvalSignif_y <- profile[qvalSignif_x] + max(profile)/20
        max_y <- max(c(
            object$alignment$df$lTsh,
            object$alignment$df$uTsh,
            object$alignment$df$profile*1.1
            ), na.rm=TRUE)
        blackProfile <- vapply(seq_len(length(profile)),function(i) {
          min(profile[i], upperThreshold[i])
          } , numeric(1))
        orangeProfile <- profile - blackProfile
        mp <- barplot(rbind(blackProfile, orangeProfile), 
            col=c('black', 'orange'), ylim=c(0, max_y), ylab='Mutation density'
            )
        axis(1, at=mp, labels=windowlimits, las=2)
        lines(mp, upperThreshold, lty=2, lwd=2, col='blue')
        qvalSignif_x <- mp[qvalue < 5e-2]

        ### plot an asterisk above the residues which are significant in 
        # terms of the pvalue
        if( length(qvalSignif_x) > 0 ) {
            text(qvalSignif_x, qvalSignif_y, '*', cex=2, col='red')
        }
    }
    }

# setGeneric('lmPlotSimple', function(object, conservation=NULL , splitLen=NULL) standardGeneric('lmPlotSimple'))
.lmPlotSimple <- function(object, conservation=NULL , splitLen=NULL) {

    opar <- par("mfrow" , "mar")
    on.exit(par(opar))
    if( length(object$alignment)==0 )
        stop('Perform alignment on the object before plotting.')

    if( length(object$mutations)==0 )
        stop('Retrieve mutations for the object before plotting.')

    if( length(object$entropy)==0 )
        stop('Perform entropy calculation on the object before plotting.')

    if( nrow(object$mutations$data)==0 ) {
        message('No mutations available for this object.')
    } else {
        if(is.null(conservation))
            conservation <- object$entropy$conservation_thr
        if( is.null(splitLen) )
            splitLen <- ncol(object$mutations$aligned)

        windowlimits <- seq_len(ncol(object$mutations$aligned))
        windowlimitsSplit <- split(windowlimits, ceiling(windowlimits/splitLen))

        mode <- object$arguments$mode
        nObj <- nrow(object$arguments$input)

        par(mar=c(2,4,2,2))
        # if( nObj > 1 ) {
        #     layoutMatrix <- as.matrix(c(1,1,2,2,3,4,4))
        #     multipleMat <- as.list(rep(NA), length(windowlimitsSplit))
        #     for( i in 1:length(windowlimitsSplit) ) 
        #         multipleMat[[i]] <- layoutMatrix+(max(layoutMatrix)*(i-1))
        #     layoutMatrix <- do.call('rbind', multipleMat)
        #     plotType <- 4
        # } else if( mode == 'genes' ) {
            layoutMatrix <- as.matrix(c(1,1,1,1,3,2,2,2,2))
            multipleMat <- as.list(rep(NA), length(windowlimitsSplit))
            for( i in seq_len(length(windowlimitsSplit) )){ 
                multipleMat[[i]] <- layoutMatrix+(max(layoutMatrix)*(i-1))
            }
            layoutMatrix <- do.call('rbind', multipleMat)
            plotType <- 3
        # } else {
        #     layoutMatrix <- as.matrix(c(1,1,2,2))
        #     multipleMat <- as.list(rep(NA), length(windowlimitsSplit))
        #     for( i in 1:length(windowlimitsSplit) ) 
        #         multipleMat[[i]] <- layoutMatrix+(max(layoutMatrix)*(i-1))
        #     layoutMatrix <- do.call('rbind', multipleMat)
        #     plotType <- 2
        # }
        layout(layoutMatrix)

        # for( i in 1:length(windowlimitsSplit) ) {

        #     windowlimits <- windowlimitsSplit[[i]]

            # origMAlign <- object$alignment$CLUSTAL
            # m <- consensusMatrix(origMAlign)
            # motif <- pcm2pfm(m[, windowlimits])
            # motif <- new('pfm', mat=motif, name='', color=colorset(alphabet='AA'))
            # motif$color <- c("-"="#FFFFFF" , motif$color)
            #m <- consensusMatrix(origMAlign)
            #motif <- pcm2pfm(m[, windowlimits])
            #motif <- new('pfm', mat=motif, name='', color=colorset(alphabet='AA'))
            log10pval <- object$entropy$log10pval
            pvalue <- object$entropy$pvalue

            ## plot 1
            if( plotType == 3) par(mar=c(0,4,4,2))
            par(mar=c(2,4,4,2))
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
            lenAln <- ncol(object$mutations$aligned)
            mut_aligned <- object$mutations$aligned
            # if( plotType == 3 ) colnames(mut_aligned) <- 1:length(mut_aligned)
            #print(nObj)
            mp <- barplot(
                mut_aligned
                , col=myPalette(nrow(object$mutations$aligned))
                , border=if(lenAln<50) 'black' else NA
                , main=paste("Log10 P-Value:" , round(log10pval , 2) , 
                    "P-Value:" , signif( pvalue , 2 ) , "Bw:" , signif(object$entropy$bw,3) 
                    )
                , ylab='Mutation #'
                , ylim=c(0, max(colSums(object$mutations$aligned)))
                )
            axis(1, at=mp, labels=windowlimits, las=2)
            par(mar=c(2,4,2,2))
            ## plot 2
            
            # if( plotType == 3) par(mar=c(2,4,0,2))
            # nullProfileSimple(object, conservation=conservation , windowlimits)
            .nullProfileSimple(object)
            # if( plotType == 3) {
            #     pSeq <- round(seq(1, nchar(object$arguments$input$AMINO_SEQ), length.out=20))
            #     axis(1, at=pSeq, labels=pSeq)
            # }
            
            ## if there is more than one object also plot 
            ## conservation plots    
            # if( nObj > 1 ) {
            #     ## plot 3
            #     mp <- barplot(object$alignment$df$conservation#[windowlimits]
            #         , col='darkgoldenrod1', ylab='Conservation', ylim=c(0,1))
            #     axis(1, at=mp, labels=windowlimits, las=2)
            #     ## plot 4
            #     plotMotifLogo(motif, p=motif$background
            #         , colset=motif$color[rownames(motif$mat)]
            #         , ylab='Logo', xaxis=FALSE
            #         )
            #     axis(1,at=1/length(windowlimits)*(seq_along(windowlimits)-.5)
            #         , labels=windowlimits,las=2)
            #     #plot(motif, ylab='Logo')
            #     # plot.new()
            # }

            #############
            ## in case the analysis is on a single gene 
            ## plot its domains
            ############### 
            if( plotType == 3 ) {

                myPfam <- LowMACAAnnotation::getMyPfam()
                gene <- as.character(object$arguments$input$Gene_Symbol)
                domains <- myPfam[myPfam$Gene_Symbol==gene
                    , c("Envelope_Start" , "Envelope_End" , "Pfam_Name")]
                ## create empty plot
                plot.new()
                plot.window(
                    xlim=range(windowlimits)
                    , ylim=c(0,0.05)
                    )
                par(mar=c(0,0,0,0))
                ## plot domains
                if(nrow(domains)>0) {
                    for( i in seq_len(nrow(domains)) ) {
                        int <- intersect(domains$Envelope_Start[i]:domains$Envelope_End[i], windowlimits)
                        if( length(int)>0 ) {
                            domains$Envelope_Start[i] <- min(int)
                            domains$Envelope_End[i] <- max(int)
                        } else {
                            domains <- domains[-i,]
                        }
                    }
                    for(i in seq_len(nrow(domains))) {
                        xleft=as.numeric(domains[i , "Envelope_Start"])
                        xright=as.numeric(domains[i , "Envelope_End"])
                        ytop=0.05
                        ybottom=0
                        col=topo.colors(nrow(domains) , alpha=0.5)[i]
                        characters <- nchar(domains[i , "Pfam_Name"])
                        rect(xleft=xleft , xright=xright 
                            , ytop=ytop , ybottom=ybottom 
                            , col=col )
                        if(characters<=3){
                            text(x=(xright+xleft)/2 , y=0.025 
                                , domains[i , "Pfam_Name"] 
                                , font=2)
                        } else {
                            if((xright-xleft)<=100){
                                text(x=(xright+xleft)/2 , y=0.025 
                                    , domains[i , "Pfam_Name"] 
                                    , font=2 , cex=0.8)
                            } else {
                                text(x=(xright+xleft)/2 , y=0.025 
                                    , domains[i , "Pfam_Name"] 
                                    , font=2)
                            }
                        }
                    }
                } else {
                    text(x=length(windowlimits)/2, y=0.025
                        , 'no pfam domains within gene sequence', cex=1.5)
                }
                par(mar=c(2,4,4,2))
            }
        # }

    }

    }