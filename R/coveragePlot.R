.tumorFreqPlotter <- function(matToPlot , colNum , cex.main , tumor.freqs){
    opar <- par("mfrow" , "mar")
    on.exit(par(opar))
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    n_plots <- nrow(matToPlot)
    used_colors <- c()
    if(is.null(colNum)){
        # Automatic detection of plot best layout
        par(mfrow=.mfrow(n_plots))
    } else {
        # Check if colNum is valid (less or equal to plots to display)
        if(colNum <= n_plots)
            par(mfrow=c( ceiling(n_plots/colNum) , colNum))
        else
            stop("Number of columns asked are more than plots to display")
    }
    sampsTitle <- paste("Freqs:" 
                , paste( paste(names(tumor.freqs) , round(tumor.freqs,3) 
                    , sep=" - ") , collapse=" , ")
                )
    for(i in seq_len(nrow(matToPlot))){
        tabToPlot <- matToPlot[i , , drop=FALSE]
        cols <- if(i<=2){
                    myPalette(i)[i]
                } else {
                    setdiff( myPalette(i) , used_colors )[1]
                }
        used_colors <- c(used_colors , cols)
        if(cex.main=="auto"){
            mycexmain <- if(n_plots>10){
                        .5
                    }else if(n_plots>5){
                        .8
                    } else {
                        1
                    }
        } else {
            mycexmain <- cex.main
        }
        tryCatch({
            xx <- barplot(tabToPlot
                , ylim=c(0,1) 
                , col=cols
                , main=paste( rownames(matToPlot)[i] , sampsTitle , sep="\n" )
                , cex.main=mycexmain
                , xlab="Incremental number of alteration per sample"
                # , ylab="Number of samples with at least # alterations"
                # ,beside=TRUE
                )}
            , error=function(e) {
                if(e$message=="figure margins too large")
                    stop(paste("Too many plots to display. ",
                        "Reduce the margins or better,",
                        " redirect the output to a pdf/png file"))
                else
                    stop(e)
            }
            )
        mycextext <- if(n_plots>20){
                        .5
                    }else if(n_plots>10){
                        .7
                    } else {
                        1
                    }
        # Add text at top of bars
        text(x = xx+0.4 , y = tabToPlot+(1/10)
            , label = paste0(100*(round(tabToPlot , 3)) , "%")
            , pos = 3, cex = mycextext
            , col = "red" , srt=90)
    }
}

# Barplot of number of samples covered by at least 1, at least 2 etc.
# This barplot can be customized including:
    # custom data type (muts, expr, cna)
    # custom subset of the data 
        # (by drug, by group, by gene, by data type , by tumor and by nothing)
    # set the maximum number of alteration to display
    # if you don't want the plot, you get the data
setGeneric('coveragePlot', function(object
    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
    , grouping=c(NA , "drug" , "group" , "gene_symbol" 
        , "alteration_id" , "tumor_type")
    , tumor_type=NULL
    , alterationType.agg=TRUE
    , collapseMutationByGene=TRUE
    , collapseByGene=FALSE
    , tumor.weights=NULL
    , tumor.freqs=NULL
    , maxNumAlt=10
    , colNum=NULL
    , cex.main="auto"
    , noPlot=FALSE) {
    standardGeneric('coveragePlot')
    })
setMethod('coveragePlot', 'CancerPanel', function(object
    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
    , grouping=c(NA , "drug" , "group" , "gene_symbol" 
        , "alteration_id" , "tumor_type")
    , tumor_type=NULL
    , alterationType.agg=TRUE
    , collapseMutationByGene=TRUE
    , collapseByGene=FALSE
    , tumor.weights=NULL
    , tumor.freqs=NULL
    , maxNumAlt=10
    , colNum=NULL
    , cex.main="auto"
    , noPlot=FALSE)
{
    # Checks
    possibleAlterations <- c("copynumber" , "expression" 
        , "mutations" , "fusions")
    possibleGrouping <- c(NA , "drug" , "group" 
        , "gene_symbol" , "alteration_id" , "tumor_type")
    # collapseByGene is more stringent than collapseMutationByGene
    # If the first is set, the other one should not be evaluated
    if(collapseByGene & collapseMutationByGene){
        collapseMutationByGene <- FALSE
    }
    if(any(alterationType %notin% possibleAlterations)){
        stop(paste("alterationType can only be one or more of the following" ,
            paste(possibleAlterations , collapse=", ")))
    }
    if( length(grouping)>1 & any(is.na(grouping)) ){
        # stop("You cannot choose grouping NA and another grouping factor")
        grouping <- NA
    }
    if(!any(is.na(grouping))){
        if(any(grouping %notin% possibleGrouping))
            stop(paste("grouping can only be one of the following:" ,
                paste(possibleGrouping , collapse=", ")))
    }
    if(("alteration_id" %in% grouping) & collapseByGene){
        warning(paste("If 'alteration_id' is in grouping variables,",
            " you cannot collapse by gene. The option was set to FALSE"))
        collapseByGene <- FALSE
    }
    if(is.null(alterationType.agg)){
        alterationType.agg <- TRUE
    } else {
        if( alterationType.agg %notin% c(TRUE,FALSE)){
            stop("alterationType.agg must be TRUE or FALSE")
        }
    }
    # Check plotting parameters
    if(!is.null(maxNumAlt)){
        if(!is.numeric(maxNumAlt)){
            stop("maxNumAlt must be a positive integer")
        }
        if(maxNumAlt<=0){
            stop("maxNumAlt must be a positive integer")
        }
        maxNumAlt <- ceiling(maxNumAlt)
    }
    if(!is.null(colNum)){
        if(!is.numeric(colNum)){
            stop("colNum must be a positive integer")
        }
        if(colNum<=0){
            stop("colNum must be a positive integer")
        }
        colNum <- ceiling(colNum)
    }
    if(is.null(cex.main)){
        cex.main <- "auto"
    } else {
        if(length(cex.main)>1 | (!is.numeric(cex.main) & cex.main!="auto")){
            stop("cex.main must be a numeric value or 'auto'")
        }
    }
    # Check tumor.weights consistency
    # tumor.freqs is a named vector of integers: e.g. c(brca=100 , luad=1000)
    if(!is.null(tumor.weights)){
        .tumor.weights.standardCheck(tumor.weights 
            , tumor.freqs , object , tumor_type)
        if(!alterationType.agg){
            warning(paste("if tumor.weights are in use,",
                " alterationType.agg is set to the default TRUE"))
            alterationType.agg <- TRUE
        }
    }
    # Check tumor.freqs consistency
    # tumor.freqs is a named vector of 0-1 coefficient that sum to 1
    # e.g. c(brca=0.1 , luad=0.9)
    if(!is.null(tumor.freqs)){
        .tumor.freqs.standardCheck(tumor.weights 
            , tumor.freqs , object , tumor_type)
        if("tumor_type" %in% grouping){
            stop(paste("If you use tumor.freqs,",
                " tumor_type cannot be in grouping parameter"))
        }
    }

    #----------------------------
    # GRAB DATA AND SAMPLES
    #----------------------------

    de <- dataExtractor(object=object , alterationType=alterationType 
                        , tumor_type=tumor_type 
                        , collapseMutationByGene=collapseMutationByGene 
                        , collapseByGene=collapseByGene 
                        , tumor.weights=tumor.weights)
    mydata <- de$data
    mysamples <- de$Samples
    tum_type_diff <- de$tumor_not_present
    rm(de)

    # If tumor.freqs is set, 
    # we basically run this method in recursive mode for each tumor type
    # and then aggregate everything
    if(!is.null(tumor.freqs)){
        covRecurse <- lapply(names(tumor.freqs) , function(tum){
            # browser()
                out <- tryCatch( suppressWarnings(
                    coveragePlot(object
                    , alterationType=alterationType
                    , grouping=grouping
                    , tumor_type=tum
                    , alterationType.agg=alterationType.agg
                    , collapseMutationByGene=collapseMutationByGene
                    , collapseByGene=collapseByGene
                    , tumor.weights=NULL
                    , tumor.freqs=NULL
                    , maxNumAlt=maxNumAlt
                    , colNum=NULL
                    , cex.main="auto"
                    , noPlot=TRUE
                        )) , error=function(e) return(NULL))
                if(is.null(out)){
                    return(NULL)
                }
                if(nrow(out$plottedTable)==1){
                    out2 <- (out$plottedTable/out$Samples)*tumor.freqs[tum]
                } else {
                    out2 <- apply(out$plottedTable , 2 , function(x) {
                        as.matrix(x)/out$Samples} ) * tumor.freqs[tum]
                }
                rownames(out2) <- rownames(out$plotted)
                return(out2)
            })
        rows <- unique( lapply(covRecurse , rownames) %>% unlist)
        cols <- unique( lapply(covRecurse , colnames) %>% unlist)
        matToPlot <- matrix(NA , nrow=length(rows) , ncol=length(cols))
        rownames(matToPlot) <- rows
        colnames(matToPlot) <- cols
        for(i in rows){
            for(j in cols){
                vals <- lapply(covRecurse , function(x) {
                    if(i %in% rownames(x) && j %in% colnames(x)){
                        out <- x[i , j]
                    } else {
                        out <- 0
                    }
                    return(out)
                }) %>% unlist
                matToPlot[i , j] <- sum(vals)
            }
        }
        if(noPlot){
            return( list(plottedTable=matToPlot , Samples=NULL) )
        } else {
            return(.tumorFreqPlotter(matToPlot 
                , colNum , cex.main , tumor.freqs))
        }

    }
    # Save current global e graphic options before 
    # modifying them and restore them upon exiting
    opar <- par("mfrow" , "mar")
    on.exit(par(opar))
    if(!any(is.na(grouping))) {
        splitterVar <- apply(mydata[ , grouping , drop=FALSE] 
            , 1 , function(x) paste(x , collapse=":"))
        mydata_split <- split(mydata , splitterVar)
        plottedDataNames <- c()
        n_pats_vec <- c()        
        if(!noPlot){
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral"))
                , space="Lab")
            n_plots <- length(mydata_split)
            used_colors <- c()
            if(is.null(colNum)){
                # Automatic detection of plot best layout
                par(mfrow=.mfrow(n_plots))
            } else {
                # Check if colNum is valid (less or equal to plots to display)
                if(colNum <= n_plots)
                    par(mfrow=c( ceiling(n_plots/colNum) , colNum))
                else
                    stop("Number of columns are more than plots to display")
            }
        }
        for(i in seq_len(length(mydata_split))){
            df <- mydata_split[[i]]
            df_name_split <- strsplit(names(mydata_split)[i] , ":") %>% unlist
            if(length(df_name_split)!=length(grouping)){
                df_name_split <- c(df_name_split , "NA")
            }
            df_name_split[df_name_split==""] <- "NA"
            myTitle <- paste(paste(grouping , df_name_split , sep=": ") 
                , collapse="  ") #%>% toupper
            plottedDataNames <- c(plottedDataNames , myTitle)
            if(nrow(df)==0) {
                tabToPlot <- rep(0 , maxNumAlt)
                names(tabToPlot) <- seq_len(maxNumAlt)
                if(i==1)
                    plottedData <- list(tabToPlot)
                else
                    plottedData <- c(plottedData , list(tabToPlot))
                if(!noPlot){
                    par(mar = c(0,0,0,0))
                    plot(c(0, 1), c(0, 1), ann = FALSE
                        , bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                    text(x = 0.5, y = 0.5
                        , paste(myTitle , "\n", "No alterations to show"), 
                                cex = 1.6, col = "black")
                }
                next
            }
            ########### PAY ATTENTION:
            # If "tumor_type" is in grouping, 
            # the reference set of samples must be the one of the tumor type
            # If "alteration_id" is in grouping, we have to choose:
                # only the samples of the specific alterationType/tumor_type
                # the samples tested for ALL alterationTypes 
                    # with data on the specific alteration_id
                # if brca was tested on 3000 pats for mutations and 1000 for cna
                    # alterationType.agg=TRUE -> common samples btw muts and cna
                    # alterationType.agg=FALSE -> 3000 for mut, 1000 for cna
            if("tumor_type" %in% grouping){
                tumTypes <- df_name_split[grouping=="tumor_type"]
            } else {
                tumTypes <- names(mysamples)[names(mysamples)!="all_tumors"]
            }
            if("tumor_type" %in% grouping & !"alteration_id" %in% grouping){
                pats <- mysamples[tumTypes] %>% unlist %>% unique
            } else if(all(c("tumor_type","alteration_id") %in% grouping)){
                if(alterationType.agg){
                    pats <- mysamples[tumTypes] %>% unlist %>% unique
                } else {
                    altID <- df_name_split[grouping=="alteration_id"]
                    altID <- .mapvalues(from=c("mut" , "cna" , "fus" , "expr") 
                        , to=c("mutations" , "copynumber" 
                            , "fusions" , "expression")
                        , altID
                        , warn_missing=FALSE)
                    pats <- object@dataFull[[altID]]$Samples[[tumTypes]]
                }
            } else if(!"tumor_type" %in% grouping & 
                "alteration_id" %in% grouping){
                if(alterationType.agg){
                    pats <- mysamples$all_tumors
                } else {
                    altID <- df_name_split[grouping=="alteration_id"]
                    altID <- .mapvalues(from=c("mut" , "cna" , "fus" , "expr") 
                        , to=c("mutations" , "copynumber" 
                            , "fusions" , "expression")
                        , altID
                        , warn_missing=FALSE)
                    pats <- unlist(object@dataFull[[altID]]$Samples) %>% unique
                }
            } else {
                pats <- mysamples$all_tumors
            }
            n_pats <- length(pats)
            n_pats_vec <- c(n_pats_vec , n_pats)
            myTitle <- c(myTitle 
                        , paste0("Number of Samples: " 
                            , n_pats
                            , " \n(" 
                            , paste(tumTypes , collapse=", ") 
                            , ")"))
            n_alts_bypat <- table(factor(df$case_id , levels=pats))
            tabToPlot <- vapply(seq_len(maxNumAlt) , function(x) {
                                length(n_alts_bypat[n_alts_bypat>=x])
                            } , numeric(1))
            names(tabToPlot) <- seq_len(maxNumAlt)
            # Collect all the data to be plotted
            if(i==1)
                plottedData <- list(tabToPlot)
            else
                plottedData <- c(plottedData , list(tabToPlot))
            if(!noPlot){
                stat1 <- paste0(round(mean(n_alts_bypat),3) 
                    , "\u00b1 " 
                    , round(sd(n_alts_bypat),3)
                    , " alterations per sample")
                myTitle <- c(myTitle 
                            , stat1
                            )
                maxylim=max(tabToPlot) + round(max(tabToPlot)/2.5)
                if(maxylim==1)
                    maxylim <- 2
                cols <- if(i<=2){
                            myPalette(i)[i]
                        } else {
                            setdiff( myPalette(i) , used_colors )[1]
                        }
                used_colors <- c(used_colors , cols)
                if(cex.main=="auto"){
                    mycexmain <- if(n_plots>10){
                                .5
                            }else if(n_plots>5){
                                .8
                            } else {
                                1
                            }
                } else {
                    mycexmain <- cex.main
                }
                tryCatch({
                    xx <- barplot(tabToPlot
                        , ylim=c(0,maxylim) 
                        , col=cols
                        , main=paste(myTitle , collapse="\n")
                        , cex.main=mycexmain
                        , xlab="Incremental number of alteration per sample"
                        )}
                    , error=function(e) {
                        if(e$message=="figure margins too large")
                            stop(paste0("Too many plots to display. ",
                                    "Reduce the margins or better, ",
                                    "redirect the output to a pdf/png file"))
                        else
                            stop(e)
                    }
                    )
                mycextext <- if(n_plots>20){
                                .5
                            }else if(n_plots>10){
                                .7
                            } else {
                                1
                            }
                # Add text at top of bars
                text(x = xx+0.4 , y = tabToPlot+(maxylim/10)
                    , label = paste0(100*(round(tabToPlot/n_pats , 3)) , "%")
                    , pos = 3, cex = mycextext
                    , col = "red" , srt=90)
            }
        }
        if(noPlot){
            out <- do.call("rbind" , plottedData)
            rownames(out) <- plottedDataNames
            return(list(plottedTable=out , Samples=n_pats_vec))
        }
    } else {
        if(nrow(mydata)==0) {
            if(!noPlot){
                par(mar = c(0,0,0,0) , mfrow=c(1,1))
                plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n'
                    , type = 'n', xaxt = 'n', yaxt = 'n')
                text(x = 0.5, y = 0.5
                    , paste(myTitle , "\n", "No alterations to show"), 
                            cex = 1.6, col = "black")
            } else {
                # A lot of useless code but it is just for coherence
                tabToPlot <- rep(0 , maxNumAlt)
                names(tabToPlot) <- seq_len(maxNumAlt)
                plottedData <- list(tabToPlot)
                out <- do.call("rbind" , plottedData)
                rownames(out) <- "noGrouping"
                #suppressWarnings(par(originalPar))
                return(out)
            }
        }
        # tumTypes <- unique(mydata$tumor_type)
        tumTypes <- names(mysamples)[names(mysamples)!="all_tumors"]
        pats <- mysamples[tumTypes] %>% unlist %>% unique
        n_pats <- length(pats)
        myTitle <- paste0("Number of Samples: " 
                          , n_pats 
                          ," \n(" 
                          , paste(tumTypes , collapse=", ") 
                          , ")")
        n_alts_bypat <- table(factor(mydata$case_id , levels=pats))
        tabToPlot <- vapply(seq_len(maxNumAlt) , function(x) {
                            length(n_alts_bypat[n_alts_bypat>=x])
                        } , numeric(1))
        names(tabToPlot) <- seq_len(maxNumAlt)
        plottedData <- list(tabToPlot)
        if(!noPlot) {
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral"))
                , space="Lab")
            stat1 <- paste0(round(mean(n_alts_bypat),3) 
                , " \u00b1 "
                , round(sd(n_alts_bypat),3)
                , "alterations per sample")
            myTitle <- c(myTitle 
                        , stat1 )
            maxylim=max(tabToPlot) + round(max(tabToPlot)/2.5)
            if(maxylim==1)
                maxylim <- 2
            cols <- myPalette(1)            
            par(mfrow=c(1,1))
            xx <- barplot(tabToPlot
                , ylim=c(0,maxylim) 
                , col=cols
                , main=paste(myTitle , collapse="\n")
                , cex.main=1
                , xlab="Incremental number of alteration per sample"
                )
            # Add text at top of bars
            text(x = xx+0.4 , y = tabToPlot+(maxylim/10)
                , label = paste0(100*(round(tabToPlot/n_pats , 3)) , "%")
                , pos = 3, cex = 1
                , col = "red" , srt=90)
        } else {
            out <- do.call("rbind" , plottedData)
            rownames(out) <- "noGrouping"
            return(list(plottedTable=out , Samples=n_pats))
        }
    }
})
