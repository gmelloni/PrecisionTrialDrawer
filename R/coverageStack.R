.tumorFreqPlotterStack <- function(finalTab , var , grouping , alterationType , tumor.freqs){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    myTitle <- paste("Number of covered samples:" , var , "BY" , grouping)
    alts <- paste("Data Type:" , paste(alterationType , collapse=", "))
    samps <- paste("Freqs:" 
                , paste( paste(names(tumor.freqs) , round(tumor.freqs , 3) , sep=" - ") , collapse=" , ")
                )
    myTitle <- paste(myTitle , alts , samps , sep="\n")
    # Calculate maxylim
    maxylim=max(colSums(finalTab)) + round(max(colSums(finalTab))/3)
    if(maxylim<1){
        maxylim <- 1
    }
    # set cex for labels and numbers on top of bars
    mycextext <- if(ncol(finalTab)>20){
                                .5
                    }else if(ncol(finalTab)>10){
                                .8
                    } else {
                                1
                    }
    border <- if(nrow(finalTab)>50){
                    NA
                } else if(any(finalTab<3)) {
                    NA
                } else {
                    NULL
                }
    labels <- colnames(finalTab)
    mycolors <- if(nrow(finalTab)>1) c(myPalette(nrow(finalTab)-1) , "dark gray") else myPalette(1)
    # legend <- if(nrow(finalTab)>1) TRUE else FALSE
    bp <- barplot(finalTab , col=mycolors 
                # , legend=legend
                , border=border , xaxt="n" , xlab="" 
                , ylim=c(0 , maxylim)
                , main=myTitle)
    if(nrow(finalTab)>1){
        legendCols <- ceiling(nrow(finalTab)/5)
        legend("topright" 
               , legend = rownames(finalTab)
               , fill = mycolors
               , ncol = legendCols
               , cex = 0.6)
    }
    # xlabels 45 degress
    text(x =  bp, y = par("usr")[3] - maxylim/40, srt = 45, adj = 1
        ,labels = labels, xpd = TRUE , cex=mycextext-.1)
    # add percentage on top of bars
    perc <- paste0(100*(round(colSums(finalTab) , 3)) , "%")
    if(!is.na(grouping)) {
        perc[seq(1 , length(perc) , 2)] <- ""
    }
    text(x = bp+0.4 
        , y = colSums(finalTab)+(maxylim/20)
        , label = perc
        , pos = 3
        , cex = mycextext
        , col = "red" 
        , srt=90)
}

.tumorFreqPlotterStackHTML <- function(finalTab , var , grouping , alterationType , tumor.freqs){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    mycolors <- if(nrow(finalTab)>1) c(myPalette(nrow(finalTab)-1) , "#A9A9A9") else myPalette(1)
    myTitle <- paste("Number of covered samples:" , var , "BY" , grouping)
    alts <- paste("Data Type:" , paste(alterationType , collapse=", "))
    # tums <- unique(names(mysamples)[names(mysamples)!="all_tumors"])
    samps <- paste("Freqs:" 
                , paste( paste(names(tumor.freqs) , round(tumor.freqs , 3) , sep=" - ") , collapse=" , ")
                )
    myTitle <- paste(myTitle , alts , samps , sep="\n")
    finalTab2 <- as.data.frame(t(finalTab))
    finalTab2$Var <- rownames(finalTab2)
    htmlColChart <- gvisColumnChart(finalTab2 
                    , xvar = "Var"
                    , yvar = colnames(finalTab2)[colnames(finalTab2)!="Var"]
                    , options=list(
                           isStacked=TRUE
                           , explorer="{actions: ['dragToZoom','rightClickToReset'], maxZoomIn:0.05 , keepInBounds: true , axis: 'both'}"
                        , title=gsub("\n" , " - " , myTitle)
                        , height=800
                        ,vAxis="{title: 'Number of Covered Samples', 
                                   titleTextStyle: {color: '#000000'}}"
                        , hAxis.viewWindowMode="pretty"
                        , vAxis.viewWindowMode="pretty"
                        , chartArea= "{width: '90%', height: '90%'}"
                        , legend= "{position: 'in' , textStyle: {color: 'black', fontSize: 10} }"
                        , titlePosition= 'in'
                        , bar.groupWidth= '100%'
                        , tooltip.trigger="both"
                        , enableScrollWheel=TRUE
                        , colors=paste0("['" , paste(mycolors , collapse="','") , "']")
           ))
    return(htmlColChart)
}

setGeneric('coverageStackPlot', function(object
                                    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                    , var=c("drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , grouping=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , tumor_type=NULL
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=TRUE
                                    , tumor.weights=NULL
                                    , tumor.freqs=NULL
                                    , plotFreq = FALSE                                
                                    , noPlot=FALSE
                                    , html=FALSE) {
    standardGeneric('coverageStackPlot')
    })
setMethod('coverageStackPlot', 'CancerPanel', function(object
                                    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                    , var=c("drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , grouping=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , tumor_type=NULL
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=TRUE
                                    , tumor.weights=NULL
                                    , tumor.freqs=NULL
                                    , plotFreq = FALSE                            
                                    , noPlot=FALSE
                                    , html=FALSE)
{
    # Checks
    possibleAlterations <- c("copynumber" , "expression" , "mutations" , "fusions")
    possibleGrouping <- c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type" )
    possibleVar <- c("drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
    if(length(var)>1){
        var <- var[1]
    } else if(length(var)==0 || is.na(var)) {
        stop("var cannot be empty or NA")
    }
    if(length(grouping)>1){
        grouping <- grouping[1]
    } else if(length(grouping)==0) {
        stop("grouping cannot be empty")
    }
    if(var %eq% grouping){
        stop("var and grouping must differ")
    }
    if(var %notin% possibleVar){
        stop(paste("var can only be one or more of the following" , paste(possibleVar , collapse=", ")))
    }
    if(any(alterationType %notin% possibleAlterations)){
        stop(paste("alterationType can only be one or more of the following" 
                   , paste(possibleAlterations , collapse=", ")))
    }
    if(!any(is.na(grouping))){
        if(any(grouping %notin% possibleGrouping))
            stop(paste("grouping can only be one of the following:" 
                       , paste(possibleGrouping , collapse=", ")))
    }
    if(("alteration_id" %in% grouping) & length(alterationType)<2){
        stop("If you select 'alteration_id' as grouping variable, you must select more than one alterationType")
    }
    if(("alteration_id" %in% grouping) & collapseByGene){
        warning("If 'alteration_id' is in grouping variables, you cannot collapse by gene. The option was set to FALSE")
        collapseByGene <- FALSE
    }
    if(("alteration_id" %in% var) & collapseByGene){
        warning("If 'alteration_id' is in var, you cannot collapse by gene. The option was set to FALSE")
        collapseByGene <- FALSE
    }
    # Check tumor.weights consistency
    # tumor.freqs is a named vector of integers: e.g. c(brca=100 , luad=1000)
    if(!is.null(tumor.weights)){
        .tumor.weights.standardCheck(tumor.weights , tumor.freqs , object , tumor_type)
    }
    # Check tumor.freqs consistency
    # tumor.freqs is a named vector of 0-1 coefficient that sum to 1: e.g. c(brca=0.1 , luad=0.9)
    if(!is.null(tumor.freqs)){
        .tumor.freqs.standardCheck(tumor.weights , tumor.freqs , object , tumor_type)
        if("tumor_type" %in% var){
            stop("If you use tumor.freqs, tumor_type cannot be in var parameter")
        }
        if("tumor_type" %in% grouping){
            stop("If you use tumor.freqs, tumor_type cannot be in grouping parameter")
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
    
    varLevels <- switch(var 
                        ,gene_symbol=unique(c(object@arguments$panel$gene_symbol , mydata$gene_symbol))
                        ,drug=unique(c(object@arguments$panel$drug , mydata$drug))
                        ,group=unique(c(object@arguments$panel$group , mydata$group))
                        ,alteration_id=c("cna" , "expr" ,"fus" , "mut")
                        ,tumor_type=unique(mydata$tumor_type))
    if(!is.na(grouping)){
        groupingLevels <- switch(grouping 
                        ,gene_symbol=unique(c(object@arguments$panel$gene_symbol , mydata$gene_symbol))
                        ,drug=unique(c(object@arguments$panel$drug , mydata$drug))
                        ,group=unique(c(object@arguments$panel$group , mydata$group))
                        ,alteration_id=c("cna" , "expr" ,"fus" , "mut")
                        ,tumor_type=unique(mydata$tumor_type))
    }
    # Reduce the dataset to a specific amount of patients sampled at random
    # among the different tumor types
    # if(!is.null(tumor.weights)){
    #     newDataAndSamps <- .tumor.weights.machine(tumor.weights , mysamples , mydata )
    #     mydata <- newDataAndSamps[[1]]
    #     mysamples <- newDataAndSamps[[2]]
    # }
    # # Collapsing Mutation by gene
    # # This option is valid for mutation. If a gene is mutated in more than 1 spot in a patient
    # # Does it count for 1 or more alterations?
    # if(collapseMutationByGene) {
    #     mydata <- unique(mydata)
    # }
    # # Collapsing by gene
    # # This option is valid for all alterations.
    # # If a gene is altered in multiple ways (e.g. both amplified and mutated) on the same patient
    # # Does it count for 1 or more alterations?
    # if(collapseByGene) {
    #     mydata <- unique(mydata[ , colnames(mydata)!="alteration_id"])
    # }
    # If tumor.freqs is set, we basically run this method in recursive mode for each tumor type
    # and then aggregate everything
    if(!is.null(tumor.freqs)){
        covRecurse <- lapply(names(tumor.freqs) , function(tum){
            # browser()
                out <- tryCatch( suppressWarnings( coverageStackPlot(object
                    , alterationType=alterationType
                    , var=var
                    , grouping=grouping
                    , tumor_type=tum
                    , collapseMutationByGene=collapseMutationByGene
                    , collapseByGene=collapseByGene
                    , tumor.weights=NULL
                    , tumor.freqs=NULL                                
                    , noPlot=TRUE
                    , html=FALSE
                        )) , error = function(e) return(NULL))
                if(is.null(out)){
                    return(NULL)
                }
                out2 <- (out$plottedTable/out$Samples['all_tumors'])*tumor.freqs[tum]
                return(out2)
            })
        rows <- unique( lapply(covRecurse , rownames) %>% unlist)
        cols <- unique( lapply(covRecurse , colnames) %>% unlist)
        finalTab <- matrix(NA , nrow=length(rows) , ncol=length(cols))
        rownames(finalTab) <- rows
        colnames(finalTab) <- cols
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
                finalTab[i , j] <- sum(vals)
            }
        }
        if(noPlot){
            return( list(plottedTable=finalTab , Samples=NULL) )
        } else {
            if(!html){
                return(.tumorFreqPlotterStack(finalTab , var , grouping , alterationType , tumor.freqs))
            } else {
                return(.tumorFreqPlotterStackHTML(finalTab , var , grouping , alterationType , tumor.freqs))
            }
        }
    }
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    if(!is.na(grouping)){
        mytab <- table( factor(mydata[ , grouping] , levels=groupingLevels) 
                        , factor(mydata[ , var] , levels=sort(varLevels)) )
        myDataForTotal <- unique(mydata[ , c(var , "case_id")])
        myTotal <- table( factor(myDataForTotal[ , var] , levels=sort(varLevels)) )
        names(myTotal) <- paste("Total" , names(myTotal))
        for(i in seq_len(ncol(mytab))){
            if(i==1){
                finalTab <- cbind(c(mytab[ , i] , 0) 
                            ,c(rep(0 , length(mytab[ , i])), myTotal[i] )
                            )
            } else {
                finalTab <- cbind(finalTab , c(mytab[ , i] , 0))
                finalTab <- cbind(finalTab , c(rep(0 , length(mytab[ , i])), myTotal[i] ) )
            }
        }
    rownames(finalTab)[nrow(finalTab)] <- "Total"
    colnames(finalTab) <- lapply(seq_len(length(myTotal)) , function(x) {
                            c(colnames(mytab)[x] , names(myTotal)[x])
                            }) %>% unlist
    } else {
        myDataForTotal <- unique(mydata[ , c(var , "case_id")])
        finalTab <- table( factor(myDataForTotal[ , var] , levels=sort(varLevels)) ) %>% as.matrix %>% t
        rownames(finalTab) <- "noGrouping"
    }
    if(noPlot){
        return(list(plottedTable=finalTab , Samples=lengths(mysamples)))
    } else {
        if(!html){
            myTitle <- paste("Number of covered samples:" , var , "BY" , grouping)
            alts <- paste("Data Type:", paste(alterationType , collapse=", "))
            tums <- unique(names(mysamples)[names(mysamples)!="all_tumors"])
            samps <- paste("Samples:" , paste( paste0(tums , "=" , lengths(mysamples[tums])) , collapse=" "))
            myTitle <- paste(myTitle , alts , samps , sep="\n")
            # set cex for labels and numbers on top of bars
            mycextext <- if(ncol(finalTab)>20){
                                        .5
                            }else if(ncol(finalTab)>10){
                                        .8
                            } else {
                                        1
                            }
            border <- if(nrow(finalTab)>50){
                            NA
                        } else if(any(finalTab<3)) {
                            NA
                        } else {
                            NULL
                        }
            labels <- colnames(finalTab)
            mycolors <- if(nrow(finalTab)>1) c(myPalette(nrow(finalTab)-1) , "dark gray") else myPalette(1)
            # Calculate frequencies
            if(var!="tumor_type"){
                perc <- paste0(100*(round(colSums(finalTab)/length(mysamples$all_tumors) , 3)) , "%")
            } else {
                Samples <- lengths(mysamples)
                if(is.na(grouping)){
                    perc <- vapply(seq_len(ncol(finalTab)) , function(x) {
                                # browser()
                                tum <- colnames(finalTab)[x]
                                paste0(100*(round(sum(finalTab[ , x])/Samples[tum] , 3)) , "%")
                        } , character(1))
                } else {
                    perc <- vapply(seq_len(ncol(finalTab)) , function(x) {
                                tum <- sub("Total " , "" , colnames(finalTab)[x])
                                paste0(100*(round(sum(finalTab[ , x])/Samples[tum] , 3)) , "%")
                        } , character(1))
                }
            }
            # add percentage on top of bars
            # perc <- paste0(100*round(myPerc , 3) , "%")
            if(!is.na(grouping)) {
                perc[seq(1 , length(perc) , 2)] <- ""
            }
            # When plotFreq is TRUE, we plot the freqeuncies rather than absolute values
            if(plotFreq){
                if(var!="tumor_type"){
                    for(i in colnames(finalTab)){
                        finalTab[ , i] <- finalTab[ , i]/length(mysamples$all_tumors)
                    }
                } else {
                    Samples <- lengths(mysamples)
                    if(is.na(grouping)){
                        for(i in colnames(finalTab)){
                            finalTab[ , i] <- finalTab[ , i]/Samples[i]
                        }
                    } else {
                        for(i in colnames(finalTab)){
                            tum <- sub("Total " , "" , i) 
                            finalTab[ , i] <- finalTab[ , i]/Samples[tum]
                        }
                    }
                }
            }
            maxylim <- max(colSums(finalTab)) + round(max(colSums(finalTab))/3)
            if(plotFreq){
                maxylim <- ceiling(maxylim)
            } else {
                if(maxylim==1){
                    maxylim <- 2
                }
                if(maxylim==0){
                    maxylim <- 1
                }
            }

            bp <- barplot(finalTab , col=mycolors 
                        # , legend=legend
                        , border=border , xaxt="n" , xlab="" 
                        , ylim=c(0 , maxylim)
                        , main=myTitle)
            if(nrow(finalTab)>1){
                legendCols <- ceiling(nrow(finalTab)/5)
                legend("topright" 
                       , legend = rownames(finalTab)
                       , fill = mycolors
                       , ncol = legendCols
                       , cex = 0.6)
            }
            # xlabels 45 degress
            text(x =  bp, y = par("usr")[3] - maxylim/40, srt = 45, adj = 1
                ,labels = labels, xpd = TRUE , cex=mycextext-.1)
            text(x = bp+0.4 
                , y = colSums(finalTab)+(maxylim/20)
                , label = perc
                , pos = 3
                , cex = mycextext
                , col = "red" 
                , srt=90)
        } else {
            mycolors <- if(nrow(finalTab)>1) c(myPalette(nrow(finalTab)-1) , "#A9A9A9") else myPalette(1)
            myTitle <- paste("Number of covered samples:" , var , "BY" , grouping)
            alts <- paste("Data Type:" , paste(alterationType , collapse=", "))
            tums <- unique(names(mysamples)[names(mysamples)!="all_tumors"])
            samps <- paste("Samples:" , paste( paste0(tums , "=" , lengths(mysamples[tums])) , collapse=" "))
            myTitle <- paste(myTitle , alts , samps , sep="\n")
            finalTab2 <- as.data.frame(t(finalTab))
            sampLen <- length(mysamples$all_tumors)
            if(var!="tumor_type"){
                finalTab2_tooltip <- lapply(colnames(finalTab2) , function(x) {
                            mynum <- round((finalTab2[ , x]/sampLen)*100 , 2)
                            out <- paste(x , paste0(mynum , "%") , sep=": ")
                            out <- paste(out , paste("Number of alterations:" , finalTab2[ , x]) , sep="\\n\\n")
                            out <- paste(out , paste("Reference Sample size:" , sampLen) , sep="\\n\\n")
                            return(out)
                            })
            } else {
                # If the variable is tumor type we must divide for the number of samples of the tumor type, not the total
                # In case there is no grouping variable, the "Total" column doesn't exist
                noTotalRowNames <- grep("^Total" , rownames(finalTab2) , invert=TRUE , value=TRUE)
                if(nrow(finalTab2)==length(noTotalRowNames)){
                    sampTum <- lengths(mysamples)[noTotalRowNames]
                    finalTab2_tooltip <- lapply(colnames(finalTab2) , function(x) {
                            mynum <- round((finalTab2[ , x]/sampTum)*100 , 2)
                            out <- paste(x , paste0(mynum , "%") , sep=": ")
                            out <- paste(out , paste("Number of alterations:" , finalTab2[ , x]) , sep="\\n\\n")
                            out <- paste(out , paste("Reference Sample size:" , sampTum) , sep="\\n\\n")
                            return(out)
                            })
                    names(finalTab2_tooltip) <- colnames(finalTab2)
                } else {
                    sampTum <- lengths(mysamples)[noTotalRowNames]
                    sampLenTum <- rep(sampTum , rep(2 , length(sampTum)))
                    finalTab2_tooltip <- lapply(colnames(finalTab2) , function(x) {
                                mynum <- round((finalTab2[ , x]/sampLenTum)*100 , 2)
                                out <- paste(x , paste0(mynum , "%") , sep=": ")
                                out <- paste(out , paste("Number of alterations:" , finalTab2[ , x]) , sep="\\n\\n")
                                out <- paste(out , paste("Reference Sample size:" , sampLenTum) , sep="\\n\\n")
                                return(out)
                                })
                    names(finalTab2_tooltip) <- colnames(finalTab2)
                }
            }
            finalTab3 <- data.frame(Var=row.names(finalTab2))
            for(i in seq_len(ncol(finalTab2))){
                nameofthecol <- colnames(finalTab2)[i]
                finalTab3[ , nameofthecol] <- finalTab2[ , nameofthecol]
                finalTab3[ , paste(nameofthecol , ".html.tooltip" , sep="")] <- paste(paste("Column:" , row.names(finalTab2)) 
                                                                                    , finalTab2_tooltip[[i]] , sep="\\n\\n")
            }
            htmlColChart <- gvisColumnChart(finalTab3 
                            , xvar = "Var"
                            , yvar = colnames(finalTab3)[colnames(finalTab3)!="Var"]
                            , options=list(
                        isStacked=TRUE
                        , explorer="{actions: ['dragToZoom','rightClickToReset'], maxZoomIn:0.05 , keepInBounds: true , axis: 'both'}"
                        , title=gsub("\n" , " - " , myTitle)
                        , height=800
                        ,vAxis="{title: 'Number of Covered Samples', 
                                   titleTextStyle: {color: '#000000'}}"
                        , hAxis.viewWindowMode="pretty"
                        , vAxis.viewWindowMode="pretty"
                        , chartArea= "{width: '90%', height: '90%'}"
                        , legend= "{position: 'in' , textStyle: {color: 'black', fontSize: 10} }"
                        , titlePosition= 'in'
                        , bar.groupWidth= '100%'
                        , tooltip.trigger="both"
                        , enableScrollWheel=TRUE
                        , colors=paste0("['" , paste(mycolors , collapse="','") , "']")
                ))
          return(htmlColChart)
      }
    }
})





