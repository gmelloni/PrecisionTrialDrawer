.saturationPlotter <- function(mydata_melt 
                                , grouping 
                                , adding 
                                , y_measure
                                , main=""
                                , legend 
                                , labelling 
                                , annotation=FALSE 
                                , mylabel = NULL 
                                , CI.show=FALSE) {
    labels <- c( levels(mydata_melt$grouping) )
    # Full global assignment of variable for R CMD check complainings
    Space=Mean=CI=Coverage=NULL
    if(y_measure[1]=="mean") {
        plotTitle <- main
        p <- ggplot(mydata_melt, aes(x=Space, y=Mean, colour=grouping, group=grouping))
        if(CI.show){
            p <- p + geom_errorbar(aes(ymin=Mean-CI, ymax=Mean+CI), colour="black", width=.1)
        }
    } else {
        plotTitle <- main
        p <- ggplot(mydata_melt, aes(x=Space, y=Coverage, colour=grouping, group=grouping)) +
            scale_y_continuous(limits=c(0,1))
    }
    # Thanks to MrFlick at stackoverflow for this function
    # We are now able to combine different aestethics: aes + aes_string
    '+.uneval' <- function(a,b) {
      'class<-'(modifyList(a,b), "uneval")
    }
    p <- p +
        geom_line() +
        geom_point(size=3 ) + 
        xlab("Genomic Space in kBase") +
        # ylab("Mean Number of variants per sample") +
        theme_bw() +
        scale_colour_hue(name=grouping,
                    breaks=levels(mydata_melt$grouping),
                    labels=labels,
                    l=40) +
        ggtitle(plotTitle)
    #   Original Version with overlapping labels
    #        geom_text(aes_string(label=adding
    #                    , hjust=.5
    #                    , vjust=-.1
    #                    , angle=45
    #                    ,show=FALSE
    #                    ,size=.5),show.legend = FALSE) +
    #        geom_label_repel(aes_string(label = adding) + aes(fill=grouping)
    #                        #,fontface = 'bold'
    #                        ,color = 'white'
    #                        ,box.padding = unit(0.6, "lines")
    #                        ,point.padding = unit(0.5, "lines")
    #                        ,force=2
    #                        ,show.legend=FALSE) +
    if(annotation){
        p <- p + annotate("text", x = Inf, y = -Inf, label = paste(mylabel , collapse="\n")
        , hjust=1.1
        , vjust=-.2
        , col="black", cex=2.5,fontface = "bold")
    }
    if(legend=="in"){
        p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))
    }
    if(labelling){
        p <- p + geom_text_repel(aes_string(label=adding) , box.padding = unit(0.6, "lines"))
    }
    return(p)
}

# scatter plot of panel size by mean coverage or by absolute coverage
setGeneric('saturationPlot', function(object
                                    , alterationType=c("copynumber" , "expression" , "mutations", "fusions")
                                    , grouping=c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
                                    , adding=c("alteration_id" , "gene_symbol" , "drug" , "group")
                                    , tumor_type=NULL
                                    , y_measure=c("mean" , "absolute")
                                    , adding.order=c("absolute" , "rate")
                                    , sum.all.feature=FALSE
                                    #, sum.by.grouping=FALSE
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=FALSE
                                    , labelling=TRUE
                                    , tumor.weights=NULL
                                    , main=""
                                    , legend=c("in" , "out")
                                    , noPlot=FALSE) {
    standardGeneric('saturationPlot')
    })
setMethod('saturationPlot', 'CancerPanel', function(object
                                    , alterationType=c("copynumber" , "expression" , "mutations", "fusions")
                                    , grouping=c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
                                    , adding=c("alteration_id" , "gene_symbol" , "drug" , "group")
                                    , tumor_type=NULL
                                    , y_measure=c("mean" , "absolute")
                                    , adding.order=c("absolute" , "rate")
                                    , sum.all.feature=FALSE
                                    #, sum.by.grouping=FALSE
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=FALSE                                
                                    , labelling=TRUE
                                    , tumor.weights=NULL
                                    , main=""
                                    , legend=c("in" , "out")
                                    , noPlot=FALSE)
{
    ## Check of parameters
    possibleAddingOrder <- c("absolute" , "rate")
    possibleAlterations <- c("copynumber" , "expression" , "mutations" , "fusions")
    possibleGrouping <- c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
    possibleAdding <- c("alteration_id" , "gene_symbol" , "drug" , "group")
    possibley_measure <- c("mean" , "absolute")
    if(any(is.na(grouping))){
        grouping <- NA
    }
    if(length(grouping)!=1)
        grouping <- grouping[1]
    if(length(adding)!=1)
        adding <- adding[1]
    if(length(adding.order)!=1)
        adding.order <- adding.order[1]
    if(any(adding.order %notin% possibleAddingOrder)){
        stop("adding.order can only be one or more of the following" %++% paste(possibleAddingOrder , collapse=", "))
    }
    if(any(alterationType %notin% possibleAlterations)){
        stop("alterationType can only be one or more of the following" %++% paste(possibleAlterations , collapse=", "))
    }
    if(any(grouping %notin% possibleGrouping)){
        stop("grouping can only be one of the following" %++% paste(possibleGrouping , collapse=", "))
    }
    if(any(adding %notin% possibleAdding)){
        stop("adding can only be one of the following" %++% paste(possibleAdding , collapse=", "))
    }
    if(any(y_measure %notin% possibley_measure)){
        stop("y_measure can only be one of the following" %++% paste(possibley_measure , collapse=", "))
    }
    if(("alteration_id" %in% grouping) & length(alterationType)<2){
        stop("If you select 'alteration_id' as grouping variable, you must select more than one alterationType")
    }
    if(!is.logical(sum.all.feature)){
      stop("sum.all.feature must be TRUE or FALSE")
    }
    # if(!is.logical(sum.by.grouping)){
    #   stop("sum.by.grouping must be TRUE or FALSE")
    # }
    legend <- legend[1]
    if(legend %notin% c("in" , "out")){
        stop("legend can only be 'in' or 'out' of the plotting area")
    }
    if(!is.null(tumor_type)){
        if(!all(tumor_type %in% object@arguments$tumor_type)){
            stop("You selected a tumor_type that has no data in this CancerPanel object")
        }
    }
    # Check tumor.weights consistency
    if(!is.null(tumor.weights)){
        if( any(is.na(tumor.weights)) || any(!is.numeric(tumor.weights)) ){
            stop("tumor.weights must be an integer vector")
        }
        if( any(names(tumor.weights) %notin% object@arguments$tumor_type) ){
            stop("The tumor.weights names do not correspond to the tumor types present in this object")
        }
        if(grouping %in% "tumor_type"){
            warning("If you group by tumor type, tumor.weights are ignored. Every tumor has its own plot")
            tumor.weights <- NULL
        } else {
            if(!is.null(tumor_type)){
                if( any(names(tumor.weights) %notin% tumor_type) ){
                    stop("Weights and tumor_type selected do not match")
                }
            }
        }
    }

    #----------------------------
    # GRAB DATA AND SAMPLES
    #----------------------------

    myenv <- new.env()
    dataExtractor(object=object , alterationType=alterationType , tumor_type=tumor_type 
            , collapseMutationByGene=collapseMutationByGene , collapseByGene=collapseByGene 
            , myenv=myenv , tumor.weights=tumor.weights)
    mydata <- get("mydata" , envir=myenv)
    mysamples <- get("mysamples" , envir=myenv)
    tum_type_diff <- get("tum_type_diff" , envir=myenv)
    # rm(list=ls(myenv) , envir=myenv)
    rm(myenv)
    # First, split by grouping variable
    if(is.na(grouping)){
        mydata_split <- list("NA"=mydata)
    } else {
        mydata_split <- split(mydata , mydata[[grouping]])
    }
    mydata_split <- lapply(mydata_split , function(x) {
                        # browser()
                        tums <- unique(x$tumor_type)
                        x$case_id <- factor(x$case_id , levels=unlist(mysamples[tums]))
                        return(x)
        })
    # Aggregate the panel according to the adding variable
    panel_transition <- object@arguments$panel
    # panel_transition <- panel
    panel_transition$alteration_id <- .mapvalues(from=c("SNV" , "CNA" , "fusion" , "expression")
                                                        ,to=c("mut" , "cna" , "fus" , "expr")
                                                        ,panel_transition$alteration
                                                        ,warn_missing=FALSE)
    if(is.na(grouping)){
        panel_transition_split <- list("NA"=panel_transition)
    } else if(grouping=="tumor_type"){
        panel_transition_split <- lapply(names(mysamples) , function(x) panel_transition)
        names(panel_transition_split) <- names(mysamples)
    } else {
        panel_transition_split <- split(panel_transition , panel_transition[[grouping]])
    }
    # Here we calculate space. We developed 4 possible scenarios of space calculation:
    # 1) sum.by.grouping=FALSE the space is the same across grouping variable
    # 2) sum.by.grouping=TRUE the space is calculate by splitting the panel by grouping variable
    # 3) sum.all.feature=TRUE every row of the panel is summed up by adding variable
    # 4) sum.all.feature=FALSE every row of the panel is summed up by adding variable and gene

    # New version 07/10/2017
    # sum by grouping makes no sense. If I put a grouping variable I want it to be divided in as many strata as grouping has
    # sum.all.feature should stay
    # sum.all.feature=TRUE. In case of a panel of cna and snv, a gene is counted twice
    # sum.all.feature=FALSE same gene counts 1

    # browser()
    # if(sum.by.grouping){
    #     panel_agg <- lapply(panel_transition_split , function(df) {
    #       if(sum.all.feature){
    #         return(aggregate(as.formula(paste("variation_len~" , adding)) , df , sum))
    #       } else {
    #         df$full <- ifelse(df$alteration %in% c("CNA" , "fusion" , "expression") , "yes" , 
    #                           ifelse(df$exact_alteration=="" , "yes" , "no"))
    #         genefull <- unique(df$gene_symbol[df$full=="yes"])
    #         df_small <- df[ , unique(c(adding , "gene_symbol" , "variation_len" , "full"))]
    #         df_small2 <- df_small[ !(df_small$gene_symbol %in% genefull & df_small$full=="no") , ]
    #         return(aggregate(as.formula(paste("variation_len~" , adding)) , df_small2 , sum))
    #       }
    #     })
    # } else {
        .spacing <- function(df , sum.all.feature){
            if(sum.all.feature){
                   out <- aggregate(as.formula(paste("variation_len~" , adding)) , df , sum)
                # return(lapply(names(panel_transition_split) , function(x) out))
             } else {
                   df$full <- ifelse(df$alteration %in% c("CNA" , "fusion" , "expression") , "yes" , 
                          ifelse(df$exact_alteration=="" , "yes" , "no"))
                genefull <- unique(df$gene_symbol[df$full=="yes"])
                df_small <- df[ , unique(c(adding , "gene_symbol" , "variation_len" , "full"))]
                # If I ask for a full gene in CNA and some mutations, remove the mutations when counting length
                df_small2 <- df_small[ !(df_small$gene_symbol %in% genefull & df_small$full=="no") , ]
                df_small2 <- unique(df_small2)
                out <- aggregate(as.formula(paste("variation_len~" , adding)) , df_small2 , sum)
               # return(lapply(names(panel_transition_split) , function(x) out))
            }
            return(out)
        }
        # panel_agg <- .spacing(panel_transition , sum.all.feature)
        panel_agg <- lapply( panel_transition_split , function(x) .spacing(x , sum.all.feature))
    # }
    names(panel_agg) <- names(panel_transition_split)
    # Cast the splitted data by the adding variable
    caster <- lapply(1:length(mydata_split) , function(i) {
        # browser()
        z <- names(mydata_split)[i]
        x <- mydata_split[[i]]
        x_cast <- reshape2::dcast(formula=as.formula( paste("case_id~" , adding)) 
                        , data=x 
                        , fun.aggregate=length 
                        , drop=FALSE 
                        , value.var="case_id")
        rownames(x_cast) <- as.character(x_cast$case_id)
        x_cast$case_id <- NULL
        x_cast[ is.na(x_cast) ] <- 0
        x_cast[ x_cast>1 ] <- 1
        # assign length to the adding variable (columns)
        x_cast_colnames <- colnames(x_cast)
        # Problem. The space can be NA if we consider translocation in every gene of a panel
        # The panel doesn't know what translocation could appear
        # In the example, RET was included for every translocation, but ERC1__RET is not in the panel
        if(adding=="gene_symbol" & any(panel_transition[ , "exact_alteration"]=="gene_fusion")){
            gene_fusions <- unique(panel_transition[ panel_transition$exact_alteration=="gene_fusion", "gene_symbol"])
            subVec <- paste0(gene_fusions , "__")
            subVec2 <- paste0("__",gene_fusions)
            searchPattern <- paste( c(subVec2 , subVec) , collapse="|")
            # df$gene_symbol2 <- sub(searchPattern , "" , df$gene_symbol)
            x_cast_colnames2 <- str_extract(x_cast_colnames , searchPattern) %>% sub("__" , "" , .)
            x_cast_colnames2[is.na(x_cast_colnames2)] <- x_cast_colnames[is.na(x_cast_colnames2)]
            space <- panel_agg[[z]][ match(x_cast_colnames2 , panel_agg[[z]][,adding]), 'variation_len']
            space <- ifelse(is.na(space) , min(space , na.rm=TRUE) , space)
            space <- space/1000
            len_df <- data.frame(gene_symbol2=x_cast_colnames2 
                                , gene_symbol=x_cast_colnames 
                                , Space=space 
                                , stringsAsFactors=FALSE)
            if(!grouping %in% colnames(panel_transition)){
                total <- unique(panel_transition[ , adding])
                missing <- setdiff( total , len_df[ , "gene_symbol2"] )
            }else{
                total <- unique(panel_transition[panel_transition[ , grouping]==z , adding])
                missing <- setdiff( total , len_df[ , "gene_symbol2"] )
            }
            # df$gene_symbol2 <- NULL
        } else {
            space <- panel_agg[[z]][ match(x_cast_colnames , panel_agg[[z]][,adding]), 'variation_len']
            space <- ifelse(is.na(space) , min(space , na.rm=TRUE) , space)
            space <- space/1000
            len_df <- data.frame(addingVar=x_cast_colnames , Space=space , stringsAsFactors=FALSE)
            colnames(len_df) <- c(adding , "Space")
            if(!grouping %in% colnames(panel_transition)){
                total <- unique(panel_transition[ , adding])
                missing <- setdiff( total , len_df[ , adding] )
            }else{
                total <- unique(panel_transition[panel_transition[ , grouping]==z , adding])
                missing <- setdiff( total , len_df[ , adding] )
            }
        }
        # calculate number of variations per adding variable
        adding_sum <- colSums(x_cast) %>% sort %>% rev
        # If the order is based just on number of variants:
        if(adding.order=="absolute"){
            x_cast <- x_cast[ , names(adding_sum) , drop=FALSE]
            df <- data.frame(adding_var=names(adding_sum) , stringsAsFactors=FALSE)
        } else {
            # If the order is based on rate of variants per length
            # things are a little more complicated
            # create a rate dafarame first and order it by rate
            adding_sum_df <- data.frame(addingVar=names(adding_sum) 
                                        , var_num=unname(adding_sum) 
                                        , stringsAsFactors=TRUE)
            colnames(adding_sum_df) <- c(adding , "var_num")
            len_df <- merge(len_df , adding_sum_df , all.y=TRUE)
            len_df$num_of_variants_per_KB <- (len_df$var_num/len_df$Space)
            len_df$var_num <- NULL
            len_df <- len_df[ order(len_df$num_of_variants_per_KB , decreasing=TRUE) , ]
            x_cast <- x_cast[ , len_df[ , adding] , drop=FALSE]
            df <- data.frame(adding_var=len_df[ , adding] , stringsAsFactors=FALSE)
        }
        x_cast2 <- sapply(1:ncol(x_cast), function(i) rowSums(x_cast[,1:i,drop=FALSE]))
        # create the dataset to plot
        df <- rename(df, c("adding_var" = adding))
        df$grouping <- z
        df$Mean <- apply(x_cast2 , 2 , mean)
        df$Coverage <- apply(x_cast2 , 2 , function(x) length(x[x!=0])/length(x))
        df$SD <- apply(x_cast2 , 2 , sd)
        df$SE <- df$SD / sqrt(nrow(x_cast2))  # Calculate standard error of the mean
        confInt <- .95
        ciMult <- qt(confInt/2 + .5, nrow(x_cast2)-1)
        df$CI <- df$SE * ciMult
        df <- .mergeOrder(df , len_df , by=adding , all.x=TRUE)
        df$Space <- cumsum(df$Space)
        return(list(df , c(length(missing) , length(total) )))
    })
    # mydata_melt <- do.call("rbind" , lapply(caster , '[[' , 1))
    mydata_melt <- as.data.frame(data.table::rbindlist(lapply(caster , '[[' , 1)) , stringsAsFactors=FALSE)
    mydata_melt$grouping <- factor(mydata_melt$grouping , levels=unique(mydata_melt$grouping))
    
    labels <- c( levels(mydata_melt$grouping) )
    mylabel <- sapply(1:length(caster) , function(x) {
                    "Missing" %++% adding %++% "in '" %++% 
                    names(mydata_split)[x] %+% "' :" %++% 
                    caster[[x]][[2]][1] %+% "/" %+% caster[[x]][[2]][2]
                    })
    if(noPlot){
        return(mydata_melt)
    } else {
        .saturationPlotter(mydata_melt 
                            , grouping 
                            , adding 
                            , y_measure 
                            , main                            
                            , legend
                            , labelling 
                            , annotation=TRUE 
                            , mylabel = mylabel 
                            , CI.show=TRUE)
    }
})