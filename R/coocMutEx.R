# Multiple ggplot in one figure
.multiplot <- function(plotlist=NULL, plot_layout=NULL) {
    # Modified from R cookbook
    # require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- plotlist

    numPlots = length(plots)

    # Make the panel
    plotCols = plot_layout[2]                          # Number of columns of plots
    plotRows = plot_layout[1] # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow,  curCol))
    }
  }

# Given a casted matrix, calculated pvalue for cooccurence
# and pvalue for mutual exclusivity with montercarlo or fisher test
# mutuated from package cooccur
.coocMutEx <- function (mat , thresh = FALSE , prob = c("hyper" , "firth")) {
    prob <- prob[1]
    spp_site_mat <- t(mat)
    spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
    spp_site_mat[spp_site_mat > 0] <- 1
    nsite <- ncol(spp_site_mat)
    nspp <- nrow(spp_site_mat)
    spp_pairs <- choose(nspp, 2)
    obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs, ncol = 3)
    prob_occur <- incidence <- cbind(1:nspp, rowSums(spp_site_mat, na.rm = TRUE))
    prob_occur[ , 2] <- prob_occur[ , 2]/nsite
    spp_pairs_mat <- t(combn(nspp,2))
    pairs_observations <- apply(spp_pairs_mat, 1, function(x) {
        spp <- x[1]
        spp_next <- x[2]
        pairs <- sum( spp_site_mat[spp, ] * spp_site_mat[spp_next, ] )
        pairs_prob <- sum( prob_occur[spp, 2] * prob_occur[spp_next, 2] )
        pairs_exp <- pairs_prob * nsite
        return(c(pairs, pairs_prob, pairs_exp))
        })

    obs_cooccur <- cbind(spp_pairs_mat, pairs_observations[1,])
    prob_cooccur <- cbind(spp_pairs_mat, pairs_observations[2,])
    exp_cooccur <- cbind(spp_pairs_mat, pairs_observations[3,])
    if (thresh) {
        # n_pairs <- nrow(prob_cooccur)
        prob_cooccur <- prob_cooccur[exp_cooccur[, 3] >= 1, ,drop=FALSE]
        obs_cooccur <- obs_cooccur[exp_cooccur[, 3] >= 1, ,drop=FALSE]
        exp_cooccur <- exp_cooccur[exp_cooccur[, 3] >= 1, ,drop=FALSE]
        # n_omitted <- n_pairs - nrow(prob_cooccur)
    }
    output <- t(sapply(1:nrow(obs_cooccur) , function(row) {
        # browser()
        sp1 <- obs_cooccur[row, 1]
        sp2 <- obs_cooccur[row, 2]
        # sp1_inc <- incidence[incidence[, 1] == sp1, 2]
        # sp2_inc <- incidence[incidence[, 1] == sp2, 2]
        # max_inc <- max(sp1_inc, sp2_inc)
        # min_inc <- min(sp1_inc, sp2_inc)
        prob_share_site <- rep(0, (nsite + 1))
        if (prob == "hyper") {
            # all.probs <- phyper(0:min_inc, min_inc, nsite - min_inc, max_inc)
            # prob_share_site[1] <- all.probs[1]
            # prob_share_site[2:length(all.probs)] <- diff(all.probs)
            p_lt <- table(mat[ , sp1] , mat[, sp2]) %>% 
                    fisher.test(. , alternative = "less")
            p_lt <- p_lt$p.value                    
            p_gt <- table(mat[ , sp1] , mat[, sp2]) %>% 
                    fisher.test(. , alternative = "greater")
            p_gt <- p_gt$p.value
        }
        # if (prob == "comb") {
        #     coprob <- function(max_inc,j,min_inc,nsite){
        #         as.numeric(round(choose(max_inc,j) * 
        #           choose(nsite - max_inc, min_inc - j),0) / 
        #             round(choose(nsite,min_inc),0))
        #     }
        #     ix1 <- (sp1_inc + sp2_inc) <= (nsite + 0:nsite) 
        #     ix2 <- 0:nsite <= min_inc
        #     prob_share_site[ix1 & ix2] <- coprob(max_inc = max_inc, 
        #             j = (0:nsite)[ix1 & ix2], min_inc = min_inc, nsite = nsite)
        #     ix <- 0:nsite <= obs_cooccur[row, 3]
        #     p_lt <- sum(prob_share_site[ix])
        #     ix <- (0:nsite >= obs_cooccur[row, 3])
        #     p_gt <- sum(prob_share_site[ix])
        # }
        if( prob == "firth"){
          firth <- tryCatch({
                    brglm::brglm( mat[ , sp1] ~ mat[, sp2] 
                        , family = binomial("logit") 
                        , method="brglm.fit")
                    } , warning=function(w) {
                          warning("Penalization couldn't reach convergence, using regular glm")
                          firth <- brglm::brglm( mat[ , sp1] ~ mat[, sp2] 
                            , family = binomial("logit") 
                            , method="glm.fit")
                          return(firth)
                    })
          firthSum <- brglm::summary.brglm(firth)
          firthP <- firthSum$coefficients[2 , "Pr(>|z|)"]
          if(firthSum$coefficients[2 , "z value"]>=0){
            p_gt <- firthP
            p_lt <- 1
          } else {
            p_lt <- firthP
            p_gt <- 1
          }
        }
        # p_exactly_obs <- prob_share_site[max(which(ix))]
        # p_lt <- round(p_lt, 5)
        # p_gt <- round(p_gt, 5)
        # p_exactly_obs <- round(p_exactly_obs, 5)
        # prob_cooccur[row, 3] <- round(prob_cooccur[row, 3], 3)
        # exp_cooccur[row, 3] <- round(exp_cooccur[row, 3], 1)
        # return(c(sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row, 
        #     3], prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt, 
        #     p_gt))
        # return(c(sp1, sp2, sp1_inc, sp2_inc, p_lt, p_gt))
        return(c(sp1, sp2, p_lt, p_gt))
        }))
    output <- as.data.frame(output)
    # colnames(output) <- c("sp1", "sp2", "sp1_inc", "sp2_inc", 
    #     "obs_cooccur", "prob_cooccur", "exp_cooccur", "pVal.MutEx", 
    #     "pVal.Cooc")
    # colnames(output) <- c("sp1", "sp2", "sp1_inc", "sp2_inc", "pVal.MutEx", "pVal.Cooc")
    colnames(output) <- c("sp1", "sp2", "pVal.MutEx", "pVal.Cooc")
    # Substitue NaNs and NA with 1
    output$pVal.Cooc <- ifelse( is.nan(output$pVal.Cooc) , 1 , output$pVal.Cooc)
    output$pVal.MutEx <- ifelse( is.nan(output$pVal.MutEx) , 1 , output$pVal.MutEx)
    # Substitute 0s with a very very small pvalue
    output$pVal.Cooc <- ifelse( output$pVal.Cooc==0 , 10^-100 , output$pVal.Cooc)
    output$pVal.MutEx <- ifelse( output$pVal.MutEx==0 , 10^-100 , output$pVal.MutEx)
    sp1_name <- merge(x = data.frame(order = 1:length(output$sp1), 
            sp1 = output$sp1), y = spp_key, by.x = "sp1", by.y = "num", 
            all.x = TRUE, sort = FALSE)
    sp2_name <- merge(x = data.frame(order = 1:length(output$sp2), 
            sp2 = output$sp2), y = spp_key, by.x = "sp2", by.y = "num", 
            all.x = TRUE, sort = FALSE)
    output$sp1_name <- sp1_name[with(sp1_name, order(order)), 
            "spp"]
    output$sp2_name <- sp2_name[with(sp2_name, order(order)), 
            "spp"]
    n_samp <- nrow(mat)
    TableOR_corr <- function(x , y){
          xtab <- table(x , y)+.5
          n00 <- xtab[1,1]
          n01 <- xtab[1,2]
          n10 <- xtab[2,1]
          n11 <- xtab[2,2]
          OR <- (n00 * n11)/(n01 * n10)
          return(OR)
        }
    stats <- t(apply(output , 1 , function(x) {
                    # browser()
                    pair <- c(x["sp1_name"] , x["sp2_name"])
                    mat_pair <- mat[ , pair]
                    freq_1 <- sum(mat_pair[ , 1])/n_samp
                    freq_2 <- sum(mat_pair[ , 2])/n_samp
                    combFreq <- mat_pair[ , 1] + mat_pair[ , 2]
                    combFreq[combFreq>1] <- 1
                    combFreq <- sum(combFreq)/n_samp
                    or <- TableOR_corr(mat_pair[ , 1] , mat_pair[ , 2])
                    c(or , freq_1 , freq_2 , combFreq)
                }))
    colnames(stats) <- c("OR" 
                        , "freq_1" 
                        , "freq_2" 
                        , "combFreq")
    finalOutput <- cbind(output , stats)
    # finalOutput <- finalOutput[ , c("sp1_name" , "sp2_name" 
    #                               , "pVal.MutEx", "pVal.Cooc"
    #                               ,"OR" , "freq_1" , "freq_2" , "combFreq")]
    finalOutput <- finalOutput[ , c("sp1_name" , "sp2_name" 
                                  , "pVal.MutEx", "pVal.Cooc"
                                  ,"OR")]
    finalOutput <- .changeFactor(finalOutput)
    return(finalOutput)
}

# When data are all not significant or there are no data at all,
# Plot a cool empty plot with title and text
.emptyGGplot <- function(label , title){
  df <- data.frame()
  empty <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) + xlab("") + ylab("") + 
            theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank()) +
            annotate("text"
                    , label = label
                    , x = 5, y = 5, size = 4, colour = "black") +
            ggtitle(title) +
            theme(plot.margin=unit(c(1,1,1,1), "cm"),legend.justification=c(1,0), legend.position=c(1,0))
  return(empty)
}

# cool upper matrix plot, mutuated from package cooccur
.coocMutEx_melter <- function(co_tab, plotrandom , pvalthr){
  if(is.null(co_tab)){
    return(NULL)
  }
  all_elements <- unique(c(co_tab$sp1_name , co_tab$sp2_name))
  mat <- matrix(0 , nrow=length(all_elements) , ncol=length(all_elements))
  dimnames(mat) <- list(all_elements, all_elements)
  co_tab$pVal.MutEx <- -co_tab$pVal.MutEx
  for(i in 1:nrow(co_tab)){
    x <- co_tab[i , ,drop=TRUE]
    sp1 <- x[["sp1_name"]]
    sp2 <- x[["sp2_name"]]
    if(abs(x[["pVal.Cooc"]])<abs(x[["pVal.MutEx"]])) {
      if(x[["pVal.Cooc"]]<=pvalthr){
          mat[sp1,sp2] <- -log10(x[["pVal.Cooc"]])
          mat[sp2,sp1] <- -log10(x[["pVal.Cooc"]])
      }
    } else {
      if(abs(x[["pVal.MutEx"]])<=pvalthr){
          mat[sp1,sp2] <- log10(abs(x[["pVal.MutEx"]]))
          mat[sp2,sp1] <- log10(abs(x[["pVal.MutEx"]]))
      }
    }
  }
  rmrandom <- function(mat){
    rowidx <- apply(mat, 1 , function(x) all(x==0))
    colidx <- apply(mat, 2 , function(x) all(x==0))
    idx <- rowidx & colidx
    return(mat[!idx , !idx])
  }
  if(plotrandom){
    cormat <- mat
  } else {
    cormat <- rmrandom(mat)
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
  #if(all(melted_cormat$value==0)){
  #    return(NULL
        # .emptyGGplot(label= "No Significant pair"
        #                       , title=title
        #                       )
  #    )
  #}
  melted_cormat$value <- ifelse(melted_cormat$value==0 , NA , melted_cormat$value)
  return(melted_cormat)
}

.plot_coocMutEx <- function(melted_cormat , title , minMelt , maxMelt , pvalthr){
  # minMelt <- min(melted_cormat$value , na.rm=TRUE)
  # maxMelt <- max(melted_cormat$value , na.rm=TRUE)
  log10pvalthr <- - log10(pvalthr)
  # fakedf <- data.frame(Random=NA)
  # Avoid complaints from check
  Var1=Var2=value=NULL
  ggheatmap <- ggplot(melted_cormat, aes(Var1, Var2)) + 
              # geom_tile(aes(fill = factor(value,levels=c(-1,0,1))), colour ="white")+
              # scale_fill_manual(values = c("#FFCC66","dark gray","light blue")
              #               , name = ""
              #               , labels = c("MutEx","random","Cooc")
              #               ,drop=FALSE) +
              geom_tile(aes(fill = value), colour ="white" )+
              # scale_fill_discrete(name="p-value\nyellow cooc\nblue mutex" , na.value="dark gray" , labels = )
              scale_fill_gradient2(space="Lab"
                                ,low = "navy blue"
                                ,na.value="dark gray"
                                ,high="#FFCC66"
                                ,midpoint=0
                                ,mid="dark gray"
                                # ,name="log10 Pvalue\nNegative MutEx\nPositive Cooc"
                                ,name="p-value\nyellow cooc\nblue mutex"
                                ,breaks=c(minMelt
                                          # ,-log10pvalthr
                                          , 0
                                          # ,log10pvalthr
                                          , maxMelt
                                          )
                                ,labels=c(paste("Lower mutex p:" , 
                                            if(minMelt<=-log10pvalthr) prettyNum(10^minMelt , digits=3) else "no significance")
                                          # , "MutEx Threshold"
                                          , "Random" 
                                          # , "Cooc Threshold"
                                          , paste("Lower cooc p:" , 
                                            if(maxMelt>=log10pvalthr) prettyNum(10^-maxMelt , digits=3) else "no significance")
                                          )
                                ,limits=c(minMelt,maxMelt)
                                )+
              theme_minimal()+ # minimal theme
              xlab("") + ylab("") + 
              theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                    size = 12, hjust = 1))+
              coord_fixed() +
              # geom_point(data = fakedf, aes(size="Random", shape = NA), colour = "dark gray")+
              # guides(size=guide_legend("Source", override.aes=list(shape=15, size = 10)))+
              ggtitle(title) #+
              #theme(plot.margin=unit(c(1,1,1,1), "cm"),legend.justification=c(1,0), legend.position=c(1,0))
  return(ggheatmap)
}



setGeneric('coocMutexPlot', function(object 
                                  , var=c("drug","group","gene_symbol")
                                  , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                  , grouping=c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
                                  , tumor_type=NULL
                                  , collapseMutationByGene=FALSE
                                  , collapseByGene=FALSE
                                  , tumor.weights=NULL
                                  , style=c("cooc" , "dendro")
                                  , prob = c("hyper" , "firth")
                                  , drop = FALSE
                                  , noPlot=FALSE
                                  , pvalthr=0.05
                                  , plotrandom=TRUE
                                  , ncolPlot =FALSE
                                  , ...) {
  standardGeneric('coocMutexPlot')
  })
setMethod('coocMutexPlot', 'CancerPanel', function(object 
                                  , var=c("drug","group","gene_symbol")
                                  , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                  , grouping=c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
                                  , tumor_type=NULL
                                  , collapseMutationByGene=FALSE
                                  , collapseByGene=FALSE
                                  , tumor.weights=NULL
                                  , style=c("cooc" , "dendro")
                                  , prob = c("hyper" , "firth")
                                  , drop = FALSE
                                  , noPlot=FALSE
                                  , pvalthr=0.05
                                  , plotrandom=TRUE
                                  , ncolPlot =FALSE
                                  , ...)
{
  #-------------------------------
  # CHECK PARAMETER CONSISTENCY
  #-------------------------------

  possibleAlterations <- c("copynumber" , "expression" , "mutations" , "fusions")
  possibleGrouping <- c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
  possibleProb <- c("firth" , "hyper")
  possiblestyle <- c("cooc" , "dendro")
  prob <- prob[1]
  style <- style[1]
  if(style %notin% possiblestyle){
    stop("style can only be one of the following" %++% paste(possiblestyle , collapse=", "))
  }
  if( prob %notin% possibleProb){
    stop("prob can only be one of the following" %++% paste(possibleProb , collapse=", "))
  }
  if(any(alterationType %notin% possibleAlterations)){
    stop("alterationType can only be one or more of the following" %++% paste(possibleAlterations , collapse=", "))
  }
  if( length(grouping)>1 & any(is.na(grouping)) ){
    grouping <- NA
  }
  if(!any(is.na(grouping))){
    if(any(grouping %notin% possibleGrouping))
      stop("grouping can only be one of the following:" %++% paste(possibleGrouping , collapse=", "))
  }
  if(("alteration_id" %in% grouping) & length(alterationType)<2){
    stop("If you select 'alteration_id' as grouping variable, you must select more than one alterationType")
  }
  if(!is.numeric(pvalthr)){
    stop("pvalthr must be a numeric value")
  } else {
    if(pvalthr>1 | pvalthr<=0){
      stop("pvalthr must be a number between 0 and 1")
    }
  }
  if(var[1] %eq% grouping[1]){
    stop("var and grouping cannot be to the same value")
  }
  if(!is.logical(plotrandom)){
    stop("plotrandom can be either TRUE or FALSE")
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
  rm(myenv)
  
  #----------------------------
  # RESHAPE DATA
  #----------------------------
  
  if(is.na(grouping)){
    mydata_split <- list("NA"=mydata)
  } else {
    mydata_split <- split(mydata , mydata[[grouping]])
  }
  mydata_split <- lapply(mydata_split , function(x) {
            # browser()
            tums <- unique(x$tumor_type)
            if(grouping %in% "tumor_type"){
              x$case_id <- factor(x$case_id , levels=unique(unlist(mysamples[tums])))
            } else {
              x$case_id <- factor(x$case_id , levels=unique(unlist(mysamples["all_tumors"])))
            }
            return(x)
    })
  # Create a matrix on the basis of the var value
  mat_split <- lapply(mydata_split , function(x) {
              # x$case_id <- factor(x$case_id , levels=)
              if("alteration_id" %in% colnames(x)){
                mat <- suppressMessages(
                  reshape2::acast(x
                        , as.formula( paste("case_id~" , var))
                        , fun.aggregate=length
                        , value.var="alteration_id"
                        , drop=drop))
              } else {
                mat <- suppressMessages(
                  reshape2::acast(x
                          , as.formula( paste("case_id~" , var))
                          , fun.aggregate=length
                          , drop=drop))
              }
              mat[ mat> 1] <- 1
              return(mat)
    })
  # return(mat_split)
  if(style=="dendro"){
      mat_dist <- lapply(mat_split , function(x) dist(t(x) , method="binary") )
      if(!noPlot){
      opar <- par("mfrow" , "mar")
      on.exit(par(opar))
          return({
              n_plots <- length(mat_dist)
              par(mfrow=.mfrow(n_plots, ncolPlot))
              for(i in 1:length(mat_dist)){
          title <- grouping %+% ": " %+% names(mat_dist)[i]
          if(attributes(mat_dist[[i]])$Size>2){
                    plot(hclust( mat_dist[[i]]  , ...)
                        , xlab="" 
              , main=title
              )
          } else {
            plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, "Not enough elements\nto display\n(minimum 3)"
                ,cex = 1, col = "black" , main=names(mydata_split)[i])
          }
              }
          })
      } else {
      return(mat_split)
    }
  }
  # Cooccurence 
  cooc <- lapply(1:length(mat_split) , function(i) {
            x <- mat_split[[i]]
            if(ncol(x)>1){
              out <- .coocMutEx(x , thresh=FALSE , prob=prob)
              out$grouping <- names(mydata_split)[i]
            } else {
              out <- NULL
            }
            return(out)
    })
  if(noPlot){
    return(do.call("rbind" , cooc))
  } else {
    n_plots <- length(cooc)
    # Automatic detection of plot best layout
    plot_layout <- .mfrow(n_plots, ncolPlot)
    all_melted <- lapply(1:length(cooc), function(i){
      .coocMutEx_melter(cooc[[i]] , plotrandom=plotrandom , pvalthr=pvalthr)
    })
    minMax <- lapply(all_melted , function(melted_cormat) {
      if(is.null(melted_cormat)){
        return(NULL)
      }
      if(all(is.na(melted_cormat$value ))){
        return(NULL)
      } else {
        return( c(min(melted_cormat$value , na.rm=TRUE)[1] , max(melted_cormat$value , na.rm=TRUE)[1]) )
      }
    })

    minMelt <- tryCatch( min(unlist(lapply(minMax , '[' , 1)) , na.rm=TRUE) , error=function(e) return(NULL))
    maxMelt <- tryCatch( max(unlist(lapply(minMax , '[' , 2)) , na.rm=TRUE) , error=function(e) return(NULL))
    whichIsMax <- which.max(abs(c(minMelt, maxMelt)))
    if(whichIsMax==1){
      maxMelt <- - minMelt
    } else {
      minMelt <- - maxMelt
    }
    # print(minMelt)
    # print(maxMelt)
    # print(minMax)
    # print(all_melted)
    # maxMelt <- max(melted_cormat$value , na.rm=TRUE)
    all_plots <- lapply(1:length(all_melted) , function(i){
                  if( is.na(grouping) && toupper(names(mydata_split)[i])=="NA" ){
                    title <- ""
                  } else {
                    title <- toupper(grouping) %+% ": " %+% toupper(names(mydata_split)[i])
                  }
                  if(!is.null(cooc[[i]])){
                    if(all(is.na(all_melted[[i]]$value))){
                      return(.emptyGGplot(label= "No Significant pairs", title=title))
                    }
                    out <- .plot_coocMutEx(all_melted[[i]] , title=title , minMelt=minMelt , maxMelt=maxMelt , pvalthr=pvalthr)
                    # out <- .plot_coocMutEx(cooc[[i]] , title=title , pvalthr=pvalthr , plotrandom=plotrandom)
                    # out <- out + theme(plot.margin=unit(c(1,1,1,1), "cm"),legend.justification=c(1,0), legend.position=c(1,0))
                  } else {
                    out <- .emptyGGplot(label= "MutEx and Cooc Analysis require at least 2 elements\nNo possible Mutex or Cooc Analysis"
                              , title=title
                              )
                    # out <- out + theme(plot.margin=unit(c(1,1,1,1), "cm"),legend.justification=c(1,0), legend.position=c(1,0))
                  }
                  })
    return(.multiplot(all_plots , plot_layout=plot_layout))
  }
})
