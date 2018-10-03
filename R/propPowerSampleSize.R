#---------------------------------------------------------------------------
# Minimize sample size under binomial distribution by cycling alpha and beta
# Adapted from package clinfun, https://CRAN.R-project.org/package=clinfun
#
.opt2single <- function(pControl , pCase , errorItype , beta , sample.size , mode) {
  if(mode == "sample.size"){
    n0 <- 1
    incr.n0 <- TRUE
    qerrorItype <- qbinom(1-errorItype, n0, pControl)
    err2 <- pbinom(qerrorItype, n0, pCase)
    if (err2 <= beta) incr.n0 <- FALSE
    while (incr.n0) {
      n0 <- n0 + 1
      qerrorItype <- qbinom(1-errorItype, n0, pControl)
      err2 <- pbinom(qerrorItype, n0, pCase)
      if (err2 <= beta) incr.n0 <- FALSE
    }
    return(n0)
  } else if(mode=="power"){
      qerrorItype <- qbinom(1-errorItype, sample.size, pControl)
      pow <- 1 - pbinom(qerrorItype, sample.size, pCase)
      return(pow)
  } else {
      stop("mode can only be power or sample.size")
  }
}



#----------------------------------------------------------------
# FUNCTION TWO CALCULATE POWER ANALYSIS FOR PROPORTIONS
#
.prop21Samps <- function(alpha , beta = NULL , sample.size = NULL , pCase , pControl 
    , case.fraction , side = c(2 , 1) , type = c("chisquare" , "arcsin" ,"exact") , mode = c("sample.size" , "power")){
    mode = mode[1]
    type = type[1]
    errorItype = alpha/side
    kappa = case.fraction / (1 - case.fraction)
    if(mode == "sample.size"){
        if(type=="chisquare"){
            nControl <- (pCase*(1-pCase)/kappa+pControl*(1-pControl))*((qnorm(1-errorItype)+qnorm(1-beta))/(pCase-pControl))^2
            n <- nControl + kappa*nControl
        } else if(type=="arcsin") {
            # n <- pCase*(1-pCase)*((qnorm(1-errorItype)+qnorm(1-beta))/(pCase-pControl))^2
            nonctr <- asin(sqrt(pCase)) - asin(sqrt(pControl))
            sss <- (qnorm(1 - alpha / side) + qnorm(1 - beta)) / (2 * nonctr)
            n <- sss * sss + .5
        } else if(type=="exact"){
            n <- .opt2single(pControl=pControl,pCase=pCase,errorItype=errorItype,beta=beta,sample.size=NULL,mode="sample.size")
        }
        return(n)
    } else if(mode == "power"){
        if(type == "chisquare"){
            nControl <- (1-case.fraction)*sample.size
            z <- (pCase-pControl)/sqrt(pCase*(1-pCase)/nControl/kappa+pControl*(1-pControl)/nControl)
            return( pnorm(abs(z)-qnorm(1-errorItype)) )
        } else if(type=="arcsin"){
            #z <- (pCase-pControl)/sqrt(pCase*(1-pCase)/sample.size)
            #return( pnorm(z-qnorm(1-errorItype))+pnorm(-z-qnorm(1-errorItype)) )
            t2 <- 2 * abs(asin(sqrt(pCase)) - asin(sqrt(pControl))) * sqrt(sample.size)
            pow <- pnorm((qnorm(alpha / side) + t2))
            return(pow)
        } else if(type=="exact"){
            pow <- .opt2single(pControl=pControl,pCase=pCase,errorItype=errorItype,beta=NULL,sample.size=sample.size,mode="power")
            return(pow)
        }
    } else {
        stop("mode can be either sample.size or power")
    }
}


#------------------------------------------------------------------------
# GIVEN THE OUTPUT OF .prop21Samps, CALCULATE SAMPLE.SIZE AT SCREENING
#
.propCalPowerSamp <- function(finalTabFrac , alpha , power , sample.size , pCase , pControl , case.fraction , side , type , round.result){
    # print(finalTabFrac)
    if(is.null(sample.size)){
        for(i in seq_len(length(pCase))){
            pa <- pCase[i]
            pb <- pControl[i]
            beta <- 1-power
            sample.size <- lapply(beta , function(x) {
                    # ((qnorm(1-alpha/2)+qnorm(1-x))/(log(hr)-log(hr0)))^2/(case.fraction*(1-case.fraction)*p.event)
                    .prop21Samps(alpha = alpha , beta = x , pCase = pa, pControl = pb , case.fraction = case.fraction, side = side, type = type , mode = "sample.size")
                    }) %>% unlist
            sample.size <- if( round.result ) ceiling(sample.size) else sample.size
            sample.size_recruit <- lapply(finalTabFrac , function(x) {
                multiplier <- ifelse(x!=0 , 1/x , NA)
                superSample <- if( round.result ) ceiling(sample.size*multiplier) else sample.size*multiplier
                return(superSample)
                } ) %>% do.call("cbind" , .)
            if(i==1){
                toBePlot <- data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size) , length(names(finalTabFrac))))
                                    , "ScreeningSampleSize" = as.vector(sample.size_recruit)
                                    , "EligibleSampleSize" = rep(sample.size , length(finalTabFrac))
                                    , Beta = rep(beta , length(finalTabFrac))
                                    , Power = rep(power , length(finalTabFrac))
                                    , "Proportion.In.Case.Control" = as.character(paste( pa , pb , sep=" - ") )
                                    , stringsAsFactors=FALSE)
            } else {
                toBePlot <- rbind(toBePlot
                                , data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size) , length(names(finalTabFrac))))
                                    , "ScreeningSampleSize" = as.vector(sample.size_recruit)
                                    , "EligibleSampleSize" = rep(sample.size , length(finalTabFrac))
                                    , Beta = rep(beta , length(finalTabFrac))
                                    , Power = rep(power , length(finalTabFrac))
                                    , "Proportion.In.Case.Control" = as.character(paste( pa , pb , sep=" - ") )
                                    , stringsAsFactors=FALSE))
            }
        }
    }
    if(is.null(power)){
        sample.size_recruit <- sample.size
        sample.size <- lapply(finalTabFrac , function(x) {
                multiplier <- ifelse(x!=0 , x , NA)
                superSample <- sample.size_recruit*multiplier
                superSample <- if( round.result ) ceiling(superSample) else superSample
                return(superSample)
                } ) %>% do.call("cbind" , .)
        for(i in seq_len(length(pCase))){
            pa <- pCase[i]
            pb <- pControl[i]
            power <- lapply(as.vector(sample.size) , function(x) {
                # z <- (log(hr)-log(hr0))*sqrt(x*case.fraction*(1-case.fraction)*p.event)
                # pnorm(z - qnorm(1-alpha/2)) + pnorm( -z - qnorm(1-alpha/2))
                .prop21Samps(alpha = alpha , sample.size = x , pCase = pa, pControl = pb , case.fraction = case.fraction, side = side, type = type , mode = "power")
            }) %>% unlist
            # print(power)
            beta <- 1 - power
            if(i==1){
                toBePlot <- data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size_recruit) , ncol(sample.size)))
                                    , "ScreeningSampleSize" = rep(sample.size_recruit , ncol(sample.size))
                                    , "EligibleSampleSize" = as.vector(sample.size)
                                    , Beta = beta
                                    , Power = power
                                    , "Proportion.In.Case.Control" = as.character(paste( pa , pb , sep=" - ") )
                                    , stringsAsFactors=FALSE)
            } else {
                toBePlot <- rbind(toBePlot
                                , data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size_recruit) , ncol(sample.size)))
                                            , "ScreeningSampleSize" = rep(sample.size_recruit , ncol(sample.size))
                                            , "EligibleSampleSize" = as.vector(sample.size)
                                        , Beta = beta
                                        , Power = power
                                        , "Proportion.In.Case.Control" = as.character(paste( pa , pb , sep=" - ") )
                                        , stringsAsFactors=FALSE))
            }
        }
    }
    return(toBePlot)
}


setGeneric('propPowerSampleSize', function(object
                                    , var=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                    , tumor_type=NULL
                                    , stratum=NULL
                                    , tumor.weights=NULL
                                    , tumor.freqs=NULL
                                    , pCase=NULL
                                    , pControl=NULL
                                    , side = c(2,1)
                                    , type = c("chisquare" , "arcsin" , "exact")
                                    , alpha=0.05
                                    , power=NULL
                                    , sample.size=NULL
                                    , case.fraction=0.5
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=FALSE
                                    , round.result=TRUE
                                    , priority.trial=NULL
                                    , priority.trial.order=c("optimal" , "as.is")
                                    , priority.trial.verbose=TRUE
                                    , noPlot=FALSE) {
    standardGeneric('propPowerSampleSize')
    })
setMethod('propPowerSampleSize', 'CancerPanel', function(object
                                    , var=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                    , tumor_type=NULL
                                    , stratum=NULL
                                    , tumor.weights=NULL
                                    , tumor.freqs=NULL
                                    , pCase=NULL
                                    , pControl=NULL
                                    , side = c(2,1)
                                    , type = c("chisquare" , "arcsin" , "exact")
                                    , alpha=0.05
                                    , power=NULL
                                    , sample.size=NULL
                                    , case.fraction=0.5
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=FALSE
                                    , round.result=TRUE
                                    , priority.trial=NULL
                                    , priority.trial.order=c("optimal" , "as.is")
                                    , priority.trial.verbose=TRUE
                                    , noPlot=FALSE)
{
    # Check var
    possiblevar <- c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
    if(any(is.na(var))){
        var <- NA
    }
    if(length(var)!=1){
        var <- var[1]
    }
    if(any(var %notin% possiblevar)){
        stop("var can only be one of the following" %++% paste(possiblevar , collapse=", "))
    }
    # Check tumor_type
    if(!is.null(tumor_type)){
        if(!all(tumor_type %in% object@arguments$tumor_type)){
            stop("You selected a tumor_type that has no data in this CancerPanel object")
        }
    }
    # Side parameter decide to perform a one-side or two-side calculation
    side <- side[1]
    requiredParam <- list(pCase=pCase , pControl=pControl , side=side , case.fraction=case.fraction , alpha=alpha)
    # Sanity check for power calculation numerical parameters
    for(i in names(requiredParam)){
        if(is.null(requiredParam[[i]]))
            stop(i %++% "cannot be NULL")
        if(!is.numeric(requiredParam[[i]]))
            stop(i %++% "must be numeric")
        if(any(requiredParam[[i]]<=0))
            stop(i %++% "must be a positive number")
        if(i %in% c("case.fraction" , "alpha")){
            if(any(requiredParam[[i]]>=1))
                stop(i %++% "must be a number strictly between 0 and 1")
        }
        if(i %in% c("pCase" , "pControl")){
            if(any(requiredParam[[i]]>1))
                stop(i %++% "must be a number strictly between 0 and 1")
        }
        if(i=="side"){
            if(requiredParam[[i]] %notin% c(2,1))
                stop(i %++% "can only be 1 or 2")
        }
    }
    # Sanity Check for type
    type <- type[1]
    if(type %notin% c("chisquare" , "arcsin" , "exact")){
        stop("type can only be arcsin, chisquare or exact")
    }
    # Sanity check for pCase and pControl
    if(length(pCase)!=length(pControl)){
        stop("pCase and pControl must be of the same length")
    }
    # Sanity check for power and sample.size
    optionalParam <- list(power=power , sample.size=sample.size )
    if(all(vapply(optionalParam , is.null , logical(1)))) {
        stop("power and sample.size cannot be both NULL")
    }
    if(all(!vapply(optionalParam , is.null , logical(1)))) {
        stop("power and sample.size cannot be both set")
    }
    if(!is.null(power)){
        if(!is.numeric(power)){
            stop("power must be a numerical vector")
        }
        if(any(power>1) | any(power <= 0)){
            stop("power must be a numerical vector with elements between 0 and 1")
        }
    }
    if(!is.null(sample.size)){
        if(!is.numeric(sample.size)){
            stop("sample.size must be a numerical vector")
        }
        if( any(sample.size<2) ){
            stop("sample.size must be a numerical vector with elements > 2")
        }
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
        if(!is.null(priority.trial)){
            warning("tumor.freqs are ignored if priority.trial is set")
            tumor.freqs <- NULL
        }
    }
    # Check for data parameters
    possibleAlterations <- c("copynumber" , "expression" , "mutations" , "fusions")
    possibleVar <- c("drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
    if(length(var)>1){
        var <- var[1]
    } else if(length(var)==0 | is.null(var)) {
        var <- NA
    }

    if(("alteration_id" %in% var) & collapseByGene){
        warning("If 'alteration_id' is var, you cannot collapse by gene. The option was set to FALSE")
        collapseByGene <- FALSE
    }

    # Check priority.trial parameter
    if(!is.null(priority.trial)){
        if(!is.character(priority.trial)){
            stop("priority.trial must be a character vector")
        }
        if(length(priority.trial)<2){
            stop("priority.trial must be a character vector containing at least 2 elements")
        }
        if(var %notin% c("drug" , "group")){
            stop("priority.trial parameter can be set only if var is drug or group")
        }
        if(var=="drug"){
            if(!all(priority.trial %in% object@arguments$drugs)){
                stop("The following drug names are not included in the object:" %++% 
                    paste(priority.trial[priority.trial %notin% object@arguments$drugs] , collapse=", "))
            }
        }
        if(var=="group"){
            if(!all(priority.trial %in% object@arguments$panel$group)){
                stop("The following group levels are not included in the object:" %++% 
                    paste(priority.trial[priority.trial %notin% object@arguments$panel$group] , collapse=", "))
            }
        }
        priority.trial.order <- priority.trial.order[1]
        if(priority.trial.order %notin% c("optimal" , "as.is")){
            stop("priority.trial.order can be either optimal or as.is")
        }
        if(is.null(power)){
            stop("priority.trial can be used only to calculate sample size, power cannot be NULL")
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
    rm(myenv);gc()
    # Detect var levels
    varLevels <- switch(var 
                        ,gene_symbol=unique(c(object@arguments$panel$gene_symbol , mydata$gene_symbol))
                        ,drug=unique(c(object@arguments$panel$drug , mydata$drug))
                        ,group=unique(c(object@arguments$panel$group , mydata$group))
                        ,alteration_id=c("cna" , "expr" ,"fus" , "mut")
                        ,tumor_type=unique(mydata$tumor_type))
    # Filter out strata not reqeusted
    if(!is.null(stratum)){
        if(!all(stratum %in% varLevels)){
            stop("The strata requested are not present in the data under var")
        }
        mydata <- mydata[ mydata[ , var] %in% stratum , ]
        varLevels <- stratum
    }
    # browser()
    if(!is.na(var)){
        myDataForTotal <- unique(mydata[ , c(var , "case_id")])
        finalTab <- table( factor(myDataForTotal[ , var] , levels=sort(varLevels)) ) %>% as.list %>% unlist
        finalTab['Full Design'] <- length(unique(myDataForTotal$case_id))
        if(var!="tumor_type"){
            finalTabFrac <- finalTab/length(mysamples[["all_tumors"]])
        } else {
            finalTabFrac <- c()
            for( i in names(finalTab)){
                if(i=="Full Design")
                    finalTabFrac['Full Design'] <- finalTab['Full Design']/length(mysamples[["all_tumors"]])
                else
                    finalTabFrac[i] <- finalTab[i]/length(mysamples[[i]])
            }
        }
    } else {
        finalTab <- c("Full Design"=length(unique(mydata$case_id)))
        finalTabFrac <- finalTab/length(mysamples[["all_tumors"]])
    }
    # If tumor.freqs is set, we basically run this method in recursive mode for every tumor type
    # and then aggregate everything using the weights of tumor.freqs
    if(!is.null(tumor.freqs)){
        if(!is.na(var) && var=="tumor_type"){
            finalTabFrac["Full Design"] <- sum( tumor.freqs*finalTabFrac[names(tumor.freqs)] , na.rm=TRUE)
            toBePlot <- .propCalPowerSamp(finalTabFrac , alpha , power , sample.size , pCase , pControl , case.fraction , side , type , round.result)
        } else {
            survRecurse <- lapply(names(tumor.freqs) , function(tum){
                out <- tryCatch( suppressWarnings(propPowerSampleSize(object
                                        , var=var
                                        , alterationType=alterationType
                                        , tumor_type=tum
                                        , pCase=pCase
                                        , pControl=pControl
                                        , side=side
                                        , alpha=alpha
                                        , power=power
                                        , sample.size=sample.size
                                        , case.fraction=case.fraction
                                        , type=type
                                        , collapseMutationByGene=collapseMutationByGene
                                        , collapseByGene=collapseByGene
                                        , noPlot=TRUE
                                        , round.result=FALSE))
                                , error=function(e) return(NULL))
                if(is.null(out)){
                    return(NULL)
                }
                Rate <- paste0("rate" , tum)
                out[ , Rate] <- out$EligibleSampleSize/out$ScreeningSampleSize
                finalTabFracTum <- aggregate(as.formula( paste(Rate , "~ Var") ) , out , FUN=mean)
                finalTabFracTum[ , Rate] <- finalTabFracTum[ , Rate]*tumor.freqs[tum]
                return(finalTabFracTum)
            })
            # print(str(survRecurse))
            survRecurseWeighted <- .mergeThemAll(survRecurse[!vapply(survRecurse , is.null , logical(1))] 
                                                 , by="Var" , all=TRUE)
            rownames(survRecurseWeighted) <- survRecurseWeighted$Var
            survRecurseWeighted$Var <- NULL
            survRecurseWeighted_vec <- rowSums(survRecurseWeighted , na.rm=TRUE)
            # toBePlot <- .propCalPowerSamp(finalTabFrac=survRecurseWeighted_vec , sample.size , pCase , pControl , power , alpha , p.event , case.fraction , round.result=round.result)
            toBePlot <- .propCalPowerSamp(finalTabFrac=survRecurseWeighted_vec , alpha , power , sample.size , pCase , pControl , case.fraction , side , type , round.result)
        }
    } else {
        # toBePlot <- .propCalPowerSamp(finalTabFrac , sample.size , pCase , pControl , power , alpha , p.event , case.fraction , round.result=round.result)
        toBePlot <- .propCalPowerSamp(finalTabFrac , alpha , power , sample.size , pCase , pControl , case.fraction , side , type , round.result)
    }
    #---------------------
    # PRIORITY TRIAL
    #---------------------
    # subset finalTabFrac according to priority.trial variables
    if(!is.null(priority.trial)) {
        powerPropMat <- data.frame(pCase = rep(pCase , length(power))
                                , pControl = rep(pControl , length(power))
                                , Power = lapply( power , function(x) rep(x , length(pCase))) %>% unlist
                                , stringsAsFactors = FALSE)
        priorTrial <- lapply(seq_len(nrow(powerPropMat)) , function(row) {
            pCase <- powerPropMat[row , "pCase"]
            pControl <- powerPropMat[row , "pControl"]
            power <- powerPropMat[row , "Power"]
            if(priority.trial.order=="optimal"){
                priority.trial_frac <- finalTabFrac[priority.trial] %>% sort
            } else {
                priority.trial_frac <- finalTabFrac[priority.trial]
            }
            k <- length(priority.trial_frac)
            toBePlot_prior <- .propCalPowerSamp(priority.trial_frac , alpha , power , sample.size , pCase , pControl , case.fraction , side , type , round.result)
            # print(toBePlot_prior)
            TSS <- toBePlot_prior[ toBePlot_prior$Var==names(priority.trial_frac)[1] , "EligibleSampleSize"]
            if(noPlot && priority.trial.verbose){
                message("\nMinimum Eligible Sample Size to reach" %++% 
                            paste0(round(power*100),"%") %++% 
                            "of power with a pCase - pControl equal to" %++% paste(pCase,pControl,sep=" - ") %++% 
                            "for each" %++% var %++% 
                            "is equal to:" %++% TSS)
            }
            # Initialize the three matrices
            matProb <- matSamps <- matSampsAssigned <- matrix(0 , ncol=k , nrow=k)
            rownames(matProb) <- rownames(matSamps) <- rownames(matSampsAssigned) <- names(priority.trial_frac)
            colnames(matProb) <- colnames(matSamps) <- colnames(matSampsAssigned) <- names(priority.trial_frac)
            # This matrix represents the order of j to follow after the recruitment
            matForOrder <- diag(1 , nrow=k , ncol=k)
            for(i in seq_len(nrow(matForOrder))){
                matForOrder[ i , which(matForOrder[i,]!=1)] <- 2:k
            }
            # i start a new run of Screening, j take care of the discarded, to check if they can fit into another drug
            colsAndRows <- seq_len(k)
            for(i in seq_len(k)){
                # If all drugs have reached TSS, stop Screening
                if(all(colSums(matSampsAssigned)>=TSS)){
                    break
                }
                # drugsDone <- c() 
                vec <- c(i , colsAndRows[ colsAndRows!=i])
                for(j in vec) {
                    if(i == j){
                        matProb[ i , j ] <- unname(priority.trial_frac[i])
                        drugsDone <- names(priority.trial_frac)[i]
                        if(i == 1){
                            matSampsAssigned[ i , j ] <- TSS
                            matSamps[i , j] <- toBePlot_prior[ toBePlot_prior$Var==names(priority.trial_frac)[i] , "ScreeningSampleSize"] - TSS
                        } else {
                            if(sum(matSampsAssigned[ , j])>=TSS){
                                matSamps[i,j] <- 0
                            } else {
                                # Here, we deal with new screenings. A new screening X is calculated as X = (TSS - already assigned)/p
                                # Then we directly subtract what we expect to obtain from this recruitment (TSS - already assigned)
                                matSamps[i , j] <- (TSS - sum(matSampsAssigned[ , j]))/matProb[ i , j ] - (TSS - sum(matSampsAssigned[ , j]))
                                matSampsAssigned[ i , j ] <- TSS - sum(matSampsAssigned[ , j])
                            }
                        }
                    } else {
                        # If there are no leftovers, stop cycling on columns (j)
                        if(matSamps[i , j_old]>=0){
                            # current drug
                            myDrug <- names(priority.trial_frac)[ j ]
                            # Remove all the samples with a drug already taken into consideration
                            sampToRemove <- myDataForTotal[ myDataForTotal[ , var] %in% drugsDone , "case_id"] %>% unique
                            # By difference, create the sample set of the remaining samples
                            sampToKeep <- setdiff(mysamples[["all_tumors"]] , sampToRemove)
                            # Number of alterations found for my drug under 
                            altNum <- myDataForTotal[ myDataForTotal[ , var] %in% myDrug & myDataForTotal$case_id %in% sampToKeep, "case_id"] %>% unique %>% length
                            matProb[ i , j ] <- altNum/length(sampToKeep)
                            # obtainable samples are calculated as the samples remained from the previous cycle * probability of obtaining new samples
                            sampObtainable <- matProb[ i , j ] * matSamps[i , j_old]
                            # update the set of drugs already screened
                            drugsDone <- c(drugsDone , myDrug)
                            # If it is the last Allocation, take all the patient you can get
                            if(which(vec==j) == k){
                                matSampsAssigned[ i , j ] <- sampObtainable
                            } else {
                                # If all the next Allocations already reached TSS, take all the patients you can get
                                if( all( colSums(matSampsAssigned[  ,  setdiff(names(matSampsAssigned) , drugsDone) , drop=FALSE]) >= TSS )){
                                    matSampsAssigned[ i , j ] <- sampObtainable
                                } else {
                                    if(sum(matSampsAssigned[ , j])>=TSS){
                                        matSampsAssigned[ i , j ] <- 0
                                    } else {
                                        # If samples obtainable are not sufficient to reach TSS, assign 'em all
                                        # If they are sufficient, assigned just what is enough to reach TSS
                                        if(TSS > sampObtainable + sum(matSampsAssigned[ , j])){
                                            matSampsAssigned[ i , j ] <- sampObtainable
                                        } else {
                                            matSampsAssigned[ i , j ] <- TSS - sum(matSampsAssigned[ , j])
                                        }
                                    }
                                }
                            }
                            # Remaining samples are the ones from the previous cycle - what was assigned
                            matSamps[ i , j ] <- matSamps[ i , j_old] - matSampsAssigned[i , j]
                            if(matSamps[i , j]<=0){
                                break
                            }
                        } else {
                            break
                        }
                    }
                    # keep track of the index of the current cycle for the next one
                    j_old <- j
                }
            }
            # matSamps is already deprived of the samples that can go to testing, so we readd them
            matSampsComplete <- matSamps + matSampsAssigned
            matSampsAssignedWithTotal <- rbind( matSampsAssigned , colSums(matSampsAssigned))
            rownames(matSampsAssignedWithTotal)[nrow(matSampsAssignedWithTotal)] <- "Total"

            torecruit <- ceiling(c(diag(matSampsComplete) , Total = sum(diag(matSampsComplete))))
            toassign <- ceiling(c(colSums(matSampsAssigned) , Total = sum(colSums(matSampsAssigned))))
            discard <- matSampsComplete[ , ncol(matSampsComplete)] - matSampsAssigned[ , ncol(matSampsAssigned)]
            todiscard <- ceiling(c(discard , Total = sum(discard)))
            return(list( Summary = rbind( Screened = torecruit
                                        , Eligible = toassign
                                        , Not.Eligible = todiscard)
                        , Screening.scheme = ceiling(matSampsComplete)
                        , Allocation.scheme = ceiling(matSampsAssignedWithTotal)
                        , Probability.scheme = matProb
                        , Base.probabilities = priority.trial_frac))
        })
        names(priorTrial) <- apply( powerPropMat , 1 , function(x) paste( paste(c("pCase" , "pControl" , "Power") , x , sep=":") , collapse=" | "))
        if(noPlot){
            if(priority.trial.verbose){
                toBePlot <- priorTrial
            } else {
                toBePlot <- powerPropMat
                toBePlot$Var <- "Priority Design"
                toBePlot$ScreeningSampleSize <- vapply(priorTrial , function(x) x[[1]][ "Screened", "Total"] , numeric(1))
                toBePlot$EligibleSampleSize <- vapply(priorTrial , function(x) x[[1]][ "Eligible", "Total"] , numeric(1))
                toBePlot$Proportion.In.Case.Control <- paste( toBePlot$pCase , toBePlot$pControl , sep = " - ")
                toBePlot$Beta <- 1 - toBePlot$Power
                toBePlot <- toBePlot[ , c("Var","ScreeningSampleSize","EligibleSampleSize","Beta","Power","Proportion.In.Case.Control")]
            }
            return(toBePlot)
        } else {
            toBePlot <- powerPropMat
            toBePlot$Var <- "Priority Design"
            toBePlot$ScreeningSampleSize <- vapply(priorTrial , function(x) x[[1]][ "Screened", "Total"] , numeric(1))
            # toBePlot$Proportion.In.Case.Control <- as.character(toBePlot$Proportion.In.Case.Control)
            toBePlot$Proportion.In.Case.Control <- paste( toBePlot$pCase , toBePlot$pControl , sep = " - ")
            # toBePlot$pControl <- NULL
        }
    }
    # # If stratum is set, only part of the levels of var are shown
    # if(!is.null(stratum)){
    #     if(all(stratum %in% toBePlot$Var)){
    #         toBePlot <- toBePlot[ toBePlot$Var %in% stratum , ] 
    #     } else {
    #         warning("The stratum requested are not present in var, returning all the dataset")
    #     }
    # }
    # Plot or return the table
    if(noPlot){
        return(toBePlot)
    } else {
        Power=ScreeningSampleSize=Proportion.In.Case.Control=NULL
        plotTitle <- "Screening sample size by power stratified by pCase levels"
        if("Full Design" %in% toBePlot$Var){
            # Show Full Design always as the last plot
            toBePlot$Var <- factor(toBePlot$Var , levels=c( setdiff( unique(toBePlot$Var) , "Full Design") , "Full Design"))
        }
        p <- ggplot(toBePlot , aes(y=Power, x=ScreeningSampleSize, colour=Proportion.In.Case.Control , group=Proportion.In.Case.Control)) +
            facet_grid(Var ~ . , scales="free_x") +
            geom_line() +
            geom_point(size=2 ) + 
            xlab("Screening Sample Size") +
            ylab("Power") +
            theme_bw() +
            ggtitle(plotTitle)
        return(p)
    }
})

