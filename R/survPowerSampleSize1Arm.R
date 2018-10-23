#-------------------------------------------------------------
# FUNCTION TWO CALCULATE 1-SAMPLE POWER ANALYSIS FOR SURVIVAL
.surv1Samps <- function(alpha = 0.05 , beta = NULL , sample.size = NULL , med0 , med1 , mode = c("sample.size" , "power") , fu , acc = NULL , side ){
 mode = mode[1]
 ber0 = -log(.5) / med0
 ber1 = -log(.5) / med1
 # Default accrual is infinite, here expressed as NULL
 if(!is.null(acc)){
     acc0 <- (1 - exp(-ber0 * acc)) / (ber0 * acc)
     acc1 <- (1 - exp(-ber1 * acc)) / (ber1 * acc)
 } else {
     acc0 <- 1
     acc1 <- 1
 }
 pi0 <- 1 - exp(-ber0 * fu) * acc0
 pi1 <- 1 - exp(-ber1 * fu) * acc1
 psi0 <- ber0^.3333333
 psi1 <- ber1^.3333333
 var0 <- psi0 * psi0 / (9 * pi0)
 var1 <- psi1 * psi1 / (9 * pi1)
 if(mode == "sample.size"){
    nonctr <- abs(psi1 - psi0)
    sss <- (qnorm(1 - alpha / side) + qnorm(1-beta) * sqrt(var1 / var0)) / nonctr
    sss <- sss * sqrt(var0)
    n <- round(sss * sss + .5)
    return(n)
 } else if (mode == "power") {
    t2 <- (abs(psi1 - psi0) * sqrt(sample.size)) / sqrt(var1)
    power <- pnorm((sqrt(var0 / var1) * qnorm(alpha / side) + t2))
    return(power)
 } else {
    stop("mode not specified")
 }
}

#------------------------------------------------------------------------
# GIVEN THE OUTPUT OF .surv1Samps, CALCULATE SAMPLE.SIZE AT RECRUITMENT
.calPowerSamp1 <- function(finalTabFrac , sample.size , MED1 , MED0 , power , alpha , round.result , fu , acc = NULL , side){
    # print(finalTabFrac)
    if(is.null(sample.size)){
        for(i in seq_len(length(MED1))){
            med1 <- MED1[i]
            med0 <- MED0[i]
            beta <- 1-power
            sample.size <- lapply(beta , function(x) {
                    .surv1Samps(alpha = alpha , beta = x , med1 = med1 , med0 = med0 
                         , mode = "sample.size" , fu = fu , acc = acc , side=side)
                    }) %>% unlist
            sample.size <- if( round.result ) ceiling(sample.size) else sample.size
            sample.size_recruit <- lapply(finalTabFrac , function(x) {
                multiplier <- ifelse(x!=0 , 1/x , NA)
                superSample <- if( round.result ) ceiling(sample.size*multiplier) else sample.size*multiplier
                return(superSample)
                } ) %>% do.call("rbind" , .)
            if(i==1){
                toBePlot <- data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size) , length(names(finalTabFrac))))
                                    , "ScreeningSampleSize" = as.vector(sample.size_recruit)
                                    , "EligibleSampleSize" = rep(sample.size , length(finalTabFrac))
                                    , Beta = rep(beta , length(finalTabFrac))
                                    , Power = rep(power , length(finalTabFrac))
                                    , MedianSurvivalCase=as.character(med1)
                                    , stringsAsFactors=FALSE)
            } else {
                toBePlot <- rbind(toBePlot
                                , data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size) , length(names(finalTabFrac))))
                                    , "ScreeningSampleSize" = as.vector(sample.size_recruit)
                                    , "EligibleSampleSize" = rep(sample.size , length(finalTabFrac))
                                    , Beta = rep(beta , length(finalTabFrac))
                                    , Power = rep(power , length(finalTabFrac))
                                    , MedianSurvivalCase=as.character(med1)
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
        for(i in seq_len(length(MED1))){
            med1 <- MED1[i]
            med0 <- MED0[i]
            power <- lapply(as.vector(sample.size) , function(x) {
                .surv1Samps(alpha = alpha , sample.size = x , med1 = med1 , med0 = med0 
                     , mode = "power" , fu=fu , acc=acc , side=side)
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
                                    , MedianSurvivalCase=as.character(med1)
                                    , stringsAsFactors=FALSE)
            } else {
                toBePlot <- rbind(toBePlot
                                , data.frame("Var" = rep(names(finalTabFrac) 
                                            , rep(length(sample.size_recruit) , ncol(sample.size)))
                                            , "ScreeningSampleSize" = rep(sample.size_recruit , ncol(sample.size))
                                            , "EligibleSampleSize" = as.vector(sample.size)
                                        , Beta = beta
                                        , Power = power
                                        , MedianSurvivalCase=as.character(med1)
                                        , stringsAsFactors=FALSE))
            }
        }
    }
    return(toBePlot)
}


setGeneric('survPowerSampleSize1Arm', function(object
                                    , var=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                    , tumor_type=NULL
                                    , stratum=NULL
                                    , tumor.weights=NULL
                                    , tumor.freqs=NULL
                                    , MED1=NULL
                                    , MED0=NULL
                                    , fu=2
                                    , acc=NULL
                                    , alpha=0.05
                                    , power=NULL
                                    , sample.size=NULL
                                    , side = c(2,1)
                                    , collapseMutationByGene=TRUE
                                    , collapseByGene=FALSE
                                    , round.result=TRUE
                                    , priority.trial=NULL
                                    , priority.trial.order=c("optimal" , "as.is")
                                    , priority.trial.verbose=TRUE
                                    , noPlot=FALSE) {
    standardGeneric('survPowerSampleSize1Arm')
    })
setMethod('survPowerSampleSize1Arm', 'CancerPanel', function(object
                                    , var=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
                                    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
                                    , tumor_type=NULL
                                    , stratum=NULL
                                    , tumor.weights=NULL
                                    , tumor.freqs=NULL
                                    , MED1=NULL
                                    , MED0=NULL
                                    , fu=2
                                    , acc=NULL
                                    , alpha=0.05
                                    , power=NULL
                                    , sample.size=NULL
                                    , side=c(2,1)
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
        stop(paste("var can only be one of the following" , paste(possiblevar , collapse=", ")))
    }
    # Check tumor_type
    if(!is.null(tumor_type)){
        if(!all(tumor_type %in% object@arguments$tumor_type)){
            stop("You selected a tumor_type that has no data in this CancerPanel object")
        }
    }
    # Check accrual time parameter acc
    if(!is.null(acc)){
        acc <- acc[1]
        if(!is.numeric(acc))
            stop("acc must be numeric")
        if(acc<=0)
            stop("acc must be a positive number")
    }
    side <- side[1]
    requiredParam <- list(MED1=MED1 , MED0=MED0 , alpha=alpha , fu=fu , side=side)
    # Sanity check for power calculation numerical parameters
    for(i in names(requiredParam)){
        if(is.null(requiredParam[[i]])){
            stop(paste(i , "cannot be NULL"))
        }
        if(!is.numeric(requiredParam[[i]])){
            stop(paste(i , "must be numeric"))
        }
        if(any(requiredParam[[i]]<=0)){
            stop(paste(i , "must be a positive number"))
        }
        if(i %in% c("alpha")){
            if(any(requiredParam[[i]]>=1)){
                stop(paste(i , "must be a number strictly between 0 and 1"))
            }
        }
      if(i == "side"){
        if(requiredParam[[i]] %notin% c(2,1)){
          stop(paste(i , "can be either 2 or 1"))
        }
      }
    }
    # Sanity check for MED1 and MED0
    if(length(MED1)!=length(MED0)){
        # if(is.null(MED0) | is.na(MED0))
        MED0 <- rep(MED0 , length(MED1))
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
                stop(paste("The following drug names are not included in the object:" ,
                    paste(priority.trial[priority.trial %notin% object@arguments$drugs] , collapse=", ")))
            }
        }
        if(var=="group"){
            if(!all(priority.trial %in% object@arguments$panel$group)){
                stop(paste("The following group levels are not included in the object:" ,
                    paste(priority.trial[priority.trial %notin% object@arguments$panel$group] , collapse=", ")))
            }
        }
      if(any(duplicated(priority.trial))){
        stop("Duplicated values in priority.trial parameter")
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

    de <- dataExtractor(object=object , alterationType=alterationType 
                        , tumor_type=tumor_type 
                        , collapseMutationByGene=collapseMutationByGene 
                        , collapseByGene=collapseByGene 
                        , tumor.weights=tumor.weights)
    mydata <- de$data
    mysamples <- de$Samples
    tum_type_diff <- de$tumor_not_present
    rm(de)
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
            toBePlot <- .calPowerSamp1(finalTabFrac , sample.size=sample.size , MED1=MED1 , MED0=MED0
                , power=power , alpha=alpha  , round.result=round.result , fu=fu , acc=acc , side=side)
        } else {
            survRecurse <- lapply(names(tumor.freqs) , function(tum){
                out <- tryCatch( suppressWarnings(survPowerSampleSize1Arm(object
                                        , var=var
                                        , alterationType=alterationType
                                        , tumor_type=tum
                                        , MED1=MED1
                                        , MED0=MED0
                                        , fu=fu
                                        , acc=acc
                                        , side=side
                                        , alpha=alpha
                                        , power=power
                                        , sample.size=sample.size
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
            toBePlot <- .calPowerSamp1(finalTabFrac=survRecurseWeighted_vec , sample.size=sample.size , MED1=MED1 , MED0=MED0
                , power=power , alpha=alpha  , round.result=round.result , fu=fu , acc=acc , side=side)
        }
    } else {
        toBePlot <- .calPowerSamp1(finalTabFrac , sample.size=sample.size , MED1=MED1 , MED0=MED0
                , power=power , alpha=alpha  , round.result=round.result , fu=fu , acc=acc , side=side)
    }
    #---------------------
    # PRIORITY TRIAL
    #---------------------
    # subset finalTabFrac according to priority.trial variables
    if(!is.null(priority.trial)) {
        powerMEDmat <- data.frame(MedianSurvivalCase = rep(MED1 , length(power))
                                , MED0 = rep(MED0 , length(power))
                                , Power = lapply( power , function(x) rep(x , length(MED1))) %>% unlist
                                , stringsAsFactors = FALSE)
        priorTrial <- lapply(seq_len(nrow(powerMEDmat)) , function(row) {
            MED1 <- powerMEDmat[row , "MedianSurvivalCase"]
            MED0 <- powerMEDmat[row , "MED0"]
            power <- powerMEDmat[row , "Power"]
            if(priority.trial.order=="optimal"){
                priority.trial_frac <- finalTabFrac[priority.trial] %>% sort
            } else {
                priority.trial_frac <- finalTabFrac[priority.trial]
            }
            k <- length(priority.trial_frac)
            toBePlot_prior <- .calPowerSamp1(priority.trial_frac , sample.size=sample.size , MED1=MED1 , MED0=MED0
                , power=power , alpha=alpha  , round.result=round.result , fu=fu , acc=acc , side=side)
            # print(toBePlot_prior)
            TSS <- toBePlot_prior[ toBePlot_prior$Var==names(priority.trial_frac)[1] , "EligibleSampleSize"]
            if(noPlot && priority.trial.verbose){
                message(paste("\nMinimum Eligible Sample Size to reach" , 
                            paste0(round(power*100),"%") , 
                            "of power with an MedianSurvivalCase equal to" , MED1 , 
                            "and a MedianSurvivalCtrl equal to" , MED0 ,
                            "for each" , var ,
                            "is equal to:" , TSS))
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
            # matSamps is already deprived of the samples that can go to testing, so we read them
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
        names(priorTrial) <- apply( powerMEDmat , 1 , function(x) paste( paste(c("MED1" , "MED0" , "Power") , x , sep=":") , collapse=" | "))
        if(noPlot){
            if(priority.trial.verbose){
                toBePlot <- priorTrial
            } else {
                toBePlot <- powerMEDmat
                toBePlot$Var <- "Priority Design"
                toBePlot$ScreeningSampleSize <- vapply(priorTrial , function(x) x[[1]][ "Screened", "Total"] , numeric(1))
                toBePlot$EligibleSampleSize <- vapply(priorTrial , function(x) x[[1]][ "Eligible", "Total"] , numeric(1))
                toBePlot$Beta <- 1 - toBePlot$Power
                toBePlot <- toBePlot[ , c("Var","ScreeningSampleSize","EligibleSampleSize","Beta","Power","MedianSurvivalCase")]
            }
            return(toBePlot)
        } else {
            toBePlot <- powerMEDmat
            toBePlot$Var <- "Priority Design"
            toBePlot$ScreeningSampleSize <- vapply(priorTrial , function(x) x[[1]][ "Screened", "Total"] , numeric(1))
            toBePlot$MedianSurvivalCase <- as.character(toBePlot$MedianSurvivalCase)
        }
    }
    # # If stratum is set, only part of the levels of var are shown
    # if(!is.null(stratum)){
    #   if(all(stratum %in% toBePlot$Var)){
    #       toBePlot <- toBePlot[ toBePlot$Var %in% stratum , ] 
    #   } else {
    #       warning("The stratum requested are not present in var, returning all the dataset")
    #   }
    # }
    # Plot or return the table
    if(noPlot){
        return(toBePlot)
    } else {
        Power=ScreeningSampleSize=MedianSurvivalCase=NULL
        plotTitle <- "Screening sample size by power and MedianSurvivalCase levels"
        if("Full Design" %in% toBePlot$Var){
            # Show Full Design always as the last plot
            toBePlot$Var <- factor(toBePlot$Var , levels=c( setdiff( unique(toBePlot$Var) , "Full Design") , "Full Design"))
        }
        # Median survival order must be the numeric one (not the character one)
        # (10 , 9 , 8 , 0) not (9 , 8 , 10 , 0)
        toBePlot$MedianSurvivalCase <- factor(toBePlot$MedianSurvivalCase 
                            , levels = unique(toBePlot$MedianSurvivalCase) %>% 
                                        as.numeric %>%
                                        sort %>%
                                        as.character)
        p <- ggplot(toBePlot , aes(y=Power, x=ScreeningSampleSize
                            , colour=MedianSurvivalCase
                            , group=MedianSurvivalCase
                            )) +
            facet_grid(Var ~ . , scales="free_x") +
            geom_line() +
            geom_point(size=2.5 ) + 
            xlab("Screening Sample Size") +
            ylab("Power") +
            theme_bw() +
            ggtitle(plotTitle)
        return(p)
    }
})
