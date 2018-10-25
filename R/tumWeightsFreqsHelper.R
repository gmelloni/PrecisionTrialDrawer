##############################################################################
# HELPER FUNCTIONS FOR ANY METHOD THAT CONTAINS TUMOR.FREQS OR TUMOR.WEIGHTS #
##############################################################################

#---------------------------------------------------
# CHECK CONSISTENCY OF PARAMETER tumor.weights
#---------------------------------------------------
.tumor.weights.standardCheck <- function(tumor.weights 
                                         , tumor.freqs , object , tumor_type){
    if( any(is.na(tumor.weights)) || any(!is.numeric(tumor.weights)) ){
        stop("tumor.weights must be an integer vector")
    }
    if( any(names(tumor.weights) %notin% cpArguments(object)$tumor_type) ){
        stop("tumor.weights names don't match tumor_types in the object")
    }
  if( any(.superdup(names(tumor.weights))) ){
    stop("duplicated tumor.weights names")
  }
    if(any(tumor.weights<0)){
        stop("tumor.weights must be a numeric vector of positive values")
    }
    if(sum(tumor.weights)<1){
        stop("tumor.weights must contain at least one value > 0")
    }
    if(any(tumor.weights>100000)){
        message("Consider lower values for tumor.weights or use tumor.freqs")
    }
    if(!is.null(tumor.freqs)){
        stop("tumor.freqs and tumor.weights can't be set both at the same time")
    }
    if(!is.null(tumor_type)){
        if( any( sort(names(tumor.weights))!=sort(tumor_type) ) ){
            stop("tumor.weights and tumor_type selected do not match")
        }
    }

}

#---------------------------------------------------
# CHECK CONSISTENCY OF PARAMETER tumor.freqs
#---------------------------------------------------

.tumor.freqs.standardCheck <- function(tumor.weights , tumor.freqs 
                                       , object , tumor_type){
    if(!is.null(tumor.weights)){
        stop("tumor.freqs and tumor.weights can't be set both at the same time")
    }
    if( any(is.na(tumor.freqs)) || any(!is.numeric(tumor.freqs)) || 
        any(tumor.freqs<0 | tumor.freqs>1) ){
        stop("tumor.freqs must be an numeric vector between 0 and 1")
    }
    if( any(names(tumor.freqs) %notin% cpArguments(object)$tumor_type) ){
        stop("tumor.freqs names don't match tumor_types in the object")
    }
  if( any(.superdup(names(tumor.freqs))) ){
    stop("duplicated tumor.weights names")
  }
    if(!is.null(tumor_type)){
        if( any( sort(names(tumor.freqs))!=sort(tumor_type) ) ){
            stop("tumor.freqs and tumor_type selected do not match")
        }
    }
    if(!isTRUE(all.equal( sum(tumor.freqs) , 1))){
         stop("tumor.freqs must sum to 1")
    }
}


#-----------------------------------------------------
# TUMOR.WEIGHTS CALCULATOR
#-----------------------------------------------------

.tumor.weights.machine <- function(tumor.weights , mysamples , mydata ){
    # replacement in sample function is used only if 
    # we request more samples than the ones available
    useReplace <- if( any( lengths( mysamples[names(tumor.weights)] ) < 
                           tumor.weights ) ) TRUE else FALSE
    sampWeights <- lapply(names(tumor.weights) , function(x) {
        if(is.null(mysamples[[x]]) || tumor.weights[x]==0){
            return(NULL)
        }
        sample(mysamples[[x]] , size=tumor.weights[x] , replace=useReplace)
    })
    names(sampWeights) <- names(tumor.weights)
    # with the weights there is the possibility that 
    # the same sample is counted twice if replaced TRUE
    # this is a possibility we want to keep
    if(!useReplace){
        mydata <- mydata[ mydata$case_id %in% unlist(sampWeights) , ]
        mysamples <- sampWeights
        mysamples$all_tumors <- unlist(unname(sampWeights))
    } else {
        # every duplicated sample must become a new fresh sample
        sampWeights_New <- lapply(sampWeights , function(sampWeights_all){
            if(is.null(sampWeights_all)){
                return(NULL)
            }
            # empty vector that will contain the modified names
            # sampWeights_allNewNames <- c()
            sampWeights_allNewNames <- rep(NA , length(sampWeights_all))
            dups <- duplicated(sampWeights_all)
            # counter will be pasted to duplicated 
            # sample names to make them unique
            counter <- 0
            # using a while loop, we rename the sample up 
            # to the point there are no more duplicates
            # How it works? 
                # Start: c(1 , 2 , 3 , 1 , 1 , 2)
                # Round 1: c(1___0 , 2___0 , 3___0 , 1 , 1 , 2) we rename the 
                    # NON duplicated and store them in sampWeights_allNewNames
                # Round 1 remove: c(1 , 1 , 2)
                # Round 2: c(1___1 , 1 , 2___1)
                # Round 2 remove: c(1) no more dups, stop the while
                # recover the new vector: c(1___0 , 2___0 , 
                    # 3___0 , 1___1 , 2___1 , 1) now all names are unique
            while( any(dups) ){
                dups <- duplicated(sampWeights_all)
                newsamp <- paste(sampWeights_all[!dups] , counter , sep="____")
                if(counter==0){
                    sampWeights_allNewNames[seq_len(length(newsamp))] <- newsamp
                } else {
                    firstNA <- which.max(is.na(sampWeights_allNewNames))
                    sampWeights_allNewNames[
                      firstNA:(firstNA+length(newsamp)-1)] <- newsamp
                }
                sampWeights_all <- sampWeights_all[dups]
                counter <- counter + 1
            }
            if(any(is.na(sampWeights_allNewNames))){
                firstNA <- which.max(is.na(sampWeights_allNewNames))
                sampWeights_allNewNames[
                  firstNA:(firstNA+length(sampWeights_all)-1)] <- 
                  sampWeights_all
            }
            out <- strsplit(sampWeights_allNewNames , "____") %>% 
              vapply(. , '[' , character(1) , 1)
            names(out) <- sampWeights_allNewNames
            return(out)
        })
        sampWeights_New_all <- unlist(unname(sampWeights_New))
        mydata <- lapply( seq_len(length(sampWeights_New_all)) , function(x) {
            oldSamp <- sampWeights_New_all[x]
            if(is.null(oldSamp)){
                return(NULL)
            }
            newSamp <- names(sampWeights_New_all)[x]
            subData <- mydata[ mydata$case_id==oldSamp , ]
            if(nrow(subData)==0){
                return(NULL)
            }
            subData$case_id <- newSamp
            return(subData)
        # }) %>% do.call("rbind" , .)
        })
        mydata <- as.data.frame(data.table::rbindlist(mydata) 
                                , stringsAsFactors = FALSE)
        rownames(mydata) <- seq_len(nrow(mydata))
        mysamples <- lapply(sampWeights_New , names)
        names(mysamples) <- names(tumor.weights)
        mysamples$all_tumors <- unname(unlist(mysamples))
    }
    return(list(mydata , mysamples))
}


