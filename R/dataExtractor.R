
#-------------------------------------------------------------
# Extract requested data from object under user specification
#-------------------------------------------------------------

# This method is used by most of the methods of the package
# For compatibility with the rest of the code, it requires an environment in which variables are stored
# Then variables are taken with a get outside of this method and the environment destroyed if myenv is provided
setGeneric('dataExtractor', function(object , alterationType=c("copynumber" , "expression" , "mutations" , "fusions") 
                        , tumor_type=NULL , collapseMutationByGene=TRUE 
                        , collapseByGene=FALSE 
                        , myenv=NULL , tumor.weights=NULL) {
  standardGeneric('dataExtractor')
  })
setMethod('dataExtractor', 'CancerPanel', function(object , alterationType=c("copynumber" , "expression" , "mutations" , "fusions") 
                        , tumor_type=NULL , collapseMutationByGene=TRUE 
                        , collapseByGene=FALSE 
                        , myenv=NULL , tumor.weights=NULL)
{
    #------------------------------------
    # Check parameters
    #------------------------------------
    possibleAlterations <- c("copynumber" , "expression" , "mutations" , "fusions")
    if( any(alterationType %notin% possibleAlterations) ){
        stop("alterationType can only be one or more of the following" %++% paste(possibleAlterations , collapse=", "))
      }
    if( !is.null(myenv) && !is.environment(myenv) ){
        stop("myenv must be an environment")
    }
    if( !is.logical(collapseByGene) ){
        stop("collapseByGene must be a logical")
    }
    if(!is.logical(collapseMutationByGene)) {
        stop("collapseMutationByGene must be a logical")
    }
    if(!is(object , "CancerPanel")) {
        stop("object must be an instance of class CancerPanel")
    }
    if(!is.null(tumor_type) && !is.character(tumor_type)) {
        stop("tumor_type must be a character vector containing cancer types")
    }
    if(!is.null(tumor_type) && any(tumor_type %notin% object@arguments$tumor_type)){
      stop("The following tumor types are not in this object:" %++% 
             paste(setdiff(tumor_type , object@arguments$tumor_type) , collapse=", "))
    }
    #------------------------------------
    # Bind dataframes according to alteratioType specified
    #------------------------------------
    if(identical(object@dataSubset , list())){
        stop("dataSubset slot is empty. run subsetAlterations")
    }
    toBeBind <- object@dataSubset[alterationType]
    checkDataExistance <- sapply(toBeBind , is.null)
    if(any(checkDataExistance)){
        warning("No data for" %++% paste(names(toBeBind)[checkDataExistance] , collapse=", ") %+% ". They will be removed from alterationType")
        alterationType <- alterationType[!checkDataExistance]
    }
    # Define the dataset
    # mydata <- do.call("rbind" , toBeBind)
    mydata <- as.data.frame(data.table::rbindlist(toBeBind) , stringsAsFactors=FALSE)
    # Reduce the dataset is tumor_type is set
    if(!is.null(tumor_type)){
        mydata <- mydata[ mydata$tumor_type %in% tumor_type , ]
        if(nrow(mydata)==0){
            stop("No data to visualize for selected tumor type")
        }
    }
    
    #------------------------------------
    # Define reference set of samples
    #------------------------------------
    # If just one alteration type is defined, the data are stored within the slot of the alteration type
    if(length(alterationType)>1) {
        all_samples <- lapply(alterationType , function(x) object@dataFull[[x]]$Samples)
        names(all_samples) <- alterationType
        available_tumor_types <- unique(lapply(all_samples , names) %>% unlist)
        mysamples <- lapply(available_tumor_types , function(x) {
                            # browser()
                            out <- Reduce(intersect , lapply(all_samples , '[[' , x) ) %>% unique
                            if(identical(character(0) , out))
                                out <- NULL
                            return(out)
                        })
        names(mysamples) <- available_tumor_types
        if(!is.null(tumor_type)){
            if(all(tumor_type %in% names(mysamples))){
                mysamples <- mysamples[tumor_type]    
            } else {
                stop("No available samples for the tumor types selected or alteration type selected")
            }
            
        }
    } else {
        mysamples <- object@dataFull[[alterationType]][["Samples"]]
        mysamples <- lapply(mysamples , unique)
        if(!is.null(tumor_type)){
            if(all(tumor_type %in% names(mysamples))){
                mysamples <- mysamples[tumor_type]    
            } else {
                stop("No available samples for the tumor types selected or alteration type selected")
            }
            
        }
    }
    mysamples$all_tumors <- unlist(mysamples) %>% unname %>% unique
    # Reduce the dataset to only the patients included in mysamples
    mydata <- mydata[ mydata$case_id %in% mysamples$all_tumors , ]
    
    #-----------------------------------------
    # Uniquify if collapse options are set
    #-----------------------------------------
    mydata$alteration_id <- strsplit(mydata$alteration_id , "_") %>% sapply(. , '[' , 1)
    # Raise a warning if the tumor types present in the object differ from mydata$tumor_type
    if(is.null(tumor_type)){
        tum_type_diff <- setdiff(object@arguments$tumor_type , unique(mydata$tumor_type))
    } else {
        # If tumor_type is specified, report the difference between the tumors requested and present
        tum_type_diff <- setdiff(tumor_type , unique(mydata$tumor_type))
    }
    if(length(tum_type_diff)!=0){
        warning("The following tumor types have no alteration to display:" %++% paste(tum_type_diff , collapse=", "))
    }
    # Collapsing Mutation by gene
    # This option is valid for mutation. If a gene is mutated in more than 1 spot in a patient
    # Does it count for 1 or more alterations?
    if(collapseMutationByGene) {
        mydata <- unique(mydata)
    }
    # Collapsing by gene
    # This option is valid for all alterations.
    # If a gene is altered in multiple ways (e.g. both amplified and mutated) on the same patient
    # Does it count for 1 or more alterations?
    if(collapseByGene) {
        mydata <- unique(mydata[ , colnames(mydata)!="alteration_id"])
    }
    if(nrow(mydata)==0){
        allcombs <- lapply( seq.int(1 , length(alterationType) , 1) , function(x) {
            combn(alterationType , x , simplify=FALSE)
        }) %>% unlist(recursive = FALSE)
        mytums <- object@arguments$tumor_type
        sampSummary <-   lapply( allcombs , function(comb) {
          allsamps <- lapply( object@dataFull[comb] , '[[' , 'Samples')
          sapply( mytums , function(tum) {
              length(Reduce("intersect" , lapply(allsamps , '[[' , tum)))
          })
        }) %>% do.call("rbind" , .)
        sampSummary <- cbind( Combinations = sapply( allcombs , function(x) paste(x , collapse=","))
                              , sampSummary)
        message(paste0(capture.output(sampSummary), collapse = "\n"))
        stop("No alteration to display for the selected tumor types and alteration types: check the output above to see what is available")
    }

    #------------------------------------
    # Apply tumor.weights
    #------------------------------------
    if(!is.null(tumor.weights)){
        newDataAndSamps <- .tumor.weights.machine(tumor.weights , mysamples , mydata )
        mydata <- newDataAndSamps[[1]]
        mysamples <- newDataAndSamps[[2]]
      }
    #------------------------------------
    # Assign variable to the environment
    #------------------------------------
    if(is.null(myenv)){
        return(list(data=mydata , Samples=mysamples , tumor_not_present=tum_type_diff))
    } else {
        myenv$mydata <- mydata
        myenv$mysamples <- mysamples
        myenv$tum_type_diff <- tum_type_diff
    }
})

