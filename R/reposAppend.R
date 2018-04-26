setGeneric('appendRepo', function(object , repos) {
  standardGeneric('appendRepo')
})
setMethod('appendRepo', 'CancerPanel', function(object , repos){
  #init
  new_panel <- object
  
  ############################
  # CHECK panel is well formed
  # --------------------------
  # in case the user provies a panel with no @dataFull
  if(length(object@dataFull) == 0){
    stop("Nothing to append to: The CancerPanel you have entered has no data. Run getAlterations() to fetch the data.")
  }
  
  ############################
  # CHECK repos is well formed
  # --------------------------
  if(class(repos) != "list" ){stop("repos should be a list")}
  
  if(length(repos) != 4 ){stop("repos should be a list of length 4")}
  
  if(!all(names(repos) %in% c( "fusions", "mutations",  "copynumber", "expression"))){
    stop("repos is not correctly formatted. There are alteration types missing in the list")
  }
  
  if(!all(sapply(repos , function(x) names(x)) %in% c("data" , "Samples"))){
    stop("Each alteration type should contain two elements named data and Samples")
  }

  #specify format
  kColnamesIndataFull <- list(
    fusions=c("tumor_type","case_id","Gene_A","Gene_B","FusionPair","tier","frame")
    ,mutations=c("entrez_gene_id","gene_symbol","case_id","mutation_type"
                 ,"amino_acid_change","genetic_profile_id","tumor_type"
                 ,"amino_position","genomic_position")
    ,copynumber=c("gene_symbol","CNA","case_id","tumor_type","CNAvalue")
    ,expression=c("gene_symbol","expression","case_id","tumor_type","expressionValue")
  )
  
  ############################
  # APPEND Repo in panel
  # --------------------------
  
  #cycle through repos
  for(i in names(repos)){

    # if both data and Samples are NULL, skip
    if(is.null(repos[[i]]$Samples) & is.null(repos[[i]]$data)) {
      next
      # if data exists and samples is NULL, stop
    } else if( !is.null(repos[[i]]$data) & is.null(repos[[i]]$Samples) ) {
      stop(paste("In appending new data," , i , ", Samples are NULL"))
    } else {
      
      # --------------------------
      # Check over data consistency before append
      # --------------------------
      # data must be a dataframe (colnames check is performed by rbind)
      # Samples must be a list
      # elements of Samples must be character vectors
      # elements of Samples must be vectors of length > 1
      
      if(class(repos[[i]]$data)!="data.frame"){
        stop(paste("In appending new data," , i , "data is not a data.frame"))
      }
      
      if(class(repos[[i]]$Samples)!="list"){
        stop(paste("In appending new data," , i , "Samples is not a list"))
      }
      
      if(!all(sapply(repos[[i]]$Samples , class) %in% c("character" , "NULL"))){
        stop(paste("some sample names in" , i , "are not character vectors"))
      }
      # This check is useless if we know it's a character vector already
      # if(!all(lengths(repos[[i]]$Samples) > 0)){
      #   stop(paste("elements of Samples in", i ,"must be character vectors of length > 0"))
      # }
      # Last column in copynumber and expression is not necessary
      # In case missing, we add it manually
      if(i %in% c("copynumber" , "expression")){
        if( ncol(repos[[i]]$data)!=5 ){
          if(i=="copynumber"){
            repos[[i]]$data$CNAvalue <- rep(NA , nrow(repos[[i]]$data))
          } else {
            repos[[i]]$data$expressionValue <- rep(NA , nrow(repos[[i]]$data))
          }
        }
      }
      
      if(is.null(object@dataFull[[i]]$data)){
        # We can append data even if no pre existing data are present in dataFull
        # In this case, rbind cannot check for colnames and we have to add another control
        if(! identical(colnames(repos[[i]]$data) , kColnamesIndataFull[[i]])){
          colnameCollapse <- paste(kColnamesIndataFull[[i]] , collapse=", ")
          stop(paste("colnames in" 
                     , i 
                     , "do not match dataFull slot specifications. Compulsory colnames are:" 
                     , colnameCollapse))
        }
      }
      
      # --------------------------
      # Append data
      # --------------------------
      new_panel@dataFull[[i]]$data <- tryCatch({
        unique( rbind( object@dataFull[[i]]$data , repos[[i]]$data) ) 
        }, error=function(e) {
          stop(paste("colnames in" , i , "do not match dataFull slot specifications. Check str(cpData(object))"))
        })
      
      # --------------------------
      # Append sample names
      # --------------------------
      all_tumor_names <- unique( c( names(object@dataFull[[i]]$Samples) 
                                    , names(repos[[i]]$Samples)) )
      new_panel@dataFull[[i]]$Samples <- 
        lapply(all_tumor_names 
               , function(x) {
                    newsamps <- c(new_panel@dataFull[[i]]$Samples[[x]] , repos[[i]]$Samples[[x]])
                    uniquenewsamps <- unique(newsamps)
                    if(!identical(newsamps , uniquenewsamps)){
                      warning(paste("Duplicated sample names in" , i , x , "were unified"))
                    }
                    return(uniquenewsamps)
                    }
               )
      names(new_panel@dataFull[[i]]$Samples) <- all_tumor_names
    }
    #-----------------------------------------
    # CHECK TUMOR TYPE AND CASE_ID CONSISTENCY
    #
    # The user cannot put samples or tumor type in the data but not in the Samples
    if(any(new_panel@dataFull[[i]]$data$case_id %notin% unlist(new_panel@dataFull[[i]]$Samples))) {
      stop(paste( "In adding" , i , ", some case_id are in the data but not in the samples"))
    }
    if(any(new_panel@dataFull[[i]]$data$tumor_type %notin% names(new_panel@dataFull[[i]]$Samples))) {
      stop(paste( "In adding" , i , ", some tumor_type are in the data but are not listed in the samples"))
    }
  }
  ############################
  # APPEND NEW CANCER TYPES IN ARGUMENTS SLOT
  #---------------------------
  newcancertypes <- lapply(new_panel@dataFull , function(x) names(x$Samples)) %>% unlist %>% unique
  new_panel@arguments$tumor_type <- sort(newcancertypes)

  ############################
  # RETURNING
  # --------------------------
  new_panel <- subsetAlterations(new_panel)
  return(new_panel)
})