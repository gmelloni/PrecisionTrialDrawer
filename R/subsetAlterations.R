# Subset the downloaded alterations to match exactly with the panel specifications
# Add annotation on each single alteration if they are match with a drug or a user defined group
setGeneric('subsetAlterations', function(object 
                                         , rules=NULL
                                         ) {
	standardGeneric('subsetAlterations')
	})
setMethod('subsetAlterations', 'CancerPanel', function(object 
                                                       , rules=NULL
                                                       )
{
	if(length(object@dataFull)==0)
		stop("There are no data in this CancerPanel. Please, launch a getAlterations method before subsetting")
	panel <- object@arguments$panel
	kEmptyDF <- data.frame(drug=character(0)
							,group=character(0)
							,gene_symbol=character(0)
							,tumor_type=character(0)
							,case_id=character(0)
							,alteration_id=character(0)
							,stringsAsFactors=FALSE)
	#---------------
	# CHECK AND TRANSFORM RULES
	#---------------
	# Rules can be taken as argument or retrieved from the object
	if(is.null(rules)){
	  rules <- object@arguments$rules
	}
	# If we retrieve it from the object it can still be NULL
	if(!is.null(rules)){
	  # If we are checking a rules panel, we split the checks in two functions
	  # druggability is scorporated and contains only the cases in which an 
	  # entire drug is excluded/included from certain tumor types
	  druggabilityWhich <- which( apply(rules[ , c("gene_symbol" , "alteration" 
	                                           , "exact_alteration", "mutation_specification") 
	                                           , drop=FALSE] 
	                                    , 1 , function(x) all( x == "")))
	  if(length(druggabilityWhich)>0){
	    # If there are druggability rules, perform check
	    druggability <- rules[ druggabilityWhich , c("drug" , "group" , "tumor_type" , "in_out"), drop=FALSE]
	    # druggabilityCheck also turn include into exclude rules
	    druggability_full <- .druggabilityCheck(druggability , tumor_type = object@arguments$tumor_type)
	    exclude <- rules[ -druggabilityWhich , , drop=FALSE]
	  } else {
	    druggability <- NULL
	    druggability_full <- NULL
	    exclude <- rules
	  }
	  # Check on exclude panel (the one with 8 columns)
	  if(is.null(exclude) | nrow(exclude)==0){
	    excluded_full <- kEmptyDF[ , c("drug" , "group" , "case_id")]
	  } else {
	    exclude <- .panelCheck(exclude , comparison_panel=panel , tumor_type=object@arguments$tumor_type)
	    # panelCheck also transform positive rules into negative so the 8th column is actually no longer needed
	    excluded_full <- .excluder(object , exclude)
	  }
	  # Update slot rules
	  object@arguments[['rules']] <- rules
	} else {
	  exclude <- NULL
	  druggability <- NULL
	  excluded_full <- kEmptyDF[ , c("drug" , "group" , "case_id")]
	  druggability_full <- NULL
	  object@arguments['rules'] <- list(NULL)
	}

	# Subsetting fusions
  if( !is.null(object@dataFull$fusions$data) ) {
    message("Subsetting fusions...")
    if(nrow(object@dataFull$fusions$data)>0){
      fus_full <- object@dataFull$fusions$data
      panel_fus <- panel[ panel$alteration=="fusion" , ]
      myorder <- c("drug" , "group" , "FusionPair" , "tumor_type" , "case_id")
      fus_subset <- lapply(c("Gene_A" , "Gene_B" , "FusionPair") , function(x) {
        if(x %notin% c("Gene_A" , "Gene_B")){
          out <- merge(fus_full ,  panel_fus, by.x=x , by.y="gene_symbol")
        } else {
          out <- merge(panel_fus , fus_full , by.x="gene_symbol" , by.y=x)
        }
          out <- out[ , myorder]
          colnames(out)[colnames(out)=="FusionPair"] <- "gene_symbol"
          return(out)
        }) %>% do.call("rbind" , .)
      fus_subset$alteration_id <- paste("fus" , 1:nrow(fus_subset) , sep="_")
      if(nrow(excluded_full)>0){
        for( i in 1:nrow(excluded_full) ){
          samp <- excluded_full$case_id[i]
          drug <- excluded_full$drug[i]
          group <- excluded_full$group[i]
          fus_subset[ fus_subset$case_id == samp & fus_subset$drug == drug , "group"] <- group
          fus_subset[ fus_subset$case_id == samp & fus_subset$drug == drug , "drug"] <- "no_drug"
        }
      }
      if(!is.null(druggability_full)){
        for( i in 1:nrow(druggability_full)){
          drug <- druggability_full$drug[i]
          group <- druggability_full$group[i]
          tum_type <- druggability_full$tumor_type[i]
          fus_subset[ fus_subset$tumor_type == tum_type & fus_subset$drug == drug , "group"] <- group
          fus_subset[ fus_subset$tumor_type == tum_type & fus_subset$drug == drug , "drug"] <- "no_drug"
        }
      }
      object@dataSubset$fusions <- fus_subset
      } else {
        object@dataSubset[['fusions']] <- kEmptyDF
      }
    } else {
        object@dataSubset['fusions'] <- list(NULL)
    }
  # Subsetting mutations
  if( !is.null(object@dataFull$mutations$data)  ) {
    message("Subsetting mutations...")
    # Subsetting
  	if(nrow(object@dataFull$mutations$data)>0){
      muts_full <- object@dataFull$mutations$data
	    rs_df <- object@arguments$dbSNP_rs
	    muts_subset <- .subsetMutations(panel=panel , muts_full=muts_full , rs_df=rs_df)
    # Exclusion
  	  if(nrow(excluded_full)>0){
	  	  for( i in 1:nrow(excluded_full) ){
		     samp <- excluded_full$case_id[i]
			   drug <- excluded_full$drug[i]
			   group <- excluded_full$group[i]
			  muts_subset[ muts_subset$case_id == samp & muts_subset$drug == drug , "group"] <- group
			  muts_subset[ muts_subset$case_id == samp & muts_subset$drug == drug , "drug"] <- "no_drug"
	      }
  	  }
	    if(!is.null(druggability_full)){
	      for( i in 1:nrow(druggability_full)){
	        drug <- druggability_full$drug[i]
	        group <- druggability_full$group[i]
	        tum_type <- druggability_full$tumor_type[i]
	        muts_subset[ muts_subset$tumor_type == tum_type & muts_subset$drug == drug , "group"] <- group
	        muts_subset[ muts_subset$tumor_type == tum_type & muts_subset$drug == drug , "drug"] <- "no_drug"
	      }
	    }
	   object@dataSubset$mutations <- muts_subset
	  } else {
	   object@dataSubset[['mutations']] <- kEmptyDF
	 }
  } else {
    	object@dataSubset['mutations'] <- list(NULL)
  }
	# Subsetting CNVs
  if( !is.null(object@dataFull$copynumber$data) ) {
    message("Subsetting copynumber...")
    if(nrow(object@dataFull$copynumber$data)>0){
      cna_full <- object@dataFull$copynumber$data
			cna_subset <- .subsetCNA(panel=panel , cna_full=cna_full )
			if(nrow(excluded_full)>0){
			  for( i in 1:nrow(excluded_full) ){
			    samp <- excluded_full$case_id[i]
			    drug <- excluded_full$drug[i]
			    group <- excluded_full$group[i]
			    cna_subset[ cna_subset$case_id == samp & cna_subset$drug == drug , "group"] <- group
			    cna_subset[ cna_subset$case_id == samp & cna_subset$drug == drug , "drug"] <- "no_drug"
			  }
			}
			if(!is.null(druggability_full)){
			  for( i in 1:nrow(druggability_full)){
			    drug <- druggability_full$drug[i]
			    group <- druggability_full$group[i]
			    tum_type <- druggability_full$tumor_type[i]
			    cna_subset[ cna_subset$tumor_type == tum_type & cna_subset$drug == drug , "group"] <- group
			    cna_subset[ cna_subset$tumor_type == tum_type & cna_subset$drug == drug , "drug"] <- "no_drug"
			  }
			}
        	object@dataSubset$copynumber <- cna_subset
    	} else {
	    	object@dataSubset[['copynumber']] <- kEmptyDF
	    }
    } else {
    	object@dataSubset['copynumber'] <- list(NULL)
    }
	  # Subsetting expression
    if( !is.null(object@dataFull$expression$data) ) {
        message("Subsetting expression...")
    	if(nrow(object@dataFull$expression$data)>0){
        	expr_full <- object@dataFull$expression$data
			    expr_subset <- .subsetExpression(panel=panel , expr_full=expr_full)
			    if(nrow(excluded_full)>0){
			      for( i in 1:nrow(excluded_full) ){
			        samp <- excluded_full$case_id[i]
			        drug <- excluded_full$drug[i]
			        group <- excluded_full$group[i]
			        expr_subset[ expr_subset$case_id == samp & expr_subset$drug == drug , "group"] <- group
			        expr_subset[ expr_subset$case_id == samp & expr_subset$drug == drug , "drug"] <- "no_drug"
			      }
			    }
			    if(!is.null(druggability_full)){
			      for( i in 1:nrow(druggability_full)){
			        drug <- druggability_full$drug[i]
			        group <- druggability_full$group[i]
			        tum_type <- druggability_full$tumor_type[i]
			        expr_subset[ expr_subset$tumor_type == tum_type & expr_subset$drug == drug , "group"] <- group
			        expr_subset[ expr_subset$tumor_type == tum_type & expr_subset$drug == drug , "drug"] <- "no_drug"
			      }
			    }
        	object@dataSubset$expression <- expr_subset
    	} else {
	    	object@dataSubset[['expression']] <- kEmptyDF
	    }
    } else {
    	object@dataSubset['expression'] <- list(NULL)
    }
	  object@dataSubset$excluded <- excluded_full
  	return(object)
})