# Perform several checks on druggability rules
.druggabilityCheck <- function(panel , tumor_type=NULL) {
  # Since tibbles are now a thing, we don't require a strict data.frame, as long as data.frame
  # class is among the classes of the panel
  if("data.frame" %notin% class(panel))
    stop('The panel must be a data.frame')
  if(class(panel)[1]!="data.frame"){
    panel <- as.data.frame(panel , stringsAsFactors=FALSE)
  }
  ##################################################
  # PANEL DATAFRAME INPUT VALIDATION
  # ------------------------------------------------
  # check that columns are ok
  # ------------------------------------------------
  # Convert factors to characters
  panel <- .changeFactor(panel)
  colnames(panel) <- tolower(colnames(panel))
  # Convert all NAs into empty strings ""
  for(i in colnames(panel)){
    panel[ , i] <- .noNA(panel[ , i])
  }
  
  # check colnames
  expected_colnames <- c("drug" , "group" , "tumor_type" , "in_out")
  if( !identical(tolower(colnames(panel)) , expected_colnames) ) {
    print("Expected colnames: ")
    print(expected_colnames)
    print("Provided colnames: ")
    print(colnames(panel)[1:5])
    stop('The druggability panel dataframe is malformed.')
  }
  
  # ------------------------------------------------
  # Check drug column
  # ------------------------------------------------
  # check if drug names are uniquely identified, also for case sensitivity
  drugDups <- .superdup(panel$drug)
  drugDupsLower <- .superdup(tolower(panel$drug))
  if(!identical(drugDups , drugDupsLower)){
    mywarn <- paste(unique( panel$drug[xor(drugDups , drugDupsLower)]) , collapse=", ")
    warning("Some drug names seems to be case sensitive. Same drug but upper/lower case?:" %++% mywarn)
  }
  # Change drug names from "" to a reserved value "no_drug"
  if("no_drug" %in% panel$drug | "" %in% panel$drug ){
    warning(paste("'no_drug' is a reserved value for empty drug value."
                  , "Druggability panel is specifically designed to assign a drug to specific tumor types." 
                  , "Lines containing no_drug or empty values were removed" , sep="\n"))
    panel <- panel[ -which(panel$drug %in% c("no_drug" , "")) , , drop=FALSE]
  }
  
  # ------------------------------------------------
  # Check group column
  # ------------------------------------------------
  # Change group names from "" to a reserved value "no_group"
  if("no_group" %in% panel$group & any(panel$group=="")){
    warning(paste("'no_group' is a reserved value for empty group value."
                  , "Since empty values are not accepted, they have been changed to 'no_group'." 
                  , "Either call them all 'no_group' or leave all the cells empty" , sep="\n"))
    
  }
  panel$group[panel$group==""] <- "no_group"
  
  # ------------------------------------------------
  # Check in_out column
  # ------------------------------------------------
  # Only include and exclude are accepted
  if( ! all( panel$in_out %in% c("include" , "exclude") )){
    stop("Druggability panel column in_out can contain only include and exclude")
  }
  
  # ------------------------------------------------
  # Check rows are correctly formatted
  # ------------------------------------------------
  # look for non duplicated rows
  if(nrow(panel)!=nrow(unique(panel)))
    warning("There are duplicated rows that have been removed")
  panel <- unique(panel)
  
  # ------------------------------------------------
  # Check that rules do not contradict
  # ------------------------------------------------
  # Same drug - tumor_type is not allowed (cannot be both excluded and included)
  panelSub <- panel[ , c("drug" , "tumor_type")]
  if(any(duplicated(panelSub))){
    print(panel[ duplicated(panel) , ])
    stop("The same drug-tumor cannot be both included and excluded")
  }
  # The include rules cannot have conflicting groups (by drug)
  if( any( panel$in_out == "include" ) ){
    panelSub <- panel[ panel$in_out == "include" , c("drug" , "group") , drop=FALSE] %>% unique
    rulesPos <- sapply( unique(panelSub$drug) , function(x) nrow(unique(panelSub[ panelSub$drug==x , ])) ) 
    if(any(rulesPos != 1)){
      stop(paste("Panel druggability has conflicting include rules in the group column."
           , "Multiple include lines must have the same group entry for each drug" , sep="\n"))
    }
  }
  
  # Convert positive rules (include) into negative rules
  if(!is.null(tumor_type)){
    if(any(panel$in_out == "include")){
      panelPos <- panel[ panel$in_out == "include" , , drop=FALSE]
      transformPositives <- lapply( unique(panelPos$drug) , function(d) {
        positives <- panelPos[ panelPos$drug == d , , drop=FALSE]
        otherTumors <- setdiff( tumor_type , positives$tumor_type)
        if(length(otherTumors)==0){
          return(
            structure(
              list(drug=character()
                   ,group=character()
                   ,tumor_type=character()
                   ,in_out=character()) , class = "data.frame"))
        }
        negatives <- data.frame(
          drug=rep(positives$drug[1] , length(otherTumors))
          , group=rep(positives$group[1] , length(otherTumors))
          , tumor_type=otherTumors
          , in_out=rep("exclude" , length(otherTumors))
          , stringsAsFactors=FALSE)
        return(negatives)
      }) %>% do.call("rbind" , .)
      panel <- rbind( transformPositives , panel[ panel$in_out != "include" , ] )
    }
  }
  return(panel)
}
  
# Perform several checks on panel format
.panelCheck <- function(panel , comparison_panel=NULL , tumor_type=NULL) {
  # Since tibbles are now a thing, we don't require a strict data.frame, as long as data.frame
  # class is among the classes of the panel
  if("data.frame" %notin% class(panel))
    stop('The panel must be a data.frame')
  if(class(panel)[1]!="data.frame"){
    panel <- as.data.frame(panel , stringsAsFactors=FALSE)
  }
  # The presence of comparison panel trigger the signal that .panelCheck
  # is in the presence of a rules panel and not the main panel
  if(!is.null(comparison_panel)){
    forExclude <- TRUE
  } else {
    forExclude <- FALSE
  }
  ##################################################
  # PANEL DATAFRAME INPUT VALIDATION
  # ------------------------------------------------
  # check that columns are ok
  # ------------------------------------------------
  # Convert factors to characters
  panel <- .changeFactor(panel)
  colnames(panel) <- tolower(colnames(panel))
  # Convert all NAs into empty strings ""
  for(i in colnames(panel)){
    panel[ , i] <- .noNA(panel[ , i])
  }
  
  ### DEPRECATED. WE WANT TO BE MORE STRICT
  # Column 6th is optional. We have to check if there is the 6th column, the
  # colnames and if missing, add it.
  # input_ncol <- ncol(panel)
  # # Check 6th column
  # if(input_ncol %in% c(6,7,8)) {
  #   if( !identical( colnames(panel)[6] , c("group") ) ) {
  #     stop('The sixth column name should be group')
  #   }
  # }
  # # If we have only 5 columns...add the 6th with a fixed value
  # if(input_ncol==5) {
  #   panel$group <- "no_group"
  #   # If it's the dataframe for rules add also the tumor_type column
  #   if(forExclude){
  #     panel$tumor_type <- ""
  #   }
  # } else if(input_ncol==6) {
  #     if(forExclude){
  #       panel$tumor_type <- ""
  #     }
  # } else {
  #   if(forExclude){
  #     if( !identical( colnames(panel)[7] , c("tumor_type") ) ) {
  #       stop('The seventh column name should be tumor_type')
  #     }
  #   }
  # }
  # 09/03/2018 columns are now mandatory
  if(forExclude){
    # check colnames
    expected_colnames <- c("drug" , "gene_symbol" , "alteration" 
                           , "exact_alteration", "mutation_specification"
                           , "group" , "tumor_type" , "in_out")
  } else {
    expected_colnames <- c("drug" , "gene_symbol" , "alteration" 
                           , "exact_alteration", "mutation_specification"
                           , "group")
  }
  if( !identical(tolower(colnames(panel)) , expected_colnames) ) {
    print("Expected colnames: ")
    print(expected_colnames)
    print("Provided colnames: ")
    print(colnames(panel)[1:5])
    stop('The rules panel dataframe is malformed')
  }
  
  # ------------------------------------------------
  # Check drug column
  # ------------------------------------------------
  # check if drug names are uniquely identified, also for case sensitivity
  drugDups <- .superdup(panel$drug)
  drugDupsLower <- .superdup(tolower(panel$drug))
  if(!identical(drugDups , drugDupsLower)){
    mywarn <- paste(unique( panel$drug[xor(drugDups , drugDupsLower)]) , collapse=", ")
    warning("Some drug names seems to be case sensitive. Same drug but upper/lower case?:" %++% mywarn)
  }
  # Change drug names from "" to a reserved value "no_drug"
  if("no_drug" %in% panel$drug & any(panel$drug=="")){
    warning(paste("'no_drug' is a reserved value for empty group value."
                  , "Since empty values are not accepted, they have been changed to 'no_drug'." 
                  , "Either call them all 'no_drug' or leave all the cell empty" , sep="\n"))
  } else {
    if(any(panel$drug=="")){
      panel$drug[panel$drug==""] <- "no_drug"
    }
  }
  # ------------------------------------------------
  # Check group column
  # ------------------------------------------------
  # Change group names from "" to a reserved value "no_group"
  if("no_group" %in% panel$group & any(panel$group=="")){
    warning(paste("'no_group' is a reserved value for empty group value."
                  , "Since empty values are not accepted, they have been changed to 'no_group'." 
                  , "Either call them all 'no_group' or leave all the cell empty" , sep="\n"))
  } else {
    if(any(panel$group=="")){
      panel$group[panel$group==""] <- "no_group"
    }
  }
  
  # ------------------------------------------------
  # Check in_out column
  # ------------------------------------------------
  # Only include and exclude are accepted
  if( ! all( panel$in_out %in% c("include" , "exclude") )){
    stop("Druggability panel column in_out can contain only include and exclude")
  }
  
  # ------------------------------------------------
  # Check rows are correctly formatted
  # ------------------------------------------------
  # look for non duplicated rows
  if( any(duplicated(panel)) )
    warning("There are duplicated rows that have been removed")
  panel <- unique(panel)
  # check gene presence at every row
  if(any( panel$gene_symbol=="" | panel$gene_symbol=="."))
    stop("The panel must contain a gene for every row")
  # Put every gene in uppercase. This will help a lot in retrieving mutations
  # panel$gene_symbol <- toupper(panel$gene_symbol)
  # Trim trailing spaces in gene symbols
  panel$gene_symbol <- .myTrimmer(panel$gene_symbol)
  
  # ------------------------------------------------
  # Check content is correctly formatted
  # ------------------------------------------------
  # check alteration format
  if(!all( panel$alteration %in% c("SNV" , "CNA" , "expression" , "fusion")))
    stop('The panel can accept only SNV, CNA fusion or expression in the alteration column')
  
  # check CNA format
  cna_allowed_values <- c("amplification" , "deletion" , "normal" , "")
  idx <- which(panel$alteration=="CNA")
  cna_comparison <- panel[ idx , "exact_alteration"] %in% cna_allowed_values
  if(!all(cna_comparison)) {
    message("Copy number alterations at lines" %++% 
              paste(idx[!cna_comparison] , collapse=", ") %++% 
              "are written in the wrong way:\n")
    print(panel[ idx[!cna_comparison] , ])
    message("The accepted values for alteration and exact alterations are:\n")
    print(data.frame(alteration=rep("CNA" , length(cna_allowed_values))
                     ,exact_alteration=cna_allowed_values))
    stop("Panel is malformed")
  }
  
  # check SNV format
  snv_allowed_values <- c(""
                          ,"mutation_type"
                          ,"amino_acid_position"
                          ,"amino_acid_variant"
                          ,"genomic_position"
                          ,"genomic_variant"
                          ,"dbSNP_rs")
  idx2 <- which(panel$alteration=="SNV")
  snv_comparison <- panel[ idx2 , "exact_alteration"] %in% snv_allowed_values
  if(!all(snv_comparison)) {
    message("Mutation alterations at lines" %++% 
              paste(idx2[!snv_comparison] , collapse=", ") %++% 
              "are written in the wrong way:\n")
    print(panel[ idx2[!snv_comparison] , ])
    message("The accepted values for alteration and exact_alterations are:\n")
    print(data.frame(alteration=rep("SNV" , length(snv_allowed_values))
                     ,exact_alteration=snv_allowed_values))
    stop("Panel is malformed")
  }
  
  # TO DO: a deep check on alteration format for mutations
  # shallow check for exact_alteration/mutation_specification formats
  # make sure that for each row, if either one of the two column values is specified, 
  # also the value of the other column is.
  muts_specs_check <- (panel[idx2 ,"exact_alteration"]=="" & panel[idx2 ,"mutation_specification"]=="") |
    (panel[idx2 ,"exact_alteration"]!="" & panel[idx2 ,"mutation_specification"]!="")
  if(!all(muts_specs_check)) {
    message("Mutation alterations at lines" %++% 
              paste(idx2[!muts_specs_check] , collapse=", ") %++% 
              "are written in the wrong way:\n")
    print(panel[ idx2[!muts_specs_check] , ])
    message("If exact_alteration is empty, mutation_specification must be empty too\n")
    stop("Panel is malformed")
  }
  
  # Check expression
  expression_allowed_values <- c("up" , "normal" , "down" , "")
  idx3 <- which(panel$alteration=="expression")
  expression_comparison <- panel[ idx3 , "exact_alteration"] %in% expression_allowed_values
  if(!all(expression_comparison)) {
    message("Mutation alterations at lines" %++% 
              paste(idx3[!expression_comparison] , collapse=", ") %++% 
              "are written in the wrong way:\n")
    print(panel[ idx3[!expression_comparison] , ])
    message("The accepted values for alteration and exact_alterations are:\n")
    print(data.frame(alteration=rep("expression" , length(expression_allowed_values))
                     ,exact_alteration=expression_allowed_values))
    stop("Panel is malformed")
  }
  
  # If this is an exclusion panel for drug resistance
  # Check that every gene in the exclusion panel is also present in the main panel
  # The missing genes throw a warning, rather than a stop because the rules potentially can involve 
  # genes that are not in the panel
  if(forExclude){
    if(all(panel$mutation_type != "fusion")){
      if( !all( panel$gene_symbol %in% comparison_panel$gene_symbol)){
        warning("The exclusion panel contains genes that are not present in the main panel")
      }
    } else {
      genepanel <- strsplit(panel$gene_symbol , split = "__") %>% unlist
      genecomparison <- strsplit(comparison_panel$gene_symbol , split = "__") %>% unlist
      if( !all(  genepanel %in% genecomparison )){
        warning("The exclusion panel contains genes that are not present in the main panel")
      }
    }
  }
  
  #--------------------------------
  # Convert positive rules (include) into negative rules
  #--------------------------------
  if(forExclude & !is.null(tumor_type)){
    if(any(panel$in_out == "include")){
      panelPos <- panel[ panel$in_out == "include" , , drop=FALSE]
      panelPos_split <- split(panelPos 
                              , f = paste(
                                panelPos$drug
                                ,panelPos$gene_symbol
                                ,panelPos$alteration
                                ,panelPos$exact_alteration
                                ,panelPos$mutation_specification
                                ,panelPos$group
                                ,panelPos$tumor_type
                                ,panelPos$in_out))
      transformPositives <- lapply(panelPos_split , function(df){
        myrule <- lapply( unique(df$drug) , function(d) {
          positives <- df[ df$drug == d , , drop=FALSE]
          otherTumors <- setdiff( tumor_type , positives$tumor_type)
          if(length(otherTumors)==0){
              negatives <- structure(
                list(drug=character()
                ,gene_symbol=character()
                ,alteration=character()
                ,exact_alteration=character()
                ,mutation_specification=character()
                ,group=character()
                ,tumor_type=character()
                ,in_out=character()) , class = "data.frame")
          } else {
            negatives <- data.frame(
              drug=rep(positives$drug[1] , length(otherTumors))
              , gene_symbol=rep(positives$gene_symbol[1] , length(otherTumors))
              , alteration=rep(positives$alteration[1] , length(otherTumors))
              , exact_alteration=rep(positives$exact_alteration[1] , length(otherTumors))
              , mutation_specification=rep(positives$mutation_specification[1] , length(otherTumors))
              , group=rep(positives$group[1] , length(otherTumors))
              , tumor_type=otherTumors
              , in_out=rep("exclude" , length(otherTumors))
              , stringsAsFactors=FALSE)
          }
          return(negatives)
          }) %>% do.call("rbind" , .)
        return(myrule)
      }) %>% do.call("rbind" , .)
      panel <- rbind( transformPositives , panel[ panel$in_out != "include" , ] )
    }
  }
  return(panel)
}

# This function creates a data.frame with sample names, drug they should be excluded for and new group assigned
.excluder <- function(object 
                      , exclude
                      ) {
  emptyExcludedDF <- data.frame(drug=character(0)
                                ,group=character(0)
                                ,case_id=character(0)
                                ,stringsAsFactors=FALSE)
  # Derive the set of tumor types to apply exclude
  # "" is equivalent to ALL tumor types
  exclude_tums <- unique(exclude$tumor_type)
  applicableTumors <- intersect(exclude_tums , object@arguments$tumor_type)
  applicableTumors <- c(applicableTumors , if("" %in% exclude_tums) "" else NULL)
  if(length(applicableTumors)==0){
    return(emptyExcludedDF)
  }
  
  byTumor <- lapply( applicableTumors , function(tum) {
    if(tum == ""){
      tum <- object@arguments$tumor_type
    }
    exclude_fus <- exclude[ exclude$alteration=="fusion" , ]  %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    fus_full <- object@dataFull$fusions$data %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    if( !is.null(fus_full) & nrow(fus_full)>0 & nrow(exclude_fus)>0) {
      myorder <- c("drug" , "group" , "FusionPair" , "tumor_type" , "case_id")
      fus_subset <- lapply(c("Gene_A" , "Gene_B" , "FusionPair") , function(x) {
        if(x %notin% c("Gene_A" , "Gene_B")){
          out <- merge(fus_full ,  exclude_fus, by.x=x , by.y="gene_symbol")
        } else {
          out <- merge(exclude_fus , fus_full , by.x="gene_symbol" , by.y=x)
        }
        out <- out[ , myorder]
        colnames(out)[colnames(out)=="FusionPair"] <- "gene_symbol"
        return(out)
      }) %>% do.call("rbind" , .)
      excluded_fus <- unique(fus_subset[ , c("drug" , "group" , "case_id")])
    } else {
      excluded_fus <- emptyExcludedDF
    }
    exclude_mut <- exclude[ exclude$alteration=="SNV" , ]  %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    muts_full <- object@dataFull$mutations$data %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    if( !is.null(muts_full) & nrow(muts_full)>0 & nrow(exclude_mut)>0) {
      rs_df <- object@arguments$dbSNP_rs
      muts_subset <- .subsetMutations(panel=exclude_mut , muts_full=muts_full , rs_df=rs_df)
      excluded_mut <- unique(muts_subset[ , c("drug" , "group" , "case_id")])
    } else {
      excluded_mut <- emptyExcludedDF
    }
    
    exclude_cna <- exclude[ exclude$alteration=="CNA" , ] %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    cna_full <- object@dataFull$copynumber$data %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    if( !is.null(cna_full) & nrow(cna_full)>0 & nrow(exclude_cna)>0) {
      cna_subset <- .subsetCNA(panel=exclude_cna , cna_full=cna_full)
      excluded_cna <- unique(cna_subset[ , c("drug" , "group" , "case_id")])
    } else {
      excluded_cna <- emptyExcludedDF
    }
    
    exclude_expr <- exclude[ exclude$alteration=="expression" , ] %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    expr_full <- object@dataFull$expression$data %>% .[ .$tumor_type %in% tum , , drop=FALSE]
    if( !is.null(expr_full) & nrow(expr_full)>0 & nrow(exclude_expr)>0) {
      expr_subset <- .subsetExpression(panel=exclude_expr , expr_full=expr_full)
      excluded_expr <- unique(expr_subset[ , c("drug" , "group" , "case_id")])
    } else {
      excluded_expr <- emptyExcludedDF
    }
    excluded_full <- unique( rbind(excluded_fus , excluded_mut , excluded_cna , excluded_expr) )
    return(excluded_full)
    }) %>% do.call("rbind" , .)
    return(byTumor)
}