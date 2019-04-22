#-------------------------------
# CHECK CUSTOM DATA CONSISTENCY
#
.checkRepos <- function(repos) {
  if (!is.list(repos)) {
    stop("repos should be a list")
  }
  if (length(repos) != 4) {
    stop("repos should be a list of length 4")
  }
  if (!all(names(repos) %in% c("fusions", "mutations"
                               ,  "copynumber", "expression"))) {
    stop("repos is not correctly formatted. alteration type missing")
  }
  reposNames <- lapply(repos , function(x)
    names(x))
  if (!all(vapply(reposNames , function(k)
    identical(k , c("data" , "Samples")) , logical(1)))) {
    stop(paste("Each alteration type should contain"
               , "two elements named data and Samples"))
  }
  kColnamesIndataFull <- list(
    fusions = c(
      "tumor_type"
      ,"case_id"
      ,"Gene_A"
      ,"Gene_B"
      ,"FusionPair"
      ,"tier"
      ,"frame"
    ),mutations = c(
      "entrez_gene_id"
      ,"gene_symbol"
      ,"case_id"
      ,"mutation_type"
      ,"amino_acid_change"
      ,"genetic_profile_id"
      ,"tumor_type"
      ,"amino_position"
      ,"genomic_position")
    ,copynumber = c("gene_symbol", "CNA", "case_id", "tumor_type", "CNAvalue")
    ,expression = c(
      "gene_symbol"
      ,"expression"
      ,"case_id"
      ,"tumor_type"
      ,"expressionValue")
  )
  for (i in names(kColnamesIndataFull)) {
    # if both data and Samples are NULL, skip
    if (is.null(repos[[i]]$Samples) & is.null(repos[[i]]$data)) {
      next
      # if data exists and samples is NULL, stop
    } else if (!is.null(repos[[i]]$data) &
               is.null(repos[[i]]$Samples)) {
      stop(paste("In" , i , ", Samples are NULL"))
    } else {
      # --------------------------
      # Check over data consistency before append
      # --------------------------
      # data must be a dataframe
      # colnames must be consistent
      # Samples must be a list
      # elements of Samples must be character vectors
      # elements of Samples must be vectors of length > 1
      
      if (!is.data.frame(repos[[i]]$data)) {
        stop(paste("In" , i , "data is not a data.frame"))
      }
      # Last column in copynumber and expression is not necessary
      # In case missing, we add it manually
      if (i %in% c("copynumber" , "expression")) {
        if (ncol(repos[[i]]$data) != 5) {
          if (i == "copynumber") {
            repos[[i]]$data$CNAvalue <- rep(NA , nrow(repos[[i]]$data))
          } else {
            repos[[i]]$data$expressionValue <- rep(NA , nrow(repos[[i]]$data))
          }
        }
      }
      if (!identical(colnames(repos[[i]]$data) ,  kColnamesIndataFull[[i]])) {
        stop(paste("colnames in repos are wrong for the alteration type" , i))
      }
      
      if (!is.list(repos[[i]]$Samples)) {
        stop(paste("In data," , i , ", Samples is not a list"))
      }
      
      if (!all(vapply(repos[[i]]$Samples , class 
                      , character(1)) == "character")) {
        stop(paste("some sample names in" , i , "are not character vectors"))
      }
      
      if (!all(lengths(repos[[i]]$Samples) > 0)) {
        stop(paste(
          "elements of Samples in",
          i ,
          "must be character vectors of length > 0"
        ))
      }
      #-----------------------------------------
      # CHECK TUMOR TYPE AND CASE_ID CONSISTENCY
      #
      # The user cannot put samples or tumor type 
      # in the data but not in the Samples
      if (any(repos[[i]]$data$case_id %notin% unlist(repos[[i]]$Samples))) {
        stop(paste(
          "In" ,
          i ,
          ", some case_id are in the data but not in the samples"
        ))
      }
      if (any(repos[[i]]$data$tumor_type %notin% names(repos[[i]]$Samples))) {
        stop(paste(
          "In" ,
          i ,
          ", some tumor_type are in the data but are not listed in the samples"
        ))
      }
    }
  }
}
################################################################################
# CLASS METHODS DECLARATION
################################################################################
# Download the alterations of interest
setGeneric('getAlterations', function(object
                ,tumor_type = NULL
                ,repos = NULL
                ,mutation_type = c("all_nonsynonymous"
                    ,"all_mutations"
                    ,"missense"
                    ,"truncating")
                ,expr_z_score = 2
                ,BPPARAM = bpparam("SerialParam")
                ,gene_block = 50){
  standardGeneric('getAlterations')
})
setMethod('getAlterations', 'CancerPanel', function(object
                ,tumor_type = NULL
                ,repos = NULL
                ,mutation_type = c("all_nonsynonymous"
                    ,"all_mutations"
                    ,"missense"
                    ,"truncating")
                ,expr_z_score = 2
                ,BPPARAM = bpparam("SerialParam")
                ,gene_block = 50)
{
  if (is.null(object)) {
    stop("No CancerPanel object provided")
  }
  ##############################################################################
  # IF repos is provided, get the data from it into the class
  #-----------------------------------------------------------------------------
  # The patients sets are retrieved from the data
  # Freq is not calculated
  if (!is.null(repos)) {
    .checkRepos(repos)
    object@dataFull$fusions <- list(data = repos$fusions$data)
    object@dataFull$mutations <- list(data = repos$mutations$data)
    object@dataFull$copynumber <-
      list(data = repos$copynumber$data)
    object@dataFull$expression <-
      list(data = repos$expression$data)
    tumor_type <- c()
    for (i in c("mutations" , "copynumber" , "expression" , "fusions")) {
      # If no Samples are provided, they are calculated from the data
      if (is.null(repos[[i]]$Samples)) {
        if (!is.null(repos[[i]])) {
          pats <- lapply(unique(repos[[i]]$tumor_type) , function(x) {
            repos[[i]][repos[[i]]$tumor_type == x , "case_id"] %>%
              unique
          })
          names(pats) <- unique(repos[[i]]$tumor_type)
          tumor_type <-
            unique(c(tumor_type , unique(repos[[i]]$tumor_type)))
          object@dataFull[[i]]$Samples <- pats
        } else {
          object@dataFull[[i]]['Samples'] <- list(NULL)
        }
      } else {
        object@dataFull[[i]]['Samples'] <- list(repos[[i]]$Samples)
      }
    }
    object@arguments$tumor_type <-
      lapply(object@dataFull , function(x) {
        if (!is.null(x$Samples))
          names(x$Samples)
      }) %>% unlist %>% unique
    return(object)
  }
  
  ##############################################################################
  # If no custom data are provided, you should specify tumor types
  #-----------------------------------------------------------------------------
  #Initial checks
  
  if (is.null(tumor_type) || tumor_type == "") {
    stop("You must specify one or more tumor types or tumor studies")
  }
  # Check consistency of tumor_type parameter
  if ("all_tumors" %in% tumor_type & length(tumor_type) > 1) {
    warning("tumor_type was set to all_tumor")
    tumor_type <- "all_tumors"
  }
  # Check format of expr_z_score
  if (!is.numeric(expr_z_score)) {
    stop("Expression z score should be a numeric value (e.g. 2 or 3)")
  }
  
  if (!is.numeric(gene_block)) {
    stop("gene_block must be a number between 1 and 100")
  } else {
    if (gene_block < 1 | gene_block > 100) {
      stop("gene_block must be a number between 1 and 100")
    }
  }
  
  #-----------------------------------------------------------------------------
  # Integrity check: All tumor types and tumor studies must be in cbioportal
  # You cannot mix tumor studies and tumor types
  #-----------------------------------------------------------------------------
  if (!is.null(tumor_type) && tumor_type[1] != "all_tumors") {
    #get list of tumor type available
    showTum <- showTumorType()
    availtumtype <- showTum$tumor_type
    names(availtumtype) <- showTum$name
    #get list of cancer studies available
    availcancerstudy <- showCancerStudy()[, 1]
    
    #Looking for any missing match. If a missmatch is found, print an error
    if (any(tumor_type %notin% availtumtype) &
        any(tumor_type %notin% availcancerstudy)) {
      notavail <- setdiff(tumor_type , c(availtumtype , availcancerstudy))
      if (length(notavail) == 0)
        stop("You probably mixed cancer studies and tumor types")
      else
        stop(
          paste(
            "The following tumor types or cancer studies are not available:" ,
            paste(notavail , collapse = ", ") ,
            ". Check with showTumorType() or showCancerStudy()"
          )
        )
    }
  }
  
  ##############################################################################
  # DOWNLOAD DATA
  #-----------------------------------------------------------------------------
  noDataMex <- paste("\nNo mutation data available for the specified" 
      ,"tumor types or cancer studies...")
  if (tumor_type[1] == "all_tumors") {
    derived_tumor_type <- showTumorType()$tumor_type %>% unique
  } else {
    derived_tumor_type <-
      unname(unique(vapply(tumor_type, function(x)
        strsplit(x , "_")[[1]][1] , character(1))))
  }
  object@arguments$tumor_type <- derived_tumor_type
  panel <- cpArguments(object)$panel
  
  #-----------------------------------------------------------------------------
  # LOOK FOR FUSIONS
  #-----------------------------------------------------------------------------
  if (any(cpArguments(object)$panel$alteration == "fusion")) {
    message("Retrieving fusions...")
    fus <-
      readRDS(file.path(
        system.file(package = "PrecisionTrialDrawer") ,
        "extdata" ,
        "TCGAtranslocation.rds"
      ))
    fusSamples <-
      readRDS(file.path(
        system.file(package = "PrecisionTrialDrawer") ,
        "extdata" ,
        "TCGAtranslocationSamples.rds"
      ))
    if (tumor_type[1] != "all_tumors") {
      # fusion table does not include the study so we include 
      # just the tumor types of the study (e.g. brca_tcga_pub becomes brca)
      fus <- fus[fus$tumor_type %in% derived_tumor_type ,]
    }
    if (nrow(fus) == 0) {
      message("No fusion data for the specified tumor types")
      object@dataFull$fusions <- list(data = NULL
                                      , Samples = NULL)
    } else {
      samples_fus <- fusSamples[names(fusSamples) %in% derived_tumor_type]
      if (any(derived_tumor_type %notin% names(fusSamples))) {
        missingFusSamples <-
          derived_tumor_type[derived_tumor_type %notin% names(fusSamples)]
        samples_fus_miss <-
          lapply(missingFusSamples , function(x)
            NULL) %>% setNames(missingFusSamples)
        samples_fus <- c(samples_fus , samples_fus_miss)
      }
      fus <- fus[fus$tier %in% c("tier1" , "tier2") ,]
      fus <- fus[fus$frame != "Out-of-frame" ,]
      fusgenes <-
        panel[panel$alteration == "fusion" , "gene_symbol"] %>% unique
      gtOut <-
        fus[fus$Gene_A %in% fusgenes |
              fus$Gene_B %in% fusgenes | fus$FusionPair %in% fusgenes ,]
      object@dataFull$fusions <- list(data = gtOut
                                      , Samples = samples_fus)
    }
  } else {
    object@dataFull$fusions <- list(data = NULL
                                    , Samples = NULL)
  }
  
  #-----------------------------------------------------------------------------
  # LOOK FOR SNV
  #-----------------------------------------------------------------------------
  if (any(cpArguments(object)$panel$alteration == "SNV")) {
    #check if all the tumor_type selected are available
    availability <-
      .checkDataAvailability(tumor_type = tumor_type 
                             , genProfile = "mutations$") %>% lengths
    #in case there is one missing
    if (all(availability == 0)) {
      message(paste0("\n" , noDataMex))
      object@dataFull$mutations <- list(data = NULL
                                        , Samples = NULL)
    } else {
      message("\nRetrieving mutations...")
      mutgenes <-
        panel[panel$alteration == "SNV" , "gene_symbol"] %>% unique
      # There is a bug in cgdsr: when genes like 
      # C10orf12 are queried, they are put in upper case
      # The solution is to revert them to their original case
      mutgenes_translator <-
        data.frame(
          original = mutgenes,
          upper = toupper(mutgenes),
          stringsAsFactors = FALSE
        )
      # cgdsr performs poorly on large queries
      # If the number of requested genes is too large 
      # we subset the SQL query in chunks
      if (length(mutgenes) > 100) {
        #split gene list in blocks of gene_block elements
        myGenesBlocks <- .subsetter(mutgenes , gene_block)
        message(
          paste("Too many genes. We subset data retrieval in"
            ,length(myGenesBlocks)
            ,"blocks of"
            ,gene_block
            ,"genes each"))
        gmOut <-
          lapply(seq_len(length(myGenesBlocks)) , function(x) {
            geneblock <- myGenesBlocks[[x]]
            .getMutations(geneblock , tumor_type = tumor_type , block =
                            x)
          })
        gmOut2 <- bplapply(names(gmOut[[1]]), function(x) {
          box <- lapply(gmOut , function(z)
            z[[x]]$out)
          as.data.frame(rbindlist(box))
        } , BPPARAM = BPPARAM)
        gmOut2 <- as.data.frame(rbindlist(gmOut2))
        if (!is.null(gmOut2)) {
          tumor_type_vec <- unique(gmOut2$genetic_profile_id) %>%
            strsplit(. , "_") %>%
            vapply(., '[' , character(1) , 1)
          names(tumor_type_vec) <-
            unique(gmOut2$genetic_profile_id)
          gmOut2$tumor_type <-
            tumor_type_vec[gmOut2$genetic_profile_id]
          tumor_type_vec2_names <- names(gmOut[[1]])
          tumor_type_vec2 <-
            strsplit(tumor_type_vec2_names , "_") %>% vapply(. , '[' 
                                        , character(1) , 1)
          samples_mut <-
            lapply(unique(tumor_type_vec2) , function(tum) {
              # browser()
              tum <-
                tumor_type_vec2_names[tumor_type_vec2 == tum]
              lapply(gmOut[[1]][tum] , '[[' , "patients") %>% unlist %>% unique
            })
          names(samples_mut) <- unique(tumor_type_vec2)
        }
      } else {
        gmOut <-
          .getMutations(mutgenes , tumor_type = tumor_type , block = NULL)
        gmOut2 <-
          as.data.frame(data.table::rbindlist(lapply(gmOut , '[[' , 1)))
        if (is.data.frame(gmOut2)) {
          if (nrow(gmOut2) == 0) {
            gmOut2 <- NULL
          }
        }
        if (!is.null(gmOut2)) {
          tumor_type_vec <- unique(gmOut2$genetic_profile_id) %>%
            strsplit(. , "_") %>%
            vapply(., '[' , character(1) , 1)
          names(tumor_type_vec) <-
            unique(gmOut2$genetic_profile_id)
          gmOut2$tumor_type <-
            tumor_type_vec[gmOut2$genetic_profile_id]
          tumor_type_vec2_names <- names(gmOut)
          tumor_type_vec2 <-
            strsplit(tumor_type_vec2_names , "_") %>% vapply(. , '[' 
                                              , character(1) , 1)
          samples_mut <-
            lapply(unique(tumor_type_vec2) , function(tum) {
              # browser()
              tum <-
                tumor_type_vec2_names[tumor_type_vec2 == tum]
              lapply(gmOut[tum] , '[[' , "patients") %>% unlist %>% unique
            })
          names(samples_mut) <- unique(tumor_type_vec2)
        }
      }
      if (is.null(gmOut2)) {
        message(
          noDataMex
        )
        object@dataFull$mutations <- list(data = NULL
                                          , Samples = NULL)
      } else {
        # Some clean up
        gmOut2 <- gmOut2[!is.na(gmOut2$mutation_type) ,]
        # Accepted TCGA style variant classification
        # We exclude every mutation that differ from the MAF v2.4 standards
        ktcga_types <- c(
          "In_Frame_Del",
          "In_Frame_Ins",
          "Missense_Mutation",
          "Frame_Shift_Del",
          "Frame_Shift_Ins",
          "Nonsense_Mutation",
          "Splice_Site",
          "Translation_Start_Site",
          "Nonstop_Mutation",
          "Silent",
          "3'UTR",
          "3'Flank",
          "5'UTR",
          "5'Flank",
          "IGR",
          "Intron",
          "RNA",
          "Targeted_Region")
        knotTransc <-
          c(
            "3'UTR",
            "3'Flank",
            "5'UTR",
            "5'Flank",
            "IGR",
            "Intron",
            "RNA",
            "Targeted_Region")
        knonsynonymous <- setdiff(ktcga_types , c("Silent" , knotTransc))
        kmiss_type <- ktcga_types[c(1, 2, 3)]
        ktrunc_type <- ktcga_types[c(4, 5, 6, 7, 8, 9)]
        mutation_type <- mutation_type[1]
        if (mutation_type == "all_nonsynonymous")
          gmOut2 <-
          gmOut2[gmOut2$mutation_type %in% knonsynonymous ,]
        if (mutation_type == "missense")
          gmOut2 <-
          gmOut2[gmOut2$mutation_type %in% kmiss_type ,]
        if (mutation_type == "truncating")
          gmOut2 <-
          gmOut2[gmOut2$mutation_type %in% ktrunc_type ,]
        if (mutation_type == "all_mutations")
          gmOut2 <-
          gmOut2[gmOut2$mutation_type %in% ktcga_types ,]
        # Add new useful columns
        gmOut2$amino_position <-
          ifelse(gmOut2$mutation_type %notin% knotTransc
                 , as.numeric(as.character(
                   str_extract(
                     string = gmOut2$amino_acid_change,
                     pattern = "\\d+"
                   )))
                 , NaN)
        #-----------------------------
        # Chr retrieval in case is NA
        #
        if (any(is.na(gmOut2$chr))) {
          nochr <- gmOut2[is.na(gmOut2$chr) , "gene_symbol"] %>% unique
          ensembl <-
            biomaRt::useMart(host = "grch37.ensembl.org"
                             ,
                             biomart = "ENSEMBL_MART_ENSEMBL"
                             ,
                             dataset = "hsapiens_gene_ensembl")
          Chr_df <- biomaRt::getBM(
            mart = ensembl
            ,
            values = nochr
            ,
            filters = "hgnc_symbol"
            ,
            attributes = c("hgnc_symbol", "chromosome_name")
          )
          Chr_df <-
            Chr_df[Chr_df$chromosome_name %in% c(seq_len(23) 
                              , "X" , "M" , "MT" , "Y") ,] %>% unique
          if (any(nochr %notin% Chr_df$hgnc_symbol)) {
            genesout <- setdiff(nochr , Chr_df$hgnc_symbol)
            stop(paste(
              "Problem retrieving chr name from biomart in" ,
              paste(genesout , collapse = ", ")
            ))
          }
          for (i in nochr) {
            mychr <-
              Chr_df[Chr_df$hgnc_symbol == i , "chromosome_name" 
                     , drop = TRUE] %>% unique %>% unname
            if (length(mychr) != 1) {
              stop(paste(
                i ,
                "is associated with more than 1 chromosome name in biomart"
              ))
            }
            gmOut2[gmOut2$gene_symbol == i , "chr"] <-
              mychr
          }
        }
        #-------------------------------
        # substitute chrMT with chrM
        #
        if (any(gmOut2$chr == "chrMT")) {
          gmOut2[gmOut2$chr == "chrMT" , "chr"] <- "chrM"
        }
        #-----------------------------------------------------------------
        # create genomic position column and get rid of each separate col
        #
        gmOut2$genomic_position <- paste(
          gmOut2$chr
          ,
          gmOut2$start_position
          ,
          paste(gmOut2$reference_allele , gmOut2$variant_allele , sep = ",")
          ,
          sep = ":"
        )
        for (z in c("chr" ,
                    "start_position" ,
                    "reference_allele" ,
                    "variant_allele")) {
          gmOut2[, z] <- NULL
        }
        # Correct TCGA samples name (revert TCGA-11-1234-01 to TCGA-11-1234)
        gmOut2$case_id <- as.character(gmOut2$case_id)
        gmOut2$case_id[grep("^TCGA-" , gmOut2$case_id)] <-
          gmOut2$case_id[grep("^TCGA-" , gmOut2$case_id)] %>%
          strsplit(. , "-") %>%
          vapply(. , function(x)
            paste(x[seq_len(3)] , collapse = "-") , character(1)) %>%
          unlist
        # Correct gene names. If cgdsr put them 
        # upper case, revert them to their original case
        if (any(mutgenes_translator$original != mutgenes_translator$upper)) {
          mutgenes_translator <-
            mutgenes_translator[mutgenes_translator$original != 
                                  mutgenes_translator$upper , , drop=FALSE]
          gmOut2$gene_symbol <-
            .mapvalues(
              gmOut2$gene_symbol
              ,
              from = mutgenes_translator$upper
              ,
              to = mutgenes_translator$original
              ,
              warn_missing = FALSE
            )
        }
        samples_mut <- lapply(samples_mut , function(x) {
          if (is.null(x[1])) {
            return(x)
          }
          x[grep("^TCGA-" , x)] <-
            x[grep("^TCGA-" , x)] %>%
            strsplit(. , "-") %>%
            vapply(. , function(x)
              paste(x[seq_len(3)] , collapse = "-") , character(1)) %>%
            unlist
          return(unique(x))
        })
        names(samples_mut) <- unique(tumor_type_vec2)
        gmOut2 <- unique(gmOut2)
        #Calculate frequencies
        object@dataFull$mutations <- list(data = gmOut2
                                          , Samples = samples_mut)
      }
    }
  } else {
    object@dataFull$mutations <- list(data = NULL
                                      , Samples = NULL)
  }
  
  #-----------------------------------------------------------------------------
  # LOOK FOR CNA
  #-----------------------------------------------------------------------------
  # Let's deal with CNAs. CNA are in theory an easier task
  # The data come in a different form as patients on rows and genes on columns
  # To uniform the output to mutations, we have to melt this matrix
  if (any(cpArguments(object)$panel$alteration == "CNA")) {
    availability <-
      .checkDataAvailability(tumor_type = tumor_type 
                             , genProfile = "gistic|cna$") %>% lengths
    if (all(availability == 0)) {
      message(
        paste("\nNo copynumber alteration data" 
              ,"available for the specified tumor types...")
      )
      object@dataFull$copynumber <- list(data = NULL
                                         , Samples = NULL)
    } else {
      message("\nRetrieving copy number alterations...")
      cnagenes <-
        panel[panel$alteration == "CNA" , "gene_symbol"] %>% unique
      # There is a bug in cgdsr: when genes like 
      # C10orf12 are queried, they are put in upper case
      # The solution is to revert them to their original case
      cnagenes_translator <-
        data.frame(
          original = cnagenes,
          upper = toupper(cnagenes),
          stringsAsFactors = FALSE
        )
      # cgdsr performs poorly on large queries
      # If the number of requested genes is too 
      # large we subset the SQL query in chunks
      if (length(cnagenes) > 100) {
        myGenesBlocks_cna <- .subsetter(cnagenes , gene_block)
        message(
          paste(
            "Too many genes. We subset data retrieval in"
            ,
            length(myGenesBlocks_cna)
            ,
            "blocks of"
            ,
            gene_block
            ,
            "genes each"
          )
        )
        gcOut <-
          lapply(seq_len(length(myGenesBlocks_cna)) , function(x) {
            geneblock <- myGenesBlocks_cna[[x]]
            tobereturned <-
              .getCNA(geneblock , tumor_type = tumor_type , block = x)
            return(tobereturned)
          })
        gcOut_melt <-
          lapply(gcOut , function(tobereturned) {
            for (x in names(tobereturned)) {
              if (!is.null(tobereturned[[x]]$out)) {
                tobereturned[[x]]$out$case_id <-
                  gsub("\\." , "-" , rownames(tobereturned[[x]]$out))
                tobereturned[[x]]$out$genetic_profile_id <-
                  x
                tum_type <-
                  strsplit(x , "_")[[1]][1]
                tobereturned[[x]]$out$tumor_type <-
                  tum_type
              }
            }
            tobereturned <-
              as.data.frame(rbindlist(lapply(tobereturned , '[[' , 1)))
            tobereturned_melt <-
              reshape2::melt(
                tobereturned ,
                id.vars = c("case_id" , "genetic_profile_id" , "tumor_type")
                ,
                value.name = "CNA" ,
                variable.name = "gene_symbol"
                ,
                factorsAsStrings = FALSE
              )
            return(tobereturned_melt)
          })
        gcOut_melt <-
          as.data.frame(rbindlist(gcOut_melt)) %>% .changeFactor(.)
        tumor_type_vec <-
          unique(gcOut_melt$genetic_profile_id) %>%
          strsplit(. , "_") %>%
          vapply(., '[' , character(1) , 1)
        names(tumor_type_vec) <-
          unique(gcOut_melt$genetic_profile_id)
        gcOut_melt$tumor_type <-
          tumor_type_vec[gcOut_melt$genetic_profile_id]
        tumor_type_vec2_names <- names(gcOut[[1]])
        tumor_type_vec2 <-
          strsplit(tumor_type_vec2_names , "_") %>% 
          vapply(. , '[' , character(1) , 1)
        samples_cna <-
          lapply(unique(tumor_type_vec2) , function(tum) {
            # browser()
            tum <-
              tumor_type_vec2_names[tumor_type_vec2 == tum]
            lapply(gcOut[[1]][tum] , '[[' , "patients") %>% unlist %>% unique
          })
        names(samples_cna) <- unique(tumor_type_vec2)
      } else {
        gcOut <- .getCNA(cnagenes , tumor_type = tumor_type , block = NULL)
        for (x in names(gcOut)) {
          if (!is.null(gcOut[[x]]$out)) {
            gcOut[[x]]$out$case_id <-
              gsub("\\." , "-" , rownames(gcOut[[x]]$out))
            gcOut[[x]]$out$genetic_profile_id <- x
            tum_type <- strsplit(x , "_")[[1]][1]
            gcOut[[x]]$out$tumor_type <- tum_type
          }
        }
        gcOut2 <-
          as.data.frame(rbindlist(lapply(gcOut , '[[' , 1)))
        gcOut_melt <-
          reshape2::melt(
            gcOut2 ,
            id.vars = c("case_id" , "genetic_profile_id" , "tumor_type")
            ,
            value.name = "CNA" ,
            variable.name = "gene_symbol"
            ,
            factorsAsStrings = FALSE
          ) %>% .changeFactor(.)
        
        tumor_type_vec <-
          unique(gcOut2$genetic_profile_id) %>%
          strsplit(. , "_") %>%
          vapply(., '[' , character(1) , 1)
        names(tumor_type_vec) <-
          unique(gcOut2$genetic_profile_id)
        gcOut_melt$tumor_type <-
          tumor_type_vec[gcOut2$genetic_profile_id]
        tumor_type_vec2_names <- names(gcOut)
        tumor_type_vec2 <-
          strsplit(tumor_type_vec2_names , "_") %>% 
          vapply(., '[' , character(1) , 1)
        samples_cna <-
          lapply(unique(tumor_type_vec2) , function(tum) {
            # browser()
            tum <-
              tumor_type_vec2_names[tumor_type_vec2 == tum]
            lapply(gcOut[tum] , '[[' , "patients") %>% unlist %>% unique
          })
        names(samples_cna) <- unique(tumor_type_vec2)
      }
      # Correct TCGA samples name (revert TCGA-11-1234-01 to TCGA-11-1234)
      gcOut_melt$case_id <- as.character(gcOut_melt$case_id)
      gcOut_melt$case_id[grep("^TCGA-" , gcOut_melt$case_id)] <-
        gcOut_melt$case_id[grep("^TCGA-" , gcOut_melt$case_id)] %>%
        strsplit(. , "-") %>%
        vapply(. , function(x)
          paste(x[seq_len(3)] , collapse = "-") , character(1)) %>%
        unlist
      samples_cna <- lapply(samples_cna , function(x) {
        if (is.null(x[1]))
          return(x)
        x[grep("^TCGA-" , x)] <-
          x[grep("^TCGA-" , x)] %>%
          strsplit(. , "-") %>%
          vapply(. , function(x)
            paste(x[seq_len(3)] , collapse = "-") , character(1)) %>%
          unlist
        return(unique(x))
      })
      # Correct gene names. If cgdsr put them upper case, 
      # revert them to their original case
      if (any(cnagenes_translator$original != cnagenes_translator$upper)) {
        cnagenes_translator <-
          cnagenes_translator[cnagenes_translator$original != 
                                cnagenes_translator$upper , , drop =FALSE]
        gcOut_melt$gene_symbol <-
          .mapvalues(
            gcOut_melt$gene_symbol
            ,
            from = cnagenes_translator$upper
            ,
            to = cnagenes_translator$original
            ,
            warn_missing = FALSE
          )
      }
      # All this mess is because gistic applies 
      # different scores at different datasets
      # This create confusion, because if the 
      # same patient belongs to different study
      # Its gistic copynumber can be different
      gcOut_melt2 <-
        data.frame(
          tumor_type = factor(gcOut_melt$tumor_type 
                              , levels = derived_tumor_type)
          ,
          gene_symbol = factor(gcOut_melt$gene_symbol , levels = cnagenes)
          ,
          CNA = ifelse(is.na(gcOut_melt$CNA) , 0 , gcOut_melt$CNA)
          ,
          case_id = gcOut_melt$case_id
          ,
          stringsAsFactors = FALSE
        ) %>% unique
      gcOut_melt2 <-
        aggregate(CNA ~ tumor_type + gene_symbol + case_id , gcOut_melt2 , mean)
      gcOut_melt2$CNA <-
        factor(as.character(round(gcOut_melt2$CNA)) 
               , levels = c("2" , "1" , "0" , "-1" , "-2"))
      # NOTE: the format of Freq is different 
      # from mutations. It is a list (not a matrix)
      # It contains an element for each tumor type
      # Each tumor type is a dataframe with amplification, 
      # normal and deletion freqeuncy
      gcOut_melt3 <-
        data.frame(
          gene_symbol = as.character(gcOut_melt2$gene_symbol)
          ,
          CNA = ifelse(
            as.character(gcOut_melt2$CNA) %in% c("1" , "0" , "-1") ,
            "normal"
            ,
            ifelse(
              as.character(gcOut_melt2$CNA) == "2" ,
              "amplification" ,
              "deletion"
            )
          )
          ,
          case_id = gcOut_melt2$case_id
          ,
          tumor_type = as.character(gcOut_melt2$tumor_type)
          ,
          CNAvalue = as.character(gcOut_melt2$CNA)
          ,
          stringsAsFactors = FALSE
        ) %>% unique
      object@dataFull$copynumber <- list(data = gcOut_melt3
                                         , Samples = samples_cna)
    }
  } else {
    object@dataFull$copynumber <- list(data = NULL
                                       , Samples = NULL)
  }
  
  
  #-----------------------------------------------------------------------------
  # LOOK FOR EXPRESSION
  #-----------------------------------------------------------------------------
  # Let's deal with Expression. Expression comes exactly as CNA
  # To uniform the output to mutations, we have to melt this matrix
  if (any(cpArguments(object)$panel$alteration == "expression")) {
    availability <-
      .checkDataAvailability(tumor_type = tumor_type 
                , genProfile = "_rna_seq_v2_mrna_median_Zscores$") %>% lengths
    if (all(availability == 0)) {
      message(paste("\nSorry, no expression data available" 
                    ,"for the specified tumor types..."))
      object@dataFull$expression <- list(data = NULL
                                         , Samples = NULL)
    } else {
      message("\nRetrieving Expression data...")
      exprgenes <-
        panel[panel$alteration == "expression" , "gene_symbol"] %>% unique
      # There is a bug in cgdsr: when genes like 
      # C10orf12 are queried, they are put in upper case
      # The solution is to revert them to their original case
      exprgenes_translator <-
        data.frame(
          original = exprgenes,
          upper = toupper(exprgenes),
          stringsAsFactors = FALSE
        )
      # cgdsr performs poorly on large queries
      # If the number of requested genes is too large 
      # we subset the SQL query in chunks
      if (length(exprgenes) > 100) {
        myGenesBlocks_expr <- .subsetter(exprgenes , gene_block)
        message(
          paste(
            "Too many genes. We subset data retrieval in"
            ,
            length(myGenesBlocks_expr)
            ,
            "blocks of"
            ,
            gene_block
            ,
            "genes each"
          )
        )
        geOut <-
          lapply(seq_len(length(myGenesBlocks_expr)) , function(x) {
            geneblock <- myGenesBlocks_expr[[x]]
            tobereturned <-
              .getExpression(geneblock , tumor_type = tumor_type , block = x)
            return(tobereturned)
          })
        geOut_melt <-
          lapply(geOut , function(tobereturned) {
            for (x in names(tobereturned)) {
              if (!is.null(tobereturned[[x]]$out)) {
                tobereturned[[x]]$out$case_id <-
                  gsub("\\." , "-" , rownames(tobereturned[[x]]$out))
                tobereturned[[x]]$out$genetic_profile_id <-
                  x
                tum_type <-
                  strsplit(x , "_")[[1]][1]
                tobereturned[[x]]$out$tumor_type <-
                  tum_type
              }
            }
            tobereturned <-
              as.data.frame(rbindlist(lapply(tobereturned , '[[' , 1)))
            tobereturned_melt <-
              reshape2::melt(
                tobereturned ,
                id.vars = c("case_id" , "genetic_profile_id" , "tumor_type")
                ,
                value.name = "expression" ,
                variable.name = "gene_symbol"
                ,
                factorsAsStrings = FALSE
              )
            return(tobereturned_melt)
          })
        geOut_melt <-
          as.data.frame(rbindlist(geOut_melt)) %>% .changeFactor(.)
        tumor_type_vec <-
          unique(geOut_melt$genetic_profile_id) %>%
          strsplit(. , "_") %>%
          vapply(., '[' , character(1) , 1)
        names(tumor_type_vec) <-
          unique(geOut_melt$genetic_profile_id)
        geOut_melt$tumor_type <-
          tumor_type_vec[geOut_melt$genetic_profile_id]
        tumor_type_vec2_names <- names(geOut[[1]])
        tumor_type_vec2 <-
          strsplit(tumor_type_vec2_names , "_") %>% 
          vapply(., '[' , character(1) , 1)
        samples_expr <-
          lapply(unique(tumor_type_vec2) , function(tum) {
            tum <- tumor_type_vec2_names[tumor_type_vec2 == tum]
            lapply(geOut[[1]][tum] , '[[' , "patients") %>% unlist %>% unique
          })
        names(samples_expr) <- unique(tumor_type_vec2)
      } else {
        geOut <-
          .getExpression(exprgenes , tumor_type = tumor_type , block = NULL)
        for (x in names(geOut)) {
          if (!is.null(geOut[[x]]$out)) {
            geOut[[x]]$out$case_id <-
              gsub("\\." , "-" , rownames(geOut[[x]]$out))
            geOut[[x]]$out$genetic_profile_id <- x
            tum_type <- strsplit(x , "_")[[1]][1]
            geOut[[x]]$out$tumor_type <- tum_type
          }
        }
        geOut2 <-
          as.data.frame(rbindlist(lapply(geOut , '[[' , 1)))
        geOut_melt <-
          reshape2::melt(
            geOut2 ,
            id.vars = c("case_id" , "genetic_profile_id" , "tumor_type")
            ,
            value.name = "expression" ,
            variable.name = "gene_symbol"
            ,
            factorsAsStrings = FALSE
          ) %>% .changeFactor(.)
        tumor_type_vec <-
          unique(geOut2$genetic_profile_id) %>%
          strsplit(. , "_") %>%
          vapply(., '[' , character(1) , 1)
        names(tumor_type_vec) <-
          unique(geOut2$genetic_profile_id)
        geOut2$tumor_type <-
          tumor_type_vec[geOut2$genetic_profile_id]
        tumor_type_vec2_names <- names(geOut)
        tumor_type_vec2 <-
          strsplit(tumor_type_vec2_names , "_") %>% 
          vapply(., '[' , character(1) , 1)
        samples_expr <-
          lapply(unique(tumor_type_vec2) , function(tum) {
            tum <- tumor_type_vec2_names[tumor_type_vec2 == tum]
            lapply(geOut[tum] , '[[' , "patients") %>% unlist %>% unique
          })
        names(samples_expr) <- unique(tumor_type_vec2)
      }
      # Correct TCGA samples name (revert TCGA-11-1234-01 to TCGA-11-1234)
      geOut_melt$case_id <- as.character(geOut_melt$case_id)
      geOut_melt$case_id[grep("^TCGA-" , geOut_melt$case_id)] <-
        geOut_melt$case_id[grep("^TCGA-" , geOut_melt$case_id)] %>%
        strsplit(. , "-") %>%
        vapply(. , function(x)
          paste(x[seq_len(3)] , collapse = "-") , character(1)) %>%
        unlist
      samples_expr <- lapply(samples_expr , function(x) {
        if (is.null(x[1]))
          return(x)
        x[grep("^TCGA-" , x)] <-
          x[grep("^TCGA-" , x)] %>%
          strsplit(. , "-") %>%
          vapply(. , function(x)
            paste(x[seq_len(3)] , collapse = "-") , character(1)) %>%
          unlist
        return(unique(x))
      })
      # Correct gene names. If cgdsr put them upper case, 
      # revert them to their original case
      if (any(exprgenes_translator$original != exprgenes_translator$upper)) {
        exprgenes_translator <-
          exprgenes_translator[exprgenes_translator$original != 
                                 exprgenes_translator$upper , , drop =FALSE]
        geOut_melt$gene_symbol <-
          .mapvalues(
            geOut_melt$gene_symbol
            ,
            from = exprgenes_translator$upper
            ,
            to = exprgenes_translator$original
            ,
            warn_missing = FALSE
          )
      }
      geOut_melt2 <-
        unique(geOut_melt[, c("tumor_type" , "case_id" 
                              , "gene_symbol" , "expression")])
      geOut_melt2$expression <-
        .noNA(geOut_melt2$expression , subs = 0)
      geOut_melt2 <-
        aggregate(expression ~ tumor_type + gene_symbol + case_id ,
                  geOut_melt2 ,
                  mean)
      geOut_melt2$expressionValue <- geOut_melt2$expression
      expression_dichotomize <-
        ifelse(
          geOut_melt2$expression >= expr_z_score ,
          "up"
          ,
          ifelse(
            geOut_melt2$expression <= (-expr_z_score) ,
            "down"
            ,
            "normal"
          )
        )
      geOut_melt2$expression <-
        factor(expression_dichotomize , levels = c("up" , "normal" , "down"))
      # NOTE: the format of Freq is different 
      # from mutations. It is a list (not a matrix)
      # It contains a element for each tumor type
      # Each tumor type is a dataframe with amplification, 
      # normal and deletion frequency
      geOut_melt3 <-
        data.frame(
          gene_symbol = as.character(geOut_melt2$gene_symbol)
          ,
          expression = as.character(geOut_melt2$expression)
          ,
          case_id = geOut_melt2$case_id
          ,
          tumor_type = as.character(geOut_melt2$tumor_type)
          ,
          expressionValue = geOut_melt2$expressionValue
          ,
          stringsAsFactors = FALSE
        )
      object@dataFull$expression <- list(data = geOut_melt3
                                         , Samples = samples_expr)
    }
  } else {
    object@dataFull$expression <- list(data = NULL
                                       , Samples = NULL)
  }
  return(object)
})