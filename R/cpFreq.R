setGeneric('cpFreq', function(object
						, alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
						, mutations.specs=c(NA ,"mutation_type","amino_acid_change","amino_position","genomic_position")
						, fusions.specs=c("bygene" , "byfusionpair")
						, tumor_type=NULL
						, tumor.weights=NULL
						, tumor.freqs=NULL
						, collapseMutationByGene=TRUE  
						, freq=c("relative" , "absolute")
						)  {
	standardGeneric('cpFreq')
	})
setMethod('cpFreq', 'CancerPanel', function(object
						, alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
						, mutations.specs=c(NA ,"mutation_type","amino_acid_change","amino_position","genomic_position")
						, fusions.specs=c("bygene" , "byfusionpair")
						, tumor_type=NULL
						, tumor.weights=NULL
						, tumor.freqs=NULL
						, collapseMutationByGene=TRUE  
						, freq=c("relative" , "absolute")
						)
{
	# Check input parameters
	if(is.null(object)){
		stop("No CancerPanel object provided")
	}
	if(is.null(object@dataFull)){
		stop("no getAlterations has been called")	
	}
	if(is.null(alterationType)){
		stop("select one alterationType")
	} else {
		if(length(alterationType)>1){
			alterationType <- alterationType[1]
			warning("More than one alterationType selected. The first one was chosen:" %++% alterationType)
		}
		if(alterationType %notin% c("copynumber" , "expression" , "mutations" , "fusions")){
			stop('AlterationType can only be one among "copynumber" , "expression" , "mutations" , "fusions"')
		}
	}
	# Check tumor_type
	if(!is.null(tumor_type)){
		if(!all(tumor_type %in% object@arguments$tumor_type)){
			stop("You selected a tumor_type that has no data in this CancerPanel object")
		}
	}
	if(is.null(mutations.specs)){
		mutations.specs <- NA
	} else {
		mutations.specs <- mutations.specs[1]
		if(mutations.specs %notin% c(NA ,"mutation_type","amino_acid_change","amino_position","genomic_position")){
			stop("Invalid mutations.specs parameter. mutations.specs can only be NA 'mutation_type' 'amino_acid_change' 'amino_position' 'genomic_position'")
		}
	}
	if(is.null(fusions.specs)){
		fusions.specs <- "bygene"
	} else {
		fusions.specs <- fusions.specs[1]
		if(fusions.specs %notin% c("bygene" , "byfusionpair")){
			stop("Invalid fusions.specs parameter. fusions.specs can be 'bygene' , 'byfusionpair'")
		}
	}
	if(is.null(freq)){
		freq <- "absolute"
	} else {
		freq <- freq[1]
		if(freq %notin% c("absolute" , "relative")){
			stop("freq parameter can only be absolute or relative")
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
		if(freq=="absolute"){
			stop("If you set tumor.freqs, you cannot obtain absolute freq. set freq = 'relative' ")
		}
	}
	# Define the full set of gene from the panel, even if there are no alterations
	kConversionTab <- c(mutations="SNV" , copynumber="CNA" , fusions="fusion" , expression="expression")
	panel <- object@arguments$panel
	genesToCheck <- panel[ panel$alteration==kConversionTab[alterationType] , "gene_symbol"] %>% unique
	if(length(genesToCheck)==0){
		stop("No genes are in the panel for" %++% alterationType)
	}
	# Retrieve the data from the object
	mydata <- object@dataFull[[alterationType]]$data
	mysamples <- object@dataFull[[alterationType]][["Samples"]]
	mysamples <- lapply(mysamples , unique)
	if(is.null(tumor_type)){
		tumor_type <- names(mysamples)
	} else {
		if(all(tumor_type %in% names(mysamples))){
			mysamples <- mysamples[tumor_type]	
		} else {
			tumor_type_not_present <- setdiff(tumor_type ,names(mysamples))
			stop("The following tumor_type are not present" %++% paste(tumor_type_not_present , collapse=", "))
		}
		mydata <- mydata[ mydata$tumor_type %in% tumor_type , ]
		if(nrow(mydata)==0){
			stop("No alterations for the tumor_type and alterationType requested")
		}
	}
	mysamples$all_tumors <- unlist(unname(mysamples))
	# Reduce the dataset to a specific amount of patients sampled at random
	# among the different tumor types
	if(!is.null(tumor.weights)){
		newDataAndSamps <- .tumor.weights.machine(tumor.weights , mysamples , mydata )
		mydata <- newDataAndSamps[[1]]
		mysamples <- newDataAndSamps[[2]]
	}
	# Reduce the dataset to only the patients included in mysamples
	mydata <- mydata[ mydata$case_id %in% mysamples$all_tumors , ]
	# mydata$alteration_id <- strsplit(mydata$alteration_id , "_") %>% sapply(. , '[' , 1)
	# Reduce dataset according to tumor types selected
	if(!is.null(tumor_type)){
		mydata <- mydata[ mydata$tumor_type %in% tumor_type , ]
	}
	# If tumor.freqs is set, we basically run this method in recursive mode for each tumor type
	# and then aggregate everything
	if(!is.null(tumor.freqs)){
		FreqRecurse <- lapply(names(tumor.freqs) , function(tum){
				out <- tryCatch( cpFreq(object
						, alterationType=alterationType
						, mutations.specs=mutations.specs
						, fusions.specs=fusions.specs
						, tumor_type=tum
						, tumor.weights=NULL
						, tumor.freqs=NULL
						, collapseMutationByGene=collapseMutationByGene
						, freq="relative"
						) , error = function(e) return(NULL))
				if(is.null(out)){
					return(NULL)
				}
				if(alterationType=="mutations"){
					freqColsInternal <- colnames(out) %notin% c("gene_symbol" , mutations.specs)
					colnames(out)[freqColsInternal] <- paste0(colnames(out)[freqColsInternal] , tum)
				}
				return(out)
			})
		names(FreqRecurse) <- names(tumor.freqs)
		notNull <- !sapply(FreqRecurse , is.null)
		FreqRecurse <- FreqRecurse[notNull]
		tumor.freqs <- tumor.freqs[notNull]
		if(alterationType=="mutations"){
			if(is.na(mutations.specs)){
				FreqRecurse_merge <- suppressWarnings(.mergeThemAll(FreqRecurse , by="gene_symbol" , all=TRUE))
			} else {
				FreqRecurse_merge <- suppressWarnings(.mergeThemAll(FreqRecurse , by=c("gene_symbol" , mutations.specs) , all=TRUE))
			}
			freqCols <- grep("^freq" , colnames(FreqRecurse_merge) , value=TRUE)
			for(i in freqCols){
				FreqRecurse_merge[is.na(FreqRecurse_merge[ , i]) , i] <- 0
			}
			FreqRecurse_weight <- cbind( FreqRecurse_merge[ , colnames(FreqRecurse_merge) %notin% freqCols , drop=FALSE]
										, apply(FreqRecurse_merge[ , freqCols , drop=FALSE] , 1 , function(x) sum(x*tumor.freqs))
										)
			if(is.na(mutations.specs)){
				colnames(FreqRecurse_weight) <- c("gene_symbol" , "freq")
			} else {
				colnames(FreqRecurse_weight) <- c("gene_symbol" , mutations.specs , "freq")
			}
		}
		if(alterationType=="copynumber" | alterationType=="expression"){
			FreqRecurse <- lapply(names(tumor.freqs) , function(tum){
				out <- FreqRecurse[[tum]]*tumor.freqs[tum]
				out$gene_symbol <- rownames(out)
				return(out)
			})
			FreqRecurse <- as.data.frame(rbindlist(FreqRecurse) , stringsAsFactors=FALSE)
			if(alterationType=="copynumber"){
				FreqRecurse_weight <- aggregate( cbind(amplification   ,  deletion     ,normal)~gene_symbol , FreqRecurse , FUN=sum)
			} else {
				FreqRecurse_weight <- aggregate( cbind(up   ,  down     ,normal)~gene_symbol , FreqRecurse , FUN=sum)
			}
			rownames(FreqRecurse_weight) <- FreqRecurse_weight$gene_symbol
			FreqRecurse_weight$gene_symbol <- NULL
		}
		if(alterationType=="fusions"){
			FreqRecurse_merge <- .mergeThemAll(FreqRecurse , by="gene_symbol" , all=TRUE)
			freqCols <- grep("^freq." , colnames(FreqRecurse_merge) , value=TRUE)
			for(i in freqCols){
				FreqRecurse_merge[is.na(FreqRecurse_merge[ , i]) , i] <- 0
			}
			FreqRecurse_weight <- cbind( FreqRecurse_merge[ , colnames(FreqRecurse_merge) %notin% freqCols , drop=FALSE]
										, apply(FreqRecurse_merge[ , freqCols , drop=FALSE] , 1 , function(x) sum(x*tumor.freqs))
										)
			colnames(FreqRecurse_weight) <- c("gene_symbol" , "freq")
		}
		return(FreqRecurse_weight)
	}
	# According to alterationType, the way the frequencies are calculated may differ
	if(alterationType=="mutations"){
		if(is.na(mutations.specs)){
			if(collapseMutationByGene) {
				mydata <- unique(mydata[ , c("gene_symbol" , "case_id")])
			}
			agg <- aggregate(case_id~gene_symbol , mydata 
				, FUN=length , simplify=TRUE)
			colnames(agg) <- c("gene_symbol" , "freq")
			missing_genes <- setdiff(genesToCheck , unique(mydata$gene_symbol))
			if(length(missing_genes)>0){
				agg <- rbind(agg , data.frame(gene_symbol=missing_genes , freq=0))
			}
		} else {
			compose <- paste("case_id ~ gene_symbol +" , mutations.specs)
			agg <- aggregate( as.formula(compose) , unique(mydata[ , c("gene_symbol" , "case_id" , mutations.specs)])
				, FUN=length , simplify=TRUE)
			colnames(agg) <- c("gene_symbol" , mutations.specs , "freq")
			missing_genes <- setdiff(levels(agg$gene_symbol) , unique(mydata$gene_symbol))
			if(length(missing_genes)>0){
				toBeRbind <- data.frame(missing_genes 
										,NA
										,0)
				colnames(toBeRbind) <- c("gene_symbol" , mutations.specs , "freq")
				agg <- rbind(agg 
							, toBeRbind)
			}
		}
		if(freq=="relative"){
			agg$freq <- agg$freq/length(mysamples$all_tumors)
		}
		return(agg)
	}

	if(alterationType=="copynumber"){
		mydata$gene_symbol <- factor(mydata$gene_symbol , levels=genesToCheck)
		out_table <- table(mydata$gene_symbol , factor(mydata$CNA , levels=c("amplification" , "deletion" , "normal"))) %>% 
					as.data.frame.matrix
		if(freq=="relative"){
			out_table <- out_table/length(mysamples$all_tumors)
		}
		return(out_table)
	}

	if(alterationType=="expression"){
		mydata$gene_symbol <- factor(mydata$gene_symbol , levels=genesToCheck)
		out_table <- table(mydata$gene_symbol , factor(mydata$expression , levels=c("up" , "down" , "normal"))) %>%
					as.data.frame.matrix
		if(freq=="relative"){
			out_table <- out_table/length(mysamples$all_tumors)
		}
		return(out_table)
	}

	if(alterationType=="fusions"){
		if(fusions.specs=="bygene"){
			agg1 <- aggregate(case_id~Gene_A , mydata , FUN=length , simplify=TRUE)
			agg2 <- aggregate(case_id~Gene_B , mydata , FUN=length , simplify=TRUE)
			colnames(agg1) <- colnames(agg2) <- c("gene_symbol" , "freq")
			agg <- rbind(agg1 , agg2)
			agg <- aggregate(freq~gene_symbol , agg , FUN=sum)
			fusiongenes <- strsplit(genesToCheck,"__") %>% unlist %>% unique
			missing_genes <- setdiff(fusiongenes  , unique(agg$gene_symbol))
			if(length(missing_genes)>0){
				agg <- rbind(agg , data.frame(gene_symbol=missing_genes , freq=0))
			}
		}
		if(fusions.specs=="byfusionpair"){
			agg <- aggregate(case_id ~ FusionPair , mydata , FUN=length , simplify=TRUE)
			colnames(agg) <- c("gene_symbol" , "freq")	
		}
		if(freq=="relative"){
			agg$freq <- agg$freq/length(mysamples$all_tumors)
		}
		return(agg)
	}
})