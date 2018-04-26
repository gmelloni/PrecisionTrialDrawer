setMethod('show', 'CancerPanel', function(object) {
	nItems <- length(object@arguments$genedata$gene_symbol)
	message(paste('\nCancerPanel object with', nItems, 'genes:'))
	items <- object@arguments$genedata$gene_symbol
	itemsText <- paste(head(items), collapse=', ')
	message(paste(itemsText, ifelse(nItems>6, ', ...\n', '\n'), sep=''))

	nItems <- length(object@arguments$drugs %>% .[.!=""])
	message(paste('and', nItems, 'drugs:'))
	if(nItems!=0){
		items <- object@arguments$drugs %>% .[.!=""]
		itemsText <- paste(head(items), collapse=', ')
		message(paste(itemsText, ifelse(nItems>6, ', ...\n', '\n'), sep=''))
	}

	message(paste("The panel contains alterations of the following types:" 
			, paste(unique(object@arguments$panel$alteration) , collapse=", ")))

	if(!is.null(object@dataFull)){
		message("\n")
		for(i in c("mutations" , "copynumber" , "fusions" , "expression")){
			if(!is.null(object@dataFull[[i]]$data)){
				message("The object contains" %++%
						i %++%
						"data for the tumor types:" %++%
						paste(sort(unique(object@dataFull[[i]]$data$tumor_type)) , collapse=", ")
						)
			} else {
				message("No" %++% i %++% "data")
			}
		}
	}
})