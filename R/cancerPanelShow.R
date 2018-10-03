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

    if(!identical(object@dataFull , list())){
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
        alterationType <- c("copynumber" , "expression" , "mutations" , "fusions")
        allcombs <- lapply( seq.int(1 , length(alterationType) , 1) , function(x) {
          combn(alterationType , x , simplify=FALSE)
        }) %>% unlist(recursive = FALSE)
        mytums <- object@arguments$tumor_type
        sampSummary <-   lapply( allcombs , function(comb) {
        allsamps <- lapply( object@dataFull[comb] , '[[' , 'Samples')
          vapply( mytums , function(tum) {
            length(Reduce("intersect" , lapply(allsamps , '[[' , tum)))
          } , numeric(1))
        }) %>% do.call("rbind" , .)
        sampSummary <- cbind( Combinations = vapply( allcombs , function(x) paste(x , collapse=",") , character(1))
                            , sampSummary)
        message("\nThe number of samples for each combination of alteration types is as follow:")
        message(paste0(capture.output(sampSummary), collapse = "\n"))
    } else {
      message("The object contains no data")
    }
})