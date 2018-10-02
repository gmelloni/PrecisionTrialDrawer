setGeneric('filterFusions', function(object , filtered , mode = c("exclude" , "keep")) {
    standardGeneric('filterFusions')
    })
setMethod('filterFusions', 'CancerPanel', function(object , filtered , mode = c("exclude" , "keep"))
{
    # browser()
    if(!is.character(filtered)){
        stop("filtered must be a character vector")
    }
    mode <- mode[1]
    if(mode %notin% c("exclude" , "keep")){
        stop("mode can only be exclude or keep")
    }
    if(identical(object@dataFull , list())){
        stop("dataFull slot is empty, no data to filter. run getAlterations")
    }
    if(is.null(object@dataFull$fusions$data)){
        stop("No fusion data to filter")
    }
    #---------------------------------
    # Check fusion format
    #
    if(any(!grepl("__" , filtered))){
        stop("Fusion format not recognized, it should be like gene1__gene2 ", paste(filtered[!grepl("__" , filtered)], collapse=", "))
    }
    #-----------------------------------
    # NOW LET'S FILTER
    #
    fus <- object@dataFull$fusions$data
    if(mode == "exclude"){
        myNewRows <- fus$FusionPair %notin% filtered
    } else {
        myNewRows <- fus$FusionPair %in% filtered
    }
    fus <- fus[ myNewRows , , drop=FALSE]
    object@dataFull$fusions$data <- fus
    object <- subsetAlterations(object)
    return(object)
})