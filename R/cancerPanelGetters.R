###############################################################################
# CLASS METHODS DECLARATION
###############################################################################
## Method cpArguments
setGeneric('cpArguments', function(object) standardGeneric('cpArguments'))
setMethod('cpArguments', 'CancerPanel', function(object) {
    return(object@arguments)
    })

# -----------------------------------------------------------------------------
## Method cpData
setGeneric('cpData', function(object) standardGeneric('cpData'))
setMethod('cpData', 'CancerPanel', function(object) {
    return(object@dataFull)
    })

# -----------------------------------------------------------------------------
## Method cpDataSubset
setGeneric('cpDataSubset', function(object) standardGeneric('cpDataSubset'))
setMethod('cpDataSubset', 'CancerPanel', function(object) {
    return(object@dataSubset)
    })