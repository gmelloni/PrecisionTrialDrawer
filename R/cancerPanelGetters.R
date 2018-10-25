###############################################################################
# CLASS METHODS DECLARATION
###############################################################################
## Method cpArguments
setGeneric('cpArguments', function(object) standardGeneric('cpArguments'))
setMethod('cpArguments', 'CancerPanel', function(object) {
    object@arguments
})
# setGeneric('cpArguments<-', function(object, value) {
# standardGeneric('cpArguments<-')})
# setReplaceMethod('cpArguments', 'CancerPanel', function(object, value) {
#   object@arguments <- value
#   validObject(object)
#   object
# })

# -----------------------------------------------------------------------------
## Method cpData
setGeneric('cpData', function(object) standardGeneric('cpData'))
setMethod('cpData', 'CancerPanel', function(object) {
    object@dataFull
})
# setGeneric('cpData<-', function(object, value) standardGeneric('cpData<-'))
# setReplaceMethod('cpData', 'CancerPanel', function(object, value) {
#   object@arguments <- value
#   validObject(object)
#   object
# })

# -----------------------------------------------------------------------------
## Method cpDataSubset
setGeneric('cpDataSubset', function(object) standardGeneric('cpDataSubset'))
setMethod('cpDataSubset', 'CancerPanel', function(object) {
    object@dataSubset
})
# setGeneric('cpDataSubset<-', function(object, value) {
# standardGeneric('cpDataSubset<-')})
# setReplaceMethod('cpDataSubset', 'CancerPanel', function(object, value) {
#   object@dataSubset <- value
#   validObject(object)
#   object
# })