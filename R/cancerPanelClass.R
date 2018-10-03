################################################################################
# CLASS DECLARATION
################################################################################
setClass('CancerPanel', 
    slots=c(
        arguments='list'
        ,dataFull='list'
        ,dataSubset='list'
        )
    ,validity=function(object) {
    if(is.list(object@arguments) || 
       is.list(object@dataFull) ||
       is.list(object@dataSubset)) {
      return("CancerPanel slots are lists")
      }
      return(FALSE)
        }
)

################################################################################
# CLASS CONSTRUCTOR
################################################################################
# Validate input
# Fetch transcript size
# Adjust it by alteration specification
# Add info to the object
newCancerPanel <- function(panel , rules=NULL 
                           , padding_length=100 , utr=FALSE 
                           , canonicalTranscript=TRUE , myhost="www.ensembl.org"
                           )
{
    message("Checking panel construction...")
  
  ##################################################
  # CHECK ARGUMENTS
  # ------------------------------------------------
  if( is.null(panel) )
      stop('You should enter a dataframe containing your panel')
  panel <- .panelCheck(panel)
  
  # rules requires a more delicate procedure for the check
  # If we retrieve it from the object it can still be NULL
  if(!is.null(rules)){
    # If we are checking a rules panel, we split the checks in two functions
    # druggability is scorporated and contain only the cases in which an 
    # entire drug is excluded/included from certain tumor types
    druggabilityWhich <- which( apply(panel[ , c("gene_symbol" , "alteration" 
                                , "exact_alteration", "mutation_specification") 
                                , drop=FALSE] 
                                , 1 , function(x) all( x == "")))
    if(length(druggabilityWhich)>0){
      # If there are druggability rules, perform check
      druggability <- rules[ druggabilityWhich 
      , c("drug" , "group" , "tumor_type" , "in_out"), drop=FALSE]
      druggability_full <- .druggabilityCheck(druggability 
        , tumor_type = object@arguments$tumor_type)
      exclude <- rules[ -druggabilityWhich , , drop=FALSE]
    } else {
      exclude <- rules
    }
    # Check on exclude panel (the one with 8 columns)
    if(! is.null(exclude)){
      if(nrow(exclude)!=0){
        exclude <- .panelCheck(exclude 
          , comparison_panel=panel , tumor_type=NULL)
      }
    }
  }
  ##################################################
  # INITIALIZE CancerPanel Object
  # ------------------------------------------------
  object <- new('CancerPanel')
  message("Calculating panel size...")
  # A precse estimate of space can be 
  # calculate but during simple panel construction
  # we calculate a variation-wise genomic length
  
  #get gene fusions and seperate the gene names
  all_genes <- unique(panel$gene_symbol %>% strsplit(. , "__") %>% unlist)
  
  # ------------------------------------------------
  # Fetch/Calculate feature size
  # ------------------------------------------------
  # get genomic space for the genes of interest
  ann_genes <- .annotateGeneLength(genes=all_genes 
    , canonicalTranscript=canonicalTranscript , myhost=myhost)
  # If we have selected an alteration, correct the size of the feature based 
  # on the selection made.
  panel <- .annotateVariationLength(panel 
    , gene_length=ann_genes , utr=utr , padding_length=padding_length)
    
  # ------------------------------------------------
  # Fetch RS coordinates 
  # ------------------------------------------------
  #If we have RS ids, fetch it
  if(any(panel$exact_alteration=="dbSNP_rs")){
      rs <- unique(panel[ panel$exact_alteration=="dbSNP_rs" 
        , "mutation_specification"])
      #get genomic position from each RS id
      rs_df <- .rsToGenomicPosition(rs)
  } else {
    #if we don't have RS ids, create an empty dataframe
      rs_df <- data.frame(rs="" , genomic_range="" , stringsAsFactors=FALSE)
  }
    
  # ------------------------------------------------
  # Fetch fusion info
  # ------------------------------------------------
  # distinguish between "gene_fusions" and "exact_fusion"
  if(any(panel$alteration=="fusion")){
      panel[ panel$alteration=="fusion" & 
        grepl("__" , panel$gene_symbol) , "exact_alteration"] <- "exact_fusion"
      panel[ panel$alteration=="fusion" & 
        !grepl("__" , panel$gene_symbol) , "exact_alteration"] <- "gene_fusion"
  }
  
  ##################################################
  # ADD INFO FETCHED TO THE newCancerPanel Object
  # ------------------------------------------------
  object@arguments$genedata <- ann_genes
  object@arguments$dbSNP_rs <- rs_df
  object@arguments$panel <- panel
  object@arguments$drugs <- panel$drug[ panel$drug!="" ] %>% unique
  object@arguments$options <- list(padding_length=padding_length 
    , utr=utr , canonicalTranscript=canonicalTranscript)
  if(is.null(rules)){
    object@arguments['rules'] <- list(NULL)
  } else {
    object@arguments[['rules']] <- rules
  }
  return(object)
}