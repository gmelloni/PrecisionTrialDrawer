\name{filterMutations}
\alias{filterMutations}
\title{
filterMutations: remove or keep specified mutations
}
\description{
This method allows to keep only or to exclude certain mutations 
according to three different filtering formats
}
\usage{
filterMutations(object 
, filtered=NULL 
, bed = NULL 
, mode = c("exclude" , "keep") 
, tumor_type=NULL)
}
\arguments{
  \item{object}{
An instance of class CancerPanel
}
  \item{filtered}{
A data.frame used as filter
}
  \item{bed}{
A data.frame in bed format (chr start end)
}
  \item{mode}{
If "exclude", mutations are removed from the object. 
If "keep", mutations specified in filtered are the only ones maintained.
}
  \item{tumor_type}{
A vector of tumor_type names. The filter will 
be active only on the specified tumor types
}
}
\details{
filtered data.frame can come in three different formats:
\enumerate{
    \item{ "gene_symbol" , "amino_position"}{ Ex. BRAF 600}
    \item{ "gene_symbol" , "amino_acid_change"}{ Ex. BRAF V600E}
    \item{ "genomic_poistion"}{ Ex. 3:1234567:A,C}
}

bed file must be composed by 3 columns: chrN, start 0-base , end 1-base

At the end of the filtering procedure, subsetAlterations is automatically run.
}
\value{
An updated instance of class CancerPanel
}
\references{
\href{www.cbioportal.org}{ data origin for mutations
, copynumber and expression data}

\href{www.tumorfusions.org}{ data origin for fusion data}
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{getAlterations}} \code{\link{subsetAlterations}} 
\code{\link{filterFusions}}
}
\examples{
# Retrieve example data
data(cpObj)
# Create a data.frame to filter
myFilter <- data.frame(gene_symbol = c("BRAF" , "PIK3CA") 
    , amino_position = c(600 , 118))
# Keep only myFilter mutations
cpObjKeep <- filterMutations(cpObj , filtered = myFilter , mode = "keep")
# Exclude myFilter mutations
cpObjExclude <- filterMutations(cpObj , filtered = myFilter , mode = "exclude")

# Create a bed file
myBed <- data.frame(chr = paste0("chr" , c(7 , 17)) 
    , start = c(140534632 , 41244326) 
    , end = c(140534732 , 41244426) , stringsAsFactors=FALSE)
# Keep only myFilter mutations
cpObjKeep <- filterMutations(cpObj , bed = myBed , mode = "keep")
# Exclude myBed mutations
cpObjExclude <- filterMutations(cpObj , bed = myBed , mode = "exclude")
}
