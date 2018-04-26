\name{filterFusions}
\alias{filterFusions}
\title{
filterFusions: remove or keep specified fusions
}
\description{
This method allows to keep only or to exclude certain fusions according to three different filtering formats
}
\usage{
filterFusions(object , filtered , mode = c("exclude" , "keep"))
}
\arguments{
  \item{object}{
An instance of class CancerPanel
}
  \item{filtered}{
A character vector used as filter
}
  \item{mode}{
If "exclude", fusions are removed from the object. 
If "keep", fusions specified in filtered are the only ones maintained.
}
}
\details{
filtered vector must be in fusion format GENE1\_\_GENE2.

At the end of the filtering procedure, subsetAlterations is automatically run.
}
\value{
An updated instance of class CancerPanel
}
\references{
\href{www.cbioportal.org}{ data origin for mutations, copynumber and expression data}

\href{www.tumorfusions.org}{ data origin for fusion data}
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{getAlterations}} \code{\link{subsetAlterations}} \code{\link{filterMutations}}
}
\examples{
# Retrieve example data
data(cpObj)
# Create a data.frame to filter
myFilter <- c("ERC1__RET", "TRIM33__RET", "EML4__ALK")
# Keep only myFilter fusions
cpObjKeep <- filterFusions(cpObj , filtered = myFilter , mode = "keep")
# Exclude myFilter fusions
cpObjExclude <- filterFusions(cpObj , filtered = myFilter , mode = "exclude")
}