\name{cpArguments}
\alias{cpArguments}
\title{
Method that returns the \code{arguments} slot in a \code{CancerPanel} object.
}
\description{
This method returns a list with the panel 
general information after the creation of the object.
}
\usage{
cpArguments(object)
}
\arguments{
  \item{object}{ a CancerPanel object}
}
\details{
The length of the output list is always defined by 5 elements. 
For example, if there are no rs numbers present in the panel, 
the \code{data.frame} in dbSNP_rs element will be simply empty. 
The first 4 elements are filled at the very creation of the 
object with \code{newCancerPanel}. The \code{tumor_type} element 
is filled after data request with \code{getAlterations}.
}
\value{
A list of 5 elements, one for each alteration type: 
\item{genedata }{a data.frame with the cds and 
cds plus utr of the requested genes.}
\item{dbSNP_rs }{ a data.frame with all the rs 
number translated in genomic hg19 positions.}
\item{panel }{ a data.frame with the original 
panel information plus the length of each alteration.}
\item{drugs }{ a character vector containing unique 
drug names as reported in the panel.}
\item{tumor_type }{ a character vector with the 
tumor types requested in the panel.}
}
\references{
\href{www.biomart.org/}{ Source of gene length and exon structure.}
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{cpData}} \code{\link{cpDataSubset}}
}
\examples{
# Load example CancerPanel object
data(cpObj2)
# Show slot dataSubset
str(cpArguments(cpObj2) )
}

