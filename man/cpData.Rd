\name{cpData}
\alias{cpData}
\title{
Method that returns the \code{dataFull} slot in a \code{CancerPanel} object.
}
\description{
This method returns a list with an element for each alteration type contained 
in the class \code{CancerPanel} slot \code{dataFull}
}
\usage{
cpData(object)
}
\arguments{
  \item{object}{ a CancerPanel object}
}
\details{
The length of the output list is always defined by 4 elements even if no data, 
for a certain alteration type, were requested.In case no data were requested or 
if there are no data for a certain tumor type, the slot is filled with \code{NULL} 
values. Every element is a list of 2 elements: 
\itemize{
 \item{\code{data}}{ A \code{data.frame} in a format specific for the alteration type.}
 \item{\code{Samples}}{ The second element is a \code{list} of vectors containing the names of the samples for each tumor type.}
}
}
\value{
 A list of 4 elements, one for each alteration type: 
 \item{mutations}{ A list of 2 elements containing the mutation alterations for each gene requested in the panel.}
 \item{copynumber }{ A list of 2 elements containing the copynumber alterations for each gene requested in the panel.}
 \item{expression }{ A list of 2 elements containing the expression alterations for each gene requested in the panel.}
 \item{fusions }{ A list of 2 elements containing the fusions alterations for each gene requested in the panel.}
}
\references{
\href{www.cbioportal.org}{ data origin for mutations, copynumber and expression data}
\href{http://54.84.12.177/PanCanFusV2}{ data origin for fusion data}
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{cpDataSubset}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Show slot dataSubset
str( cpData(cpObj) )
}

