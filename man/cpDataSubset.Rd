\name{cpDataSubset}
\alias{cpDataSubset}
\title{
Method that returns the \code{dataSubset} slot in a \code{CancerPanel} object.
}
\description{
This method returns a list with an element for each alteration type. 
It shows the content of the slot \code{dataSubset} in the \code{CancerPanel} object.
}
\usage{
cpDataSubset(object)
}
\arguments{
  \item{object}{ a CancerPanel object}
}
\details{
The slot showed with this getter is a list of length 4 even if no data 
for a certain alteration type were requested. In case no data were 
requested or if there are no data for a certain tumor type, 
the slot is filled with a \code{NULL} value. 
If data are  present but there are no alteration for the 
specified tumor types the slot is filled with a 0 rows \code{data.frame}.
}
\value{
A list of 4 elements, one for each alteration type: 
\item{mutations }{ a data.frame containing exactly the mutation alterations requested in the panel.}
\item{copynumber }{ a data.frame containing exactly the copynumber alterations requested in the panel.}
\item{expression }{ a data.frame containing exactly the expression alterations requested in the panel.}
\item{fusions }{ a data.frame containing exactly the fusions alterations requested in the panel.}
}
\references{
\href{www.cbioportal.org}{ data origin for mutations, copynumber and expression data}

\href{www.tumorfusions.org}{ data origin for fusion data}
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{cpData}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Show slot dataSubset
str( cpDataSubset(cpObj) )
}

