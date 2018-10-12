\name{cpDesign}
\alias{cpDesign}
\docType{data}
\title{
Example of a cancer panel design export
}
\description{
Bed style export of cpObj
}
\usage{
data(cpDesign)
}
\format{
A design is a list composed by four elements:
\describe{
\item{\code{GeneIntervals}:}{all CDS and CDS + UTR for all the genes in the panel}
\item{\code{TargetIntervals}:}{all requested target regions 
(specific single mutations) divided and collapsed by gene symbol}
\item{\code{FullGenes}:}{gene symbols of the genes considered for their full sequence}
\item{\code{BedStylePanel}:}{the entire panel in bed format, merged by chromosome, start and end.}
  }
}
\source{
Derived from running panelDesigner on cpObj
}
\examples{
#Load cpDesign and show its structure
data(cpDesign)
str(cpDesign)
}
\keyword{datasets}

