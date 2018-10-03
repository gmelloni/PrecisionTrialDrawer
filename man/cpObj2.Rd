\name{cpObj2}
\alias{cpObj2}
\docType{data}
\title{
Example of a cancer panel
}
\description{
A CancerPanel object used in PrecisionTrailDrawer
}
\usage{
data(cpObj2)
}
\format{
A small CancerPanel object with 3 slots and mutation data from every 
cBioportal study on Acute Myeloid Leukemia (LAML)
\describe{
    \item{\code{arguments}:}{Object of class \code{"list"} containing panel contents information.
                            it id filled at the panel construction.
                            }
    \item{\code{dataFull}:}{Object of class \code{"list"} containing a slot for every alteration type.
                            it is filled with getAlterations method and contains all the alterations 
                            for all the genes in the panel and every tumor type requested.
                            }
    \item{\code{dataSubset}:}{Object of class \code{"list"} containing a slot for every alteration type.
                            it is filled with subsetAlterations method and reports all the alterations 
                            for every requested specific alteration. All slots have the same data.frame layout.
                            }
  }
}
\source{
\href{www.cbioportal.org}{ data origin for LAML}
}
\examples{
#Load cpObj2 and show its structure
data(cpObj2)
str(cpObj2)
}
\keyword{datasets}
