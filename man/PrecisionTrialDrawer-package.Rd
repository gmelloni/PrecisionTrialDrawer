\name{PrecisionTrialDrawer-package}
\alias{PrecisionTrialDrawer-package}
\alias{PrecisionTrialDrawer}
\docType{package}
\title{
PrecisionTrialDrawer: a package to analyze and design custom cancer panels
}
\description{
PrecisionTrialDrawer is a package that a allows 
the design and study of a cancer panel 
to create pilot sets for a clinical trial.
}
\details{
It allows simulations on real TCGA data to design a possible clinical trial, 
both as basket or umbrella design. 

Furthermore, it can provide an easy function to 
calculate genomic space occupied by the cancer 
panel in a commonly used bed format 
ready to be submitted to a NGS company.
}
\author{
Giorgio Melloni, Alessandro Guida
}
\references{
\href{www.cbioportal.org}{ Cbioportal web site}
\href{http://54.84.12.177/PanCanFusV2}{ MD Anderson TCGA Fusion Portal}
}
\keyword{ package }
\examples{
# Retrieve a Cancer Panel 
data(cpObj)
# Show the number of covered samples
coveragePlot(cpObj 
            , alterationType=c("mutations" , "copynumber") 
            , grouping=c("group" , "tumor_type"))
}
