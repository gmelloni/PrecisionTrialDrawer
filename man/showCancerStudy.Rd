\name{showCancerStudy}
\alias{showCancerStudy}
\title{
List of cancer studies included in the cBioportal
}
\description{
Show all the possible cancer studies accepted by PrecisionTrialDrawer
}
\usage{
showCancerStudy(tumor_type=NULL)
}
\arguments{
\item{tumor_type}{ a valid tumor type barcode in TCGA standard coding}
}
\value{
A named vector of all the cancer studies available in cBioPortalData package 
that can be passed to the method getAlterations.
}
\details{
This method is a wrapper around \code{\link[cBioPortalData]{getStudies}} 
and show all the cancer studies included in the cBioPortal. 
If you specify the tumor_type, follow 
the list of \code{showTumorType}.
Even though the cancer studies are present, 
it doesn't mean that there are data for every alteration. 
A tumor type could have mutation data but not copynumber 
or viceversa. Furthermore, fusion data 
have a different source and could not be included in the list.
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{getAlterations}} \code{\link{showTumorType}} 
\code{\link[cBioPortalData]{getStudies}}
}
\examples{
myCanStudies <- showCancerStudy()
head(myCanStudies)
}