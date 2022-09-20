\name{showTumorType}
\alias{showTumorType}
\title{
List of tumor type barcode
}
\description{
Show all the possible tumor types accepted by PrecisionTrialDrawer
}
\usage{
showTumorType()
}
\value{
A data.frame of all the tumor types available in cBioPortalData package 
that can be passed to the method getAlterations. Every element is the 
aggregation of all the available sequenced data from all the studies 
involved in a particular tumor type.
}
\details{
This method is a wrapper around \code{\link[cBioPortalData]{getStudies}}
and show all the barcodes for the tumor 
types included in the cBioPortal. 
Even though the tumor types are present, 
it doesn't mean that there are data for every alteration. 
A tumor type could have mutation data but not copynumber or viceversa. 
Furthermore, fusion data 
have a different source and could not be included in the list.
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{getAlterations}} \code{\link{showCancerStudy}} 
\code{\link[cBioPortalData]{getStudies}}
}
\examples{
myTumTypes <- showTumorType()
head(myTumTypes)
}