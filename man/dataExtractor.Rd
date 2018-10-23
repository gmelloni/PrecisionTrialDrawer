\name{dataExtractor}
\alias{dataExtractor}
\title{
dataExtractor: extract data from a CancerPanel object
}
\description{
Given user specified options, return specific data from a CancerPanel object, 
including alteration and Samples
}
\usage{
dataExtractor(object 
    , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
    , tumor_type=NULL 
    , collapseMutationByGene=TRUE 
    , collapseByGene=FALSE 
    , tumor.weights=NULL)
}
\arguments{
    \item{object}{
An instance of class CancerPanel
}
    \item{alterationType}{
what kind of alteration to include. It can be one or more between 
"copynumber", "expression", "mutations", "fusions".
Default is to include all kind of alterations.
}
    \item{tumor_type}{
only plot one or more tumor types among 
the ones available in the object.
}
    \item{collapseMutationByGene}{
A logical that collapse all mutations on 
the same gene for a single patient as 
a single alteration.
}
    \item{collapseByGene}{
A logical that collapse all alterations 
on the same gene for a single patient as 
a single alteration. e.g. if a sample has 
TP53 both mutated and deleted as copynumber, it 
will count for one alteration only.
}
    \item{tumor.weights}{
A named vector of integer values containing an 
amount of samples to be randomly sampled from the data. 
Each element should correspond to a different 
tumor type and is named after its tumor code. See details
}
}
\details{
This function is used internally by most of the methods of 
the package and provide a common data extractor for 
a CancerPanel object. It is a low level function to retrieve 
data for other custom usages, in particular via tumor.weights.
}
\value{
A named list with data, samples and tumors not 
present in the CancerPanel object is returned.
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{getAlterations}} \code{\link{subsetAlterations}}
}
\examples{
# Retrieve example data
data(cpObj)
# Extract CNA and mutation data
mydata <- dataExtractor(cpObj 
, alterationType=c("copynumber" , "mutations") 
, tumor_type="brca")
# It is particularly useful for bootstrap simulations
# Here we extract 10 random samples composed by 30 brca and 40 luad
myboot <- replicate(10
, dataExtractor(cpObj
, alterationType="mutations"
, tumor.weights=c("brca"=30 , "luad"=40)
))
}
