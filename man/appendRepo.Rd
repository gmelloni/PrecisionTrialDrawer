\name{appendRepo}
\alias{appendRepo}
\title{
appendRepo: Add extra population data
}
\description{
Function that allows to add extra samples for running simulations. If the user 
wants to add custom data to an existing cBioportal dataset, this function allows to append 
new samples to the panel to increase its inference power. This function takes care 
of appending to the panel a list of alterations compatible with the 
@dataFull slot (the slot that stores the cBioportal downloaded info).
}
\usage{
appendRepo(object, repos)
}
\arguments{
  \item{object}{
An instance of class CancerPanel
}
  \item{repos}{
A list() in format compatible to @dataFull slot
}
}
\details{
dataFull slot is a named list of 4 elements: "fusions" , "mutations" , "copynumber" , "expression".
Each element is a list itself of two named elements: "data" (a data.frame), "Samples" (a list of character vectors). 
They must be all present, but they can be NULL if no data are provided. 
Check example to check correct format of colnames inside each data type.

At the end of the appending, subsetAlterations is automatically run.
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
\code{\link{getAlterations}} \code{\link{subsetAlterations}}
}
\examples{
# Retrieve example data
data(cpObj)
# Check format of slot dataFull
str(cpData(cpObj))
# Create two new mutations
newmutations <- data.frame(
    entrez_gene_id=c(7157 , 7157)
    ,gene_symbol=c("TP53" , "TP53")
    ,case_id=c("new_samp1" , "new_samp2")
    ,mutation_type=c("Nonsense_Mutation" , "Nonsense_Mutation")
    ,amino_acid_change=c("R306*" , "Y126*")
    ,genetic_profile_id=c("mynewbreaststudy" , "mynewbreaststudy")
    ,tumor_type=c("brca" , "brca")
    ,amino_position=c(306 , 126)
    ,genomic_position=c("17:7577022:G,A" , "17:7578552:G,C")
    ,stringsAsFactors=FALSE
    )
newsamples <- c("new_samp1" , "new_samp2")
# A dataFull slot style list should look like this
toBeAdded <- list(fusions=list(data=NULL , Samples=NULL)
                 , mutations=list(data=newmutations , Samples=list(brca=newsamples))
                 , copynumber=list(data=NULL , Samples=NULL)
                 , expression=list(data=NULL , Samples=NULL)
                 )
cpObjplus <- appendRepo(cpObj , toBeAdded)
}
