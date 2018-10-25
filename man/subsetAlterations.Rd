\name{subsetAlterations}
\alias{subsetAlterations}
\title{
Retrieve exactly the the alteration requested in 
the panel after a getAlterations method
}
\description{
As soon as all the data are retrieved from local or 
remote repositories, this method simply subsets 
the alterations as exactly requested by the panel 
and put them in a format that is common to all kinds of alterations
}
\usage{
subsetAlterations(object , rules)
}
\arguments{
  \item{object}{
a CancerPanel object
}
  \item{rules}{
a data.frame with the same columns as the panel plus 
tumor_type and in_out. Check newCancerPanel function
}
}
\details{
This method will raise an error if no data are available 
in the slot dataFull. The method is compulsory to draw all plots 
since it attaches to each alteration all the panel characteristics, 
like for example druggability and group.
}
\value{
This method returns the object itself with the slot dataSubset updated. 
This slot is formed by a list of 4 elements: 
\item{mutations }{ a data.frame containing exactly the 
mutation alterations requested in the panel}
\item{copynumber }{ a data.frame containing exactly 
the copynumber alterations requested in the panel}
\item{expression }{ a data.frame containing exactly 
the expression alterations requested in the panel}
\item{fusions }{ a data.frame containing exactly the 
fusions alterations requested in the panel}
\item{excluded }{ a data.frame containing the samples 
that are resistant to a drug, if rules is set}

Each of the first 4 element is a data.frame with the following columns:
\item{drug}{ drug associated with the alteration, if any}
\item{group}{ group associated with the alteration, if present}
\item{gene_symbol}{ gene associated with the alteration. 
In case of fusions, a fusion gene in the format gene1__gene2}
\item{tumor_type}{ tumor_type in which the alterations were found}
\item{case_id}{ sample IDs in which the alterations were found}
\item{alteration_id}{ type of alteration with autonumeric 
identifier. e.g. a specific fusion is fus_1, fus_2 etc.}

Please note that the number of alterations per alteration type 
can have more rows than the original dataFull slot. This is 
because the same alteration can be a target of multiple drugs 
or belong to more than one group and is therefore repeated.
}
\references{
\href{www.cbioportal.org}{ data origin for mutations, 
copynumber and expression data}

\href{www.tumorfusions.org}{ data origin for fusion data}
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{getAlterations}} \code{\link{newCancerPanel}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Launch subsetAlterations excluding certain mutations 
# from being considered actionable
rules <- data.frame(
    drug=c("Erlotinib" , "Erlotinib", "Erlotinib","Erlotinib","Olaparib")
    , gene_symbol=c("EGFR" , "KRAS", "" , "", "")
    , alteration=c("SNV" , "SNV", "" , "", "")
    , exact_alteration=c("amino_acid_variant" , "", "" , "", "")
    , mutation_specification=c("T790M" , "" , "" , "", "")
    , group=c("Driver" , "Driver", "Driver" , "Driver", "Driver")
    , tumor_type=c("luad" , "luad" , "luad" , "coadread","brca")
    , in_out=c("exclude" , "exclude" , "include" , "include" , "include")
    , stringsAsFactors = FALSE)
cpObj <- subsetAlterations(cpObj , rules = rules)
# See the updated slot
str(cpDataSubset(cpObj))
}

