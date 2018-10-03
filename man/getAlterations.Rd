\name{getAlterations}
\alias{getAlterations}
\title{
Retrieve genomic data for each gene in the panel
}
\description{
This method updates the CancerPanel object with data from cBioportal and MD Anderson fusion database
}
\usage{
getAlterations(object, tumor_type = NULL, repos = NULL
  , mutation_type = c("all_nonsynonymous" , "all_mutations" , "missense" , "truncating")
  , expr_z_score = 2
  , BPPARAM = bpparam("SerialParam")
  , gene_block=50)
}
\arguments{
  \item{object}{ A CancerPanel object}
  \item{tumor_type}{ a vector of tumor types from the ones available using \code{showTumorType} or in the first column of \code{showCancerStudy}. See details}
  \item{repos}{ a list containing custom data to be used in the object}
  \item{mutation_type}{ decide the kind of mutations to retrieve. Only non synonymous, only missense, only truncating or all the mutations. The default is to retrieve only non synonymous}
  \item{expr_z_score}{ a number that expresses the threshold at which a gene is considered upregulated or downregulated}
  \item{BPPARAM}{ parameter for \code{bplapply} to parallelize part of the code}
  \item{gene_block}{ Set how many genes at a time are requested to cBioPortal. Default 50}
}
\details{
This method fills the slot dataFull in a cancer panel object. It retrieves data by gene 
from cBioportal and MD Anderson fusion database according to the specifications of the panel.
This slot is composed by a list of 4 elements, one for each alterationType: 
'fusions' , 'mutations' , 'copynumber' , 'expression'. Every element is a list of 2 elements: data , Samples.
The first element is a data.frame in a format specific for the alteration type. The second element 
is a list of vectors containing the names of all the samples analyzed for each tumor type (both altered and "wild-type"). 

\code{tumor_type} parameter can be either a list of tumor types ( 'brca' , 'luad') or specific cancer studies 
( 'brca_tcga_pub2015' , 'luad_tcga_pub' ) but not the two things together. Check availability with the functions mentioned above. 
In case of cancer studies selection, fusions are retrieved from cancer type definition 
( 'brca_tcga_pub2015' becomes simply 'brca') since they come from a different source.

The expression value are expressed as up or down according to the threshold in 'expr_z_score'. A gene is considered 
upregulated if its z-score is over 2 or downregulated if is lower than -2.

The copynumber values are expressed as 'amp' or 'del' according to GISTIC definition. A gene is reported as amplified 
or deleted if its value after GISTIC evaluation is 2 or -2.

A message about the current tumor in download is prompted for every study. If the message appears more than once, 
it probably means that the genes you requested were more than 100 and so the query was actually split in two or 
more chunks to avoid an overload on cbioportal database.
}
\value{
The method returns the original CancerPanel object with the slot dataFull updated.
}
\references{
\href{www.cbioportal.org}{ data origin for mutations, copynumber and expression data}

\href{www.tumorfusions.org}{ data origin for fusion data}
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{subsetAlterations}} \code{\link{showTumorType}} \code{\link{showCancerStudy}}
}
\examples{
#Load panelexample
data(cpObj2)
# Retrieve data from AML
cpObj <- getAlterations(cpObj2 , tumor_type=c("laml"))
}

