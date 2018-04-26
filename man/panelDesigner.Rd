\name{panelDesigner}
\alias{panelDesigner}
\title{
Estimate panel design according to user specifications
}
\description{
This function takes a CancerPanel object and returns the regions ready to be submitted 
for sequencing
}
\usage{
panelDesigner(object
	, alterationType = c("copynumber", "expression", "mutations", "fusions")
    , padding_length = 100
    , merge_window = 50
    , utr = FALSE 
    , canonicalTranscript=TRUE
    , BPPARAM=bpparam("SerialParam")
    , myhost="www.ensembl.org")
}
\arguments{
  \item{object}{
a CancerPanel Object
}
  \item{alterationType}{
by default, the design of the panel is created by mixing all the different types of alterations. 
With this parameter you can separate the design by alteration type.
}
  \item{padding_length}{
elongation on both side in case of single spot genomic request
}
  \item{merge_window}{
if two ranges are very close to each other what is the minimum length accepted for them to be separated and not merged?
}
  \item{utr}{
if TRUE, the genes ranges in the panel design are taken as CDS plus utr. Default is to take just the coding sequence
}
  \item{canonicalTranscript}{
if FALSE, every exon of every transcript of the gene is taken into consideration in calculating gene length. Default to TRUE is 
to select the canonical transcript
}
  \item{BPPARAM}{
an object of class BiocParallelParam to distribute REST API queries from Ensembl and HGNC. Serialization is the default.  
}
  \item{myhost}{ In case of a biomart breakdown, choose a different host than the default ensembl.org. 
check availability on \href{http://www.ensembl.org/info/about/mirrors.html}{biomart mirrors}
}
}
\details{
In the majority of cases, copynumber and mutations data are retrieved using different technologies and the design
should be separated. Use 'alterationType' parameter to create multiple libraries. 
In case of fusions, the design will take into account all the genes that form the fusion. The technology used 
to find fusion genes can rely on RNA rather than DNA, so in this case it is better to avoid 
this function. A similar idea can be applied for expression data.
'merge_window' parameter is generally calculated by the sequencing company, so set it to 0 if you don't want to decide it upfront. 
The ability of the machine to capture a region and the cost associated with a change in this measure depends on the technology itself. 
It can be very difficult to find the proper trade-off between library size and number of ranges.
The larger is the intronic region accepted, the larger is the library size because you will accept a lot of off-targets. On the other end, 
the more regions in your library, the higher will be the number of amplicons used.
}
\value{
A list of 4 elements: 
\item{GeneIntervals}{ a data.frame containing all gene wise intervals on cds and cds plus utr}
\item{TargetIntervals}{ if the panel contains specific regions, it is a data.frame of non-full gene sequences
                      as requested in the panel }
\item{FullGenes}{ a character vector of genes that will be sequenced for their entire length}
\item{BedStylePanel}{ a bed style data.frame with chromosome start and end of the collapse of GeneIntervals and TargetIntervals}
}
\references{
\href{http://www.ensembl.org/info/website/upload/bed.html}{ bed file format according to ensembl}

\href{http://www.ensembl.org/Help/Glossary?id=346}{ canonical transcript definition according to ENSEMBL}
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{newCancerPanel}}
}
\examples{
#Load panelexample
data(panelexample)
# Create a cancer panel object
mypanel <- newCancerPanel(panelexample)
# Design your panel
# Parallelize part of the code using BiocParallel backend for Unix systems
if(tolower(Sys.info()["sysname"])=="windows"){
  mydesign <- panelDesigner(mypanel 
                          , alterationType="mutations")
} else {
  mydesign <- panelDesigner(mypanel 
                          , alterationType="mutations" 
                          , BPPARAM=BiocParallel::MulticoreParam(workers=2)
                          )
}
# Retrieve bed style sequences
head( mydesign[['BedStylePanel']] )
}

