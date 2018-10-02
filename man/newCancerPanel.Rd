\name{newCancerPanel}
\alias{newCancerPanel}
\title{
CancerPanel object constructor
}
\description{
Given a \code{data.frame} with your panel specifications, it creates a \code{CancerPanel}
object to be used for both simulations and panel design.
}
\usage{
newCancerPanel(panel
            , rules=NULL
            , padding_length = 100
            , utr = FALSE 
            , canonicalTranscript=TRUE
            , myhost="www.ensembl.org"
            )
}
\arguments{
  \item{panel}{ A data.frame describing your panel with alterations, associated drugs and relative genes involved.}
  \item{padding_length}{ An integer that defines how much to extend the targeted regions. \code{paddind_lenght} value 
                   is subtracted from the start coordinate of the targeted region in the panel and, at the same time, 
                   it is also added to the end coordinate of the panel. The result is an extension of the target region by the 
                   \code{padding_length} value in both 5' and 3' direction.}
  \item{utr}{ If \code{TRUE}, the genomic coordinates will also include the UTR regions.}
  \item{canonicalTranscript}{ if FALSE, every exon of every transcript of the gene is taken into consideration 
                   in calculating gene length. Default to TRUE is to select the canonical transcript (see references)}
  \item{myhost}{ In case of a biomart breakdown, choose a different host than the default ensembl.org. 
                 check availability on \href{http://www.ensembl.org/info/about/mirrors.html}{biomart mirrors}}
  \item{rules}{ a data.frame similar to the panel that implement a set of rules so that specific 
  associations between genes/mutations/actionability are overwritten. See Details}
}
\details{
This constructor accepts a \code{data.frame}, \code{tibble} or \code{data.table} with the following columns: 
\describe{
   \item{\code{drug} - character vector}{ A character vector of drug names or drug compounds. 
   It is required, but it can be also filled with \code{NA} if no compound is associated 
   with the alterations or the user is not interested in this feature.}
   \item{\code{gene_symbol} - character vector}{ A character vector of HGNC gene symbols. 
   In case of specific fusion gene, the format is 'gene1__gene2'.}
   \item{\code{alteration} - character vector}{ A character vector with one of the following values: 
   \code{SNV} (Single Nucleotide Variation) or \code{CNA} (Somatic Copy Number Alteration), 
   \code{expression} (up or down gene expression), \code{fusion} 
   (hybrid gene formed from two previously separated genes). 
   This represents the class of alteration. } 
   \item{\code{exact_alteration} - character vector}{ According to the alteration column, 
   it can be one of the cases described in the table below.}
   \item{\code{mutation_specification} - character vector}{ This column refines the location 
   of the alteration type defined in the previous column. 
   In case the record in \code{alteration} is set to 'SNV', the location of the mutation 
   must be specified according to the available options shown in the table below. 
   In case of an alteration type different from 'SNV' the value must be left blank ("") or NA. } 
   \item{\code{group} - character vector  (Optional)}{ A character vector describing 
    a specific group for the alteration. In the \code{\link{panelexample}}, 
    we use it to divide the alterations between druggable (Actionable) and non druggable (Driver).  
    Another possible use is to perform comparisons between different panels.}
}
Possible values for a cancer panel to specify the alterations are: 
\tabular{lll}{
   \bold{alteration} \tab \bold{exact_alteration} \tab \bold{mutation_specification}\cr 
   CNA \tab amplification \tab ""\cr 
   CNA \tab deletion \tab ""\cr 
   expression \tab up \tab ""\cr 
   expression \tab down \tab ""\cr 
   SNV \tab "" \tab ""\cr 
   SNV \tab mutation_type \tab missense\cr 
   SNV \tab mutation_type \tab truncating\cr 
   SNV \tab amino_acid_position \tab 300-300\cr 
   SNV \tab amino_acid_position \tab 300-350\cr 
   SNV \tab amino_acid_variant \tab V600E\cr 
   SNV \tab genomic_position \tab 13:20000-40000\cr 
   SNV \tab genomic_position \tab 13:20000-20000\cr 
   SNV \tab genomic_variant \tab 13:20000:A,C\cr 
   SNV \tab dbSNP_rs \tab rs1234567\cr 
   fusion \tab "" \tab ""
}
\code{drug} values NA or empty are transformed to "no_drug", that is a reserved value. 
\code{group} values NA or empty are transformed to "no_group", that is a reserved value.
\code{gene_symbol} is mandatory for any alteration type because alterations are retrieved from cBioPortal using this key. 
In case of rs IDs, the closest gene symbol is the ideal annotation, even if the position is intergenic and generally not 
associated with any gene.

\code{rules} parameter implements a set of negation rules.

A data.frame like the following can be used in the rules parameter. 
The data.frame is the same as above but it adds tumor_type and in_out columns: 
\tabular{llllllll}{
   \bold{drug} \tab \bold{gene_symbol} \tab \bold{alteration} \tab \bold{exact_alteration} \tab \bold{mutation_specification} \tab \bold{group} \tab \bold{tumor_type} \tab \bold{in_out}\cr 
   Erlotinib \tab EGFR \tab SNV \tab amino_acid_variant \tab T790M \tab Driver \tab  \tab exclude\cr
   Erlotinib \tab KRAS \tab SNV \tab  \tab  \tab Driver \tab  \tab exclude \cr
   Erlotinib \tab      \tab     \tab  \tab  \tab Driver \tab luad \tab include \cr
   Olaparib  \tab      \tab     \tab  \tab  \tab Driver \tab brca \tab include \cr
}

The new column \code{tumor_type} can contain a single tumor_type code 
or an empty string (which means that the rule is valid for any tumor type). 
The other new column is \code{in_out} which can only contained 'include' or 'exclude' values.

The first two rows implement a resistance rule. 
Any sample with a T790M mutaion on EGFR or a KRAS mutation 
cannot be associated with Erlotinib because it generates resistance to EGFR inhibitors.
The effect is that every sample with either a EGFR T790M mutation 
or any KRAS mutation will no longer be associated with Erlotinib (it will be considered 'no_drug') 
but the mutation will not be filtered out. In the group slot, 
the user can put what group the patient with that drug association 
will be changed into (in the example from Actionable to Driver). 
We set the rule as exclusion (in_out=exclude) for any tumor type.

The other two rows are more stringent and contain a drug inclusion/exclusion rule. 
Erlotinib can only be associated with luad and Olaparib only to brca (in_out=include). 
No matter what tumor type will be used in the future, these rules will always be applied so that 
any tumor type that is not luad or brca will be excluded.

}
\value{
This method returns a \code{CancerPanel} object with the slot arguments updated.
}
\references{
\href{www.biomart.org}{ source of gene length and exon structure}

\href{http://www.genenames.org}{ source of official gene symbols and mapping with ensembl}

\href{http://www.ensembl.org/Help/Glossary?id=346}{ canonical transcript definition according to ENSEMBL}

\href{http://www.ensembl.org/info/about/mirrors.html}{ biomart mirrors }
}
\author{
Giorgio Melloni , Alessandro Guida
}
\seealso{
\code{\link{panelDesigner}}
}
\examples{
# Load the panel example
data(panelexample)
# Create a CancerPanel object
mypanel <- newCancerPanel(panelexample , canonicalTranscript=FALSE)
}

