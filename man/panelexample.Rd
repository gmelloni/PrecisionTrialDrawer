\name{panelexample}
\alias{panelexample}
\docType{data}
\title{
Example of a cancer panel
}
\description{
A data frame containing a cancer panel data.frame to create a CancerPanel object
}
\usage{
data(panelexample)
}
\format{
A data.frame of 6 columns:
\enumerate{
    \item{drug: }{ character vector of drug names or any chemical identifiers}
    \item{gene_symbol: }{ character vector HGNC official gene symbol}
    \item{alteration: }{ a character vector of accepted alteration type, 
                        'SNV', 'CNA', 'expression', 'fusion'}
    \item{exact_alteration: }{ a character vector that 
    identifies the exact alteration 
    depending on alteration value. For example, alteration 'CNA' could
    correspond to exact_alteration 'amplification' or 'deletion'}
    \item{mutation_specification: }{ a character vector that 
    must be empty if alteration 
    is not 'SNV' and identifies the exact alteration in various format}
    \item{group: }{ a character vector that is useful for custom annotation of 
                    each single row in the panel. 
                    In the example is used to identify actionable 
                    from driver variants}
}
}
\examples{
#Load panelexample and show its structure
data(panelexample)
str(panelexample)
}
\keyword{datasets}
