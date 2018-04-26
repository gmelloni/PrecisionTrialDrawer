\name{cpFreq}
\alias{cpFreq}
\title{
Calculate frequencies of mutations, copynumber, fusion or expression on a CancerPanel object 
}
\description{
Given a CancerPanel object, it returns a data.frame with absolute or 
relative frequencies of alteration per gene
}
\usage{
cpFreq(object
, alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
, mutations.specs=c(NA ,"mutation_type","amino_acid_change"
                    ,"amino_position","genomic_position")
, fusions.specs=c("bygene" , "byfusionpair")
, tumor_type=NULL
, tumor.weights=NULL
, tumor.freqs=NULL
, collapseMutationByGene=TRUE 
, freq=c("relative" , "absolute")
)
}
\arguments{
  \item{object}{
A CancerPanel object filled with genomic data.
}
  \item{alterationType}{
A character vector containing one of the following: 
"copynumber", "expression", "mutations", "fusions".
}
  \item{mutations.specs}{
If alterationType is mutations, this parameters allows further subsetting of frequency of mutation. See details
}
  \item{fusions.specs}{
If "bygene", frequencies of fusion events are calculated by gene, if "byfusionpair" , the calculation are run 
by the pair gene1__gene2
}
  \item{tumor_type}{
A character vector of tumor types to include in the plot among the one included in the object
}
  \item{tumor.weights}{
A named vector of integer values containing an amount of samples to be randomly sampled from the data. 
Each element should correspond to a different tumor type and is named after its tumor code. See details
}
  \item{tumor.freqs}{
A named vector of values between 0 and 1 which sum 1. It contains the expected proportion of patients that are planned 
to be recruited. See Details
}
  \item{collapseMutationByGene}{
A logical that collapse all mutations on the same gene for a single patient as 
a single alteration.
}
  \item{freq}{
return the absolute number of samples with an alteration or the relative number according to the selected cohort
}
}
\details{
This simple frequency calculator allows an exploratory analysis on the frequency of alteration by gene. 
Using \code{mutations.specs}, the user can also calculate the mutation frequency by:
\enumerate{
  \item{mutation_type: }{ each gene is divided by missense, frameshift, splice site etc.}
  \item{amino_acid_change: }{ each gene is divided by amino acid change (V600E , R258H, etc.)}
  \item{amino_position: }{ each gene is stratified by aminoacidic position (600, 258, etc.)}
  \item{genomic_position: }{ each gene is subdivided by exact mutation at genomic level (1:10000:A,C etc.)}
}
Both \code{mutations.specs} and \code{fusions.specs} sort no effect for copynumber and expression.
The table reported for copynumber and expression is different from mutations and fusions. It reports, for every gene, 
the relative or absolute number of patients with "amplification" , "deletion" and "normal" CNA 
and "up" , "down" "normal" level of expression.

By default, \code{cpFreq} will use all the available data from the object, using all the samples for the requested alterationTypes. Nevertheless, one could 
be interested in creating a compound design that is composed by a certain number of samples per tumor type. This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and your design can be stopped when the desired sample size for a given tumor type is reached. 
By adding tumor.weights, we can achieve such target (see examples). Unfortunately, there are two main drawbacks in doing so:
\enumerate{
  \item{small sample size:}{ by selecting small random samples, the real frequency can be distorted. 
                            to avoid this, it is better to run several small samples and then bootstrap them}
  \item{recycling:}{ if the sample size for a tumor type requested by the user is above the available number of cBioportal samples, 
                  the samples are recycled. This has the effect of stabilizing the frequencies but y_measure = "absolute" will have no real meaning when 
                  the heterogeneity of the samples is lost.}
}

A user balanced design can be also obtained using \code{tumor.freqs} parameter. In this case the fraction of altered samples are 
first calculated tumor-wise and then reaggregated using the weights provided by \code{tumor.freqs}. 
If the fraction of altered samples are 0.3 and 0.4 for breast cancer and lung cancer respectively, if you set tumor.freqs = c(brca=0.9 , luad=0.1), 
the full design will have a frequency equal to 0.3*0.9 + 0.4*0.1 = 0.31, that is basically equal to the one of breast samples. 
If this parameter is not set, the total amount of samples available is used with unpredictable balancing. Note that the result can only be expressed as 
relative frequency.

Both tumor.freqs and tumor.weights can achieve a balanced design according to user specification. For having a quick idea of the sample size required, 
it is better to use the former. For having an idea about the possible distribution of sample size giving a finite number of samples
it is better to run the function with tumor.weights several times and aggregate the results. By default, \code{survPowerSampleSize} will use all the available data from the object, using all the samples for the requested alterationTypes. Nevertheless, one could 
be interested in creating a compound design that is composed by a certain number of samples per tumor type. This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and your design can be stopped when the desired sample size for a given tumor type is reached. 
By adding tumor.weights, we can achieve such target (see examples). Unfortunately, there are two main drawbacks in doing so:
\enumerate{
  \item{small sample size:}{ by selecting small random samples, the real frequency can be distorted. 
                            to avoid this, it is better to run several small samples and then bootstrap them}
  \item{recycling:}{ if the sample size for a tumor type requested by the user is above the available number of cBioportal samples, 
                  the samples are recycled. This has the effect of stabilizing the frequencies but y_measure = "absolute" will have no real meaning when 
                  the heterogeneity of the samples is lost.}
}

A user balanced design can be also obtained using \code{tumor.freqs} parameter. In this case the fraction of altered samples are 
first calculated tumor-wise and then reaggregated using the weights provided by \code{tumor.freqs}. 
If the fraction of altered samples are 0.3 and 0.4 for breast cancer and lung cancer respectively, if you set tumor.freqs = c(brca=0.9 , luad=0.1), 
the full design will have a frequency equal to 0.3*0.9 + 0.4*0.1 = 0.31, that is basically equal to the one of breast samples. 
If this parameter is not set, the total amount of samples available is used with unpredictable balancing.

Both tumor.freqs and tumor.weights can achieve a balanced design according to user specification. For having a quick idea of the sample size required, 
it is better to use the former. For having an idea about the possible distribution of sample size giving a few samples (for example a minimum and a maximum sample size) 
it is better to run the function with tumor.weights several times and aggregate the results.
}
\value{
the function returns a object of class data.frame
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{coveragePlot}} \code{\link{coverageStackPlot}}
}
\examples{
# Load example CancerPanel object
data(cpObj)

# Calculate relative frequencies of mutations by gene in breast cancer
cpFreq(cpObj , alterationType="mutations", tumor_type="brca")

# Calculate relative frequencies of mutations by gene and amino_acid_change 
# in both breast and lung cancer
cpFreq(cpObj , alterationType="mutations"
  , tumor_type=NULL 
  , mutations.specs="amino_acid_change")
  
# Calculate the absolute number of samples with amplified, deleted or normal gene
# Using lung cancer available data
cpFreq(cpObj , alterationType="copynumber" , tumor_type="luad" , freq="absolute")

# Calculate frequencies of fusion pairs in all tumor types
cpFreq(cpObj , alterationType="fusions" , fusions.specs="byfusionpair")

# Now calculate mutation freq by gene using 90\% of luad and 10\% of brca samples
cpFreq(cpObj , alterationType="mutations", tumor.freqs=c(brca=0.9 , luad=0.1))
}