\name{coveragePlot}
\alias{coveragePlot}
\title{
A series of bar charts representing the number of samples harbouring at least
one or more alterations
}
\description{
Given a CancerPanel object, it returns one or more bar charts representing the number of samples 
covered by at least 1, 2 ,3 or more alterations using specified data. Each plot is controlled by 
the grouping parameter.
}
\usage{
coveragePlot(object
, alterationType = c("copynumber", "expression", "mutations", "fusions")
, grouping = c(NA,"drug","group","gene_symbol","alteration_id","tumor_type")
, tumor_type=NULL
, alterationType.agg=TRUE
, collapseMutationByGene = TRUE
, collapseByGene = FALSE
, tumor.weights=NULL
, tumor.freqs=NULL
, maxNumAlt = 10
, colNum=NULL
, cex.main="auto"
, noPlot = FALSE)
}
\arguments{
  \item{object}{
A CancerPanel object filled with genomic data.
}
  \item{alterationType}{
A character vector containing one or more of the following: 
"copynumber", "expression", "mutations", "fusions".
}
  \item{grouping}{
A character vector containing one or more of the following: 
NA, "drug", "group", "gene_symbol", "alteration_id", "tumor_type".
}
  \item{tumor_type}{
A character vector of tumor types to include in the plot among the one included in the object
}
  \item{alterationType.agg}{
logical value. If TRUE, the default, the frequencies displayed are calculated over all the samples that were tested 
for all the alterationType requested. If FALSE all the samples tested for the specified alteration_id stratum are used. 
It sorts an effect if 'alterationType' length is > 1 and 'alteration_id' is in grouping parameter. See details.
}
  \item{collapseMutationByGene}{
A logical that collapse all mutations on the same gene for a single patient as 
a single alteration.
}
  \item{collapseByGene}{
A logical that collapse all alterations on the same gene for a single patient as 
a single alteration. e.g. if a sample has TP53 both mutated and deleted as copynumber, it 
will count for one alteration only.
}
  \item{tumor.weights}{
A named vector of integer values containing an amount of samples to be randomly sampled from the data. 
Each element should correspond to a different tumor type and is named after its tumor code. See details
}
  \item{tumor.freqs}{
A named vector of values between 0 and 1 which sum 1. It contains the expected proportion of patients that are planned 
to be recruited. See Details
}
  \item{maxNumAlt}{
This number represents the maximum number on X axis.
}
  \item{colNum}{
If set, represents the number of columns in plotting layout. If NULL, best square 
representation is chosen instead.
}
  \item{cex.main}{
a numerical value or "auto". This parameter can set the size of each plot main title. 
Default is "auto", for automatic resizing.
}
  \item{noPlot}{
If TRUE, the plot is not shown but just the data used to drawn it.
}
}
\details{
According to the chosen \code{alterationType}, the package will look 
for all the samples with available data for all the selected \code{alterationType}. 
For example, if alterationType = c( "mutations" , "copynumber"), 
only common samples with both mutation and copynumber data are used by default. If \code{alterationType.agg} is 
set to FALSE and "alteration_id" is in grouping, the default behaviour changes. "mutations" plot will be 
displayed with the frequencies relative to all the samples tested for mutations and "copynumber" with all 
the samples tested for CNA.
If "tumor_type" is in grouping variable, each plot is evaluated on the samples relative to the tumor type.
The number of plots depends on the multiplication of the levels of the grouping variable. 
If you put too many grouping variable, it is better to draw a \code{\link{coverageStackPlot}} 
or to redirect the output to a file.

By default, \code{coveragePlot} will use all the available data from the object, using all the samples for the requested alterationTypes. Nevertheless, one could 
be interested in creating a compound design that is composed by a certain number of samples per tumor type. This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and your trial can be stopped when the desired sample size for a given tumor type is reached. 
By adding tumor.weights, we can achieve such target (see examples). Unfortunately, there are two main drawbacks in doing so:
\enumerate{
  \item{small sample size:}{ by selecting small random samples, the real frequency can be distorted. 
                            to avoid this, it is better to run several small samples and then aggregate the results}
  \item{recycling:}{ if the sample size for a tumor type requested by the user is above the available number of cBioportal samples, 
                  the samples are recycled. This has the effect of stabilizing the frequencies but y_measure = "absolute" will have no real meaning when 
                  the heterogeneity of the samples is lost.}
}

A user balanced design can be also obtained using \code{tumor.freqs} parameter. In this case the fraction of altered samples are 
first calculated tumor-wise and then reaggregated using the weights provided by \code{tumor.freqs}. 
If the fraction of altered samples are 0.3 and 0.4 for breast cancer and lung cancer respectively, if you set tumor.freqs = c(brca=0.9 , luad=0.1), 
the full design will have a frequency equal to 0.3*0.9 + 0.4*0.1 = 0.31, that is basically equal to the one of breast samples. 
If this parameter is not set, the total amount of samples available is used with unpredictable balancing. In the examples, brca and luad data are used. 
Breast samples are at least twice as much as luad samples and tumor.freqs can help with a more balanced simulation.

Both tumor.freqs and tumor.weights can achieve a balanced design according to user specification. To have a quick idea of the sample size required, 
it is better to use the former. For having an idea about the possible distribution of sample size giving a few samples (for example a minimum and a maximum sample size) 
it is better to run the function with tumor.weights several times and aggregate the results to obtain mean values, confidence intervals etc.
}
\value{
If noPlot=FALSE, this method returns a bar chart or a series of bar charts. Y-axis represents the number of samples, X-axis the incremental number of alterations 
per sample. In case tumor.freqs is set, the Y-axis represents the relative frequency that is reported as text on the top of the bars.
If noPlot=TRUE, it returns a named list: 
\item{plottedTable }{a matrix with absolute number of samples plotted. 
Every column represents how many samples retain at least 1, 2, 3 ... alterations. 
Every row is a different plot for one of the specified grouping levels. If tumor.freqs is used, 
relative frequencies are reported instead.
}
\item{Samples }{ a numeric vector corresponding to the rows of 
plottedTable representing the number of reference sample for each plot.
If tumor.freqs is used, Samples is NULL.
}
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{saturationPlot}} \code{\link{coverageStackPlot}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Plot the coverage of this panel by tumor type and drug
# Using mutations and copynumber data
coveragePlot(cpObj , alterationType=c("mutations" , "copynumber") 
          , grouping=c("tumor_type" , "drug")
          , maxNumAlt=5)
# Report data used in the plot above
coveragePlot(cpObj , alterationType=c("mutations" , "copynumber") 
          , grouping=c("tumor_type" , "drug")
          , maxNumAlt=5
          , noPlot=TRUE)

# Observe the effect of alterationType.agg on frequencies
coveragePlot(cpObj , alterationType=c("mutations" , "copynumber") 
            , grouping=c("alteration_id")
            , noPlot=FALSE)
coveragePlot(cpObj , alterationType=c("mutations" , "copynumber") 
            , grouping=c("alteration_id")
            , noPlot=FALSE
            , alterationType.agg = FALSE)
# Use tumor.weights. Create a design with 40 random luad samples and 10 brca
coveragePlot(cpObj , alterationType=c("mutations" , "copynumber" , "fusions") 
          , grouping=NA
          , maxNumAlt=5
          , tumor.weights = c(brca=10 , luad=40))
# Use tumor.freqs. Calculate frequencies with 90\% weight to luad samples and 10\% to brca
coveragePlot(cpObj , alterationType=c("mutations" , "copynumber" , "fusions") 
          , grouping=NA
          , maxNumAlt=5
          , tumor.freqs = c(brca=0.1 , luad=0.9))
}