\name{coverageStackPlot}
\alias{coverageStackPlot}
\title{
A stacked and beside bar chart representing a breakdown of covered samples
}
\description{
Given a CancerPanel object, it returns one bar chart representing the number of samples 
covered by at least 1 alteration under 'var' divided by 'grouping'
}
\usage{
coverageStackPlot(object
, alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
, var=c("drug","group","gene_symbol","alteration_id","tumor_type")
, grouping=c(NA,"drug","group","gene_symbol","alteration_id","tumor_type")
, tumor_type=NULL
, collapseMutationByGene=TRUE
, collapseByGene=TRUE
, tumor.weights=NULL
, tumor.freqs=NULL
, plotFreq = FALSE
, noPlot=FALSE
, html=FALSE)
}
\arguments{
  \item{object}{
A CancerPanel object filled with genomic data.
}
  \item{alterationType}{
A character vector containing one or more of the following: 
"copynumber", "expression", "mutations", "fusions".
}
  \item{var}{
a character vector of length 1 containing one or more of the following: 
"drug", "group", "gene_symbol", "alteration_id" , "tumor_type". 
This parameter is compulsory and decide the classes of the bars.
}
  \item{grouping}{
a character vector of length 1 containing one or more of the following: 
NA, "drug", "group", "gene_symbol", "alteration_id", "tumor_type".
This parameter decide the breakdown of var. 
If not set, it is considered NA and only 'var' is plotted with no stacking.
}
  \item{tumor\_type}{
a character vector containing tumor types to be plotted
}
  \item{collapseMutationByGene}{
A logical that collapse all mutations on the same gene for a single patient as 
a single alteration.
}
  \item{collapseByGene}{
A logical that collapse all alterations on the same gene for a single patient as 
a single alteration. e.g. if a sample has TP53 both mutated and deleted as copynumber, 
it will count for one alteration only.
}
  \item{tumor.weights}{
A named vector of integer values containing an amount of samples to be randomly sampled from the data. 
Each element should correspond to a different tumor type and is named after its tumor code. See details
}
  \item{tumor.freqs}{
A named vector of values between 0 and 1 which sum 1. It contains the expected proportion of patients that are planned 
to be recruited. See Details
}
  \item{plotFreq}{
If TRUE, the plot return the relative frequencies instead of the absolute number of samples.
}
  \item{noPlot}{
If TRUE, the plot is not shown but just the data used to drawn it.
}
  \item{html}{
If TRUE, an html interactive version of the plot is reported using \code{googleVis}.
}
}
\details{
This plot is a more compact (although less informative) version of the 
\code{\link{coveragePlot}}.
According to the chosen \code{alterionType}, the package will look 
for all the samples with available data for all the selected \code{alterationType}. 
For example, if alteratonType = c( "mutations" , "copynumber"), 
only common samples with both mutation and copynumber data are used. 
If both 'var' and 'grouping' are set, the plot will show two bars for every level of 'var'. 
The first one is a breakdown by 'grouping', while the second one is the total number 
of unique samples covered by at least one alteration. The first bar of the two is 
generally higher, because the breakdown does not sum up. 
For example, if we show a coverage stack plot of "drug" divided by "gene_symbol", 
the first bar will show the number of covered samples by every gene 
(considering a sample twice if is altered in more than one gene). 
The second bar is the total number of covered samples for the drug. The legend 
is not plotted if grouping is set to NA.

By default, \code{coverageStackPlot} will use all the available data from the object, using all the samples for the requested alterationTypes. Nevertheless, one could 
be interested in creating a compound design that is composed by a certain number of samples per tumor type. This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and trial can be stopped when the desired sample size for a given tumor type is reached. 
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
If noPlot=FALSE, this method returns a bar chart. Y-axis represents the number of samples, X-axis the number of alterations 
per sample. In case tumor.freqs is set, the Y-axis represents the relative frequency that is reported as text on the top of the bars.
If noPlot=TRUE, it returns a named list: 
\item{plottedTable }{a matrix with absolute number of samples plotted. 
Every column is a level of 'var' while every row 
represents one of the possible breakdown ('grouping').
}
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{saturationPlot}} \code{\link{coveragePlot}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Plot the number of covered samples 
# Using mutations and copynumber data
coverageStackPlot(cpObj , alterationType=c("mutations" , "copynumber") 
          , var="drug"
          , grouping="gene_symbol"
          , tumor_type="brca")
# Show an interactive version of the plot
# Save the html code first
myHtmlPlot <- coverageStackPlot(cpObj 
          , alterationType=c("mutations" , "copynumber") 
          , var="drug"
          , grouping="gene_symbol"
          , tumor_type="brca"
          , noPlot=FALSE
          , html=TRUE)
# Plot the code above
plot(myHtmlPlot)
}
