\name{saturationPlot}
\alias{saturationPlot}
\title{
Plot your panel along with incremental genomic space occupied adding one piece at a time
}
\description{
This plot method returns a scatter plot with genomic space on X axis and average/absolute number of alterations on Y axis. 
The way the plot is built is incremental. We add one feature at a time starting from the most altered and we see how many samples 
we include at each step and how much space is occupied. 
}
\usage{
saturationPlot(object
  , alterationType = c( "copynumber", "expression", "mutations", "fusions")
  , grouping = c(NA, "drug", "group", "alteration_id", "tumor_type")
  , adding = c( "alteration_id", "gene_symbol", "drug", "group")
  , tumor_type = NULL
  , y_measure = c( "mean", "absolute")
  , adding.order=c( "absolute", "rate")
	, sum.all.feature=FALSE
	, collapseMutationByGene=TRUE
	, collapseByGene=FALSE
  , labelling=TRUE
  , tumor.weights=NULL
  , main=""
  , legend=c("in" , "out")
  , noPlot = FALSE)
}
\arguments{
  \item{object}{
a CancerPanel object
}
  \item{alterationType}{
what kind of alteration to include. It can be one or more between 
"copynumber", "expression", "mutations", "fusions".
Default is to include all kind of alterations.
}
  \item{grouping}{
One of the following: "drug", "group", "alteration_id", "tumor_type". This parameter 
draws a curve for every level of the chosen grouping.
if set to NA, the panel is not split and the plot is a single curve.
}
  \item{adding}{
One of the following: "alteration_id", "gene_symbol", "drug", "group". This parameter 
will set which variable is added at every point of the plot. see details
}
  \item{tumor_type}{
only plot one or more tumor types among the ones available in the object.
}
  \item{y_measure}{
if 'mean', the measure on Y axis is the mean number of 
alterations per sample. Confidence interval of the measure 
is also reported. If 'absolute', the relative frequency 
of samples covered by at least one alteration is reported, 
similarly to \code{\link{coveragePlot}.}
}
  \item{adding.order}{ This parameter modifies the order of entrance of the adding variable. 
If 'absolute', the adding variable starts from the most altered up to the less frequently altered. If 'rate', 
the order of entrance, from left to right is based on the number of alterations divided by the length in kb. 
}
  \item{sum.all.feature}{
logical. if TRUE every gene length of the panel is summed up by the adding variable. The effect is that if 
a gene is considered both for SNV and CNA, it is counted twice.
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
  \item{labelling}{
if FALSE, the dots are not labelled. It is useful for very large comparative plots. Default TRUE.
}
  \item{tumor.weights}{
A named vector of integer values containing an amount of samples to be randomly sampled from the data. 
Each element should correspond to a different tumor type and is named after its tumor code. See details
}
  \item{main}{
Set a name for the plot
}
  \item{legend}{
if 'in' the legend is plotted in the top left corner, if 'out', outside of the plotting area 
}
  \item{noPlot}{
if TRUE, the plot is not shown and data to create it are reported instead.
}
}
\details{
This plot is particularly useful to evaluate the panel piece by piece. 
At the last point, we can observe the maximum coverage or maximum mean value of the panel, as reported in the coveragePlot. 
If we go back one point at a time, we can appreciate how many samples we gained by adding a new drug or a new gene to the panel. 
It is often the case that our panel is redundant for certain drugs or genes and there is no point in wasting sequencing space for 
a gene that is poorly altered and doesn't allow further improvement to our clinical trial.

By default, \code{saturationPlot} will use all the available data from the object, 
using all the samples for the requested alterationTypes. Nevertheless, one could 
be interested in creating a compound design that is composed by a certain number 
of samples per tumor type. This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and your design 
can be stopped when the desired sample size for a given tumor type is reached. 
By adding tumor.weights, we can achieve such target (see examples). 
Unfortunately, there are two main drawbacks in doing so:
\enumerate{
  \item{small sample size:}{ by selecting small random samples, the real frequency can be distorted. 
                            to avoid this, it is better to run several small samples and then bootstrap them}
  \item{recycling:}{ if the sample size for a tumor type requested by the user is above the available number of cBioportal samples, 
                  the samples are recycled. This has the effect of stabilizing the frequencies but y_measure = "absolute" will have no real meaning when 
                  the heterogeneity of the samples is lost.}
}


If noPlot is TRUE, the method returns a data.frame with 8 or 9 columns, 
depending on how the adding.order parameter was set: 
\describe{
  \item{gene_symbol , drug , alteration_id or group}{ the adding variable chosen by the user}
  \item{grouping}{ the grouping variable chosen by the user}
  \item{Mean}{ the value plotted on Y axis if 'mean' is chosen as y_measure parameter}
  \item{Coverage}{ the value plotted on Y axis if 'absolute' is chosen as y_measure parameter}
  \item{SD}{ standard deviation of Mean}
  \item{SE}{ standard error of Mean}
  \item{CI}{ confidence interval of Mean}
  \item{Space}{ genomic space in kBases ordered by grouping variable}
  \item{num_of_variants_per_KB}{ if adding.order='rate', 
    this additional column is added. It represents the number 
    of alteration divided by the feature length}
}
}
\value{
An incremental scatter plot if noPlot is FALSE, a data.frame otherwise.
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{coveragePlot}}
}
\examples{
# Load example CancerPanel object
data(cpObj)

# Plot the saturation of this panel by tumor type adding one drug at a time
# Using mutations and copynumber data
saturationPlot(cpObj , alterationType=c( "mutations" , "copynumber")
          , adding="drug"
          , grouping="tumor_type"
          , y_measure="absolute"
          )
# Report data used in the plot above
myplotData <- saturationPlot(cpObj , alterationType=c( "mutations" , "copynumber")
          , adding="drug"
          , grouping="tumor_type"
          , y_measure="absolute"
          , noPlot=TRUE)
str(myplotData)
# Plot the saturation of the panel by tumor type adding one gene at a time
# Using mutation and copynumber data
# Modify the order of entrance 
# using the rate between number of alterations and length
saturationPlot(cpObj , alterationType=c( "mutations" , "copynumber") 
              , adding="gene_symbol" 
              , grouping="tumor_type" 
              , adding.order="rate" 
              , y_measure="absolute"
              )
# Plot with no grouping giving more weight to lung cancer samples
# Note that we ask for more samples than the availables in luad dataset
# the code will recycle the samples to account for this forced disequilibrium
saturationPlot(cpObj , alterationType=c( "mutations" , "copynumber") 
              , adding="gene_symbol"
              , y_measure="mean"
              , tumor.weights=c(brca=500 , luad=2000)
              )
}

