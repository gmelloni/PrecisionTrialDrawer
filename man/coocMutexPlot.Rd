\name{coocMutexPlot}
\alias{coocMutexPlot}
\title{
Plot cooccurence and mutual exclusivity between pairs of drugs or genes
}
\description{
This plot reports how two features, like drugs or genes, 
can be considered close or distant in terms of 
occurrence on the same set of samples. In case of genes, 
for example, it consider the number of times a gene 
is altered together with another one or in a mutual exclusive fashion. 
In case of drugs, all the drug target are pulled 
together and you can appreciate if they act on the same of different targets.
}
\usage{
coocMutexPlot(object
 , var=c("drug","group","gene_symbol")
 , alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
 , grouping=c(NA , "drug" , "group" , "alteration_id" , "tumor_type")
 , tumor_type=NULL
 , collapseMutationByGene=FALSE
 , collapseByGene=FALSE
 , tumor.weights=NULL
 , style=c("cooc" , "dendro")
 , prob = c("hyper","firth")
 , drop=FALSE
 , noPlot=FALSE
 , pvalthr=0.05
 , plotrandom=TRUE
 , ncolPlot=FALSE
 , ...)
}
\arguments{
  \item{object}{
a CancerPanel object
}
  \item{var}{
One of the following: "drug", "group", "gene_symbol". This parameter 
will set which variable is used in the plot
}
  \item{alterationType}{
what kind of alteration to include. It can be one or more between 
"copynumber", "expression", "mutations", "fusions".
Default is to include all kind of alterations.
}
  \item{grouping}{
One of the following: "drug", "group", "alteration_id", "tumor_type". This parameter 
draws a plot for every level of the chosen grouping.
if set to NA, the panel is not split and the plot is one.
}
  \item{tumor_type}{
A character vector of tumor types to include in the plot among the one included in the object
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
  \item{style}{
If 'cooc', the default, it performs pairwise cooccurence and mutual exclusivity test and the plot 
is a pvalue upper triangle heatmap. If dendro, it performs hierarchical clustering using binary distance 
between 'var' subjects.
}
  \item{prob}{
One of the following: "hyper" or "firth". Two ways of calculating cooccurrence mutex 
pvalues. The first uses the hypergeometric distribution as in the Fisher test. The second uses a penalized 
logistic regression and is particularly indicated when the alterations are rare.
}
  \item{drop}{
Logical indicating if the table of cooccurence should include (FALSE) 
or not include (TRUE) the samples that are never altered 
for any element of the var of interested. Default FALSE. See details.
}
  \item{noPlot}{
if TRUE, the plot is not shown and data to create it are reported.
}
  \item{pvalthr}{
In the plot, every square under the threshold is depicted in gray
}
  \item{plotrandom}{
If TRUE, all elements of var are reported, even if they have no significant pair
}
  \item{ncolPlot}{
Number of columns required for a multiplot. With default = FALSE, the function calculates the optimal configuration based on the number of plots that need to be printed
}
 \item{...}{
Further arguments passed to \code{hclust} function
}
}
\details{
This plot explores if there is cooccurrence or mutual 
exclusivity between features selected in the panel. 
It is particularly useful to evaluate the opportunity 
to add a new gene or a new drug in an umbrella design. 
Two drugs that acts on mutual exclusive pathways are more 
suitable for an umbrella design that seek at 
enlarging the spectrum of covered samples, even though 
one of the two drugs has few affected samples. 
On the other hand, if a drug has been proven to be more 
effective or reliable and its target are alterated 
together with another drug, there is no point in adding 
the less effective cure. Another way of seeing this 
feature is by using clustering adding option 'dendro' to 
style parameter. The reported plot will lack of pvalues 
but it is more general including distances between 
hierarchically aggregated drugs or genes.

If noPlot is TRUE, the method returns a data.frame with 6 columns in case of style 'cooc':
\describe{
  \item{sp1_name}{ the first 'var' value for the mut-ex analysis}
  \item{sp2_name}{ the second 'var' value for the mut-ex analysis}
  \item{pVal.MutEx}{ pvalue associated with mutual exclusivity evaluation}
  \item{pVal.Cooc}{ pvalue associated with cooccurence evaluation}
  \item{OR}{ corrected odds ratio, of the confusion matrix between sp1 and sp2}
  \item{grouping}{ grouping variable chosen by the user}
}
If noPlot is TRUE, the method returns a list of hclust objects in case of style 'dendro'

The option drop can completely change the results, so check exactly what is your initial question. 
Let's imagine a coocMutex plot by gene. If drop=FALSE, all the samples tested for mutations 
on those genes are included. Otherwise, only the samples with at least one mutation 
in the genes of interested will be included. The default is to keep all the samples but this procedure is biased towards 
cooccurrence. This is caused by the fact that mutations are rare, 
so the cooccurrence of no mutations is generally very high and it counts as no-mutations.

By default, \code{coocMutexPlot} will use all the available data from the object, 
using all the samples for the requested alterationTypes. Nevertheless, one could 
be interested in creating a compound design that is composed by a certain number 
of samples per tumor type. This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and your trial 
can be stopped when the desired sample size for a given tumor type is reached. 
By adding tumor.weights, we can achieve such target.
Unfortunately, there are two main drawbacks in doing so:
\enumerate{
  \item{small sample size:}{ by selecting small random samples, the real frequency can be distorted. 
                            to avoid this, it is better to run several small samples and then aggregate the results}
  \item{recycling:}{ if the sample size for a tumor type requested by the user is above the available number of cBioportal samples, 
                  the samples are recycled. This has the effect of stabilizing the frequencies but y_measure = "absolute" will have no real meaning when 
                  the heterogeneity of the samples is lost.}
}
}
\value{
In case of style 'cooc', an upper triangle discrete heatmap if noPlot is FALSE, a data.frame otherwise. 
In case of style 'dendro', a dendrogram if noPlot is FALSE, a list of patient by alteration matrices.
}
\references{
\href{https://CRAN.R-project.org/package=cooccur}{ code re-written from package cooccur}
\href{https://CRAN.R-project.org/package=brglm}{ implementation of penalized glm, here logistic regression}
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{saturationPlot}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Plot cooccurence and mutual exclusivity between pairs of genes by tumor type
coocMutexPlot(cpObj
              , var="gene_symbol" 
              , grouping="tumor_type" 
              , alterationType="mutations")
}

