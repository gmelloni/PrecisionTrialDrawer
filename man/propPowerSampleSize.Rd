\name{propPowerSampleSize}
\alias{propPowerSampleSize}
\title{
Calculate sample size or power required in a 
2-sample or 1-sample proportion equality study
}
\description{
This plot method returns a scatter plot of required sample size at 
screening by statistical power divided by 
couples of case-control probabilities.
}
\usage{
propPowerSampleSize(object
, var=c(NA , "drug" , "group" , "gene_symbol" , "alteration_id" , "tumor_type")
, alterationType=c("copynumber" , "expression" , "mutations" , "fusions")
, tumor_type=NULL
, stratum=NULL
, tumor.weights=NULL
, tumor.freqs=NULL
, pCase=NULL
, pControl=NULL
, side = c(2,1)
, type=c("chisquare" , "arcsin" , "exact")
, alpha=0.05
, power=NULL
, sample.size=NULL
, case.fraction=0.5
, collapseMutationByGene=TRUE
, collapseByGene=FALSE
, round.result=TRUE
, priority.trial=NULL
, priority.trial.order=c("optimal" , "as.is")
, priority.trial.verbose=TRUE
, noPlot=FALSE)
}
\arguments{
  \item{object}{
a CancerPanel object
}
  \item{var}{
one among NA , "drug" , "group" , "gene_symbol" 
, "alteration_id" or "tumor_type". It defines the arms 
of the studies to be projected. With var=NA, 
the projection of the entire panel is displayed.
}
  \item{alterationType}{
what kind of alteration to include. It can be one or more between 
"copynumber", "expression", "mutations", "fusions".
Default is to include all kind of alterations.
}
  \item{tumor_type}{
only plot one or more tumor types among the ones available in the object.
}
  \item{stratum}{
a character vector containing one or more specific 
elements of var to be plotted instead of all the arms of the study. 
If it is not present, a warning is raised and the full design is returned.
}
  \item{tumor.weights}{
A named vector of integer values containing an 
amount of samples to be randomly sampled from the data. 
Each element should correspond to a different tumor type 
and is named after its tumor code. See details
}
  \item{tumor.freqs}{
A named vector of values between 0 and 1 which sum 1. 
It contains the expected proportion of patients that are planned 
to be screened. See Details
}
  \item{pCase}{
a numerical vector of one or more postulated 
proportions for cases (between 0 and 1)
}
  \item{pControl}{
a numerical vector of one or more postulated 
proportions for controls of the same length of pCase (between 0 and 1)
}
  \item{side}{
perform a 2-tail or 1-tail calculation. Default 2
}
  \item{type}{
calculate sample size using chisquare, 
arcsin or exact method. chisquare is used for a 2-sample 
equality while the other two are used in case in case pControl is 
fixed (1-sample case). case.fraction has no effect using arcsin or exact type.
}
  \item{alpha}{
a numerical value between 0 and 1 that reports 
the type I error threshold. Default 0.05 (5\%)
}
  \item{power}{
a numerical vector of values between 0 and 1 
that expresses the level of type II error. It is 
used to estimate sample size
}
  \item{sample.size}{
a positive integer numerical vector that reports 
the postulated sample size at screening. It is used 
to estimate the power of the study.
}
  \item{case.fraction}{
a numerical value between 0 and 1 representing the fraction 
of total sample size allocated to group 'Case'. Group control have 
an allocation of 1 - case.fraction
}
  \item{collapseMutationByGene}{
A logical that collapse all mutations on the same gene for a single patient as 
a single alteration.
}
  \item{collapseByGene}{
A logical that collapse all alterations on 
the same gene for a single patient as 
a single alteration. e.g. if a sample has 
TP53 both mutated and deleted as copynumber, it 
will count for one alteration only.
}
  \item{round.result}{
logical indicating if the sample size should be rounded with ceiling or not.
}
  \item{priority.trial}{
A character vector of drugs or group levels to start 
the design of a priority trial. See Details.
}
  \item{priority.trial.order}{
Either "optimal" or "as.is". If "optimal" is used, 
the screening starts from the rarest drug or group level 
up to the most common to guarantee minimal sample size at screening. 
In case of "as.is", the order of priority.trial 
remains unchanged.
}
  \item{priority.trial.verbose}{
If TRUE, the result of a priority.trial will be a complete report in a 5-element list.
}
  \item{noPlot}{
if TRUE, the plot is not shown and data are reported instead.
}
}
\details{
This method estimates sample size or power on the basis 
of one of the two information. Using multiple sample 
sizes or power, power curves are reported simulating 
different scenarios. Power or sample size 
are required but not both at the same time. 
HR must be also set but if a vector is provided, the plot will 
show multiple curves according to the various hazard ratios. 
'p.event', 'alpha' and 'case.fraction' are 
instead fixed for all the arms of the study (represented by the 'var' parameter).

If noPlot=TRUE, a data.frame with 6 column is reported instead:
\describe{
    \item{Var}{ levels of chosen variable}
    \item{ScreeningSampleSize}{ total sample size estimation at 
        screening on the basis of frequency of alteration}
    \item{EligibleSampleSize}{ sample size estimated as sum 
        of cases and controls after screening}
    \item{Beta}{ tested beta values}
    \item{Power}{ tested 1 - beta values}
    \item{Proportion.In.Case.Control}{ couples of pCase - pControl tested}
}

The algorithm estimates sample size on the basis of 
no a priori probability of finding 
a case or control subject ("EligibleSampleSize" column). 
In a basket or umbrella design, 
this number must be multiplied by the frequency of 
alteration that we expect to find based on the simulation 
run on the panel. If our panel can cover the 
50\% of the samples with a target therapy 
and 100 samples are required to reach 80\% power, 
we have to screen at least 200 patients 
in order to reach the desired number of cases 
in the sample size ("ScreeningSampleSize" column).

Similarly, if you want to estimate the power of the panel 
given an estimated sample size, we first 
multiply 'sample.size' by the frequency of expected alterations 
and then perform power estimation. 
'sample.size' is therefore intended at screening.

When 'var' variable is set, the algorithm provides the estimated 
sample size for each stratum of the variable. 
For example, if we set it to 'drug', a power curve for each 
drug type is displayed, without taking into account possible 
overlaps. If a sample shows multiple targettable alterations, 
it will be reused for every drug type that targets 
those alterations.

By default, \code{propPowerSampleSize} will use all the available 
data from the object, using all the samples for the requested alterationTypes. 
Nevertheless, one could 
be interested in creating a compound design that is composed 
by a certain number of samples per tumor type. 
This is the typical situation of basket trials, where you seek for 
specific alteration, rather than specific tumor types and your design 
can be stopped when the desired sample size for a given tumor type is reached. 
By adding tumor.weights, we can achieve such target (see examples). 
Unfortunately, there are two main drawbacks in doing so:
\enumerate{
    \item{small sample size:}{ by selecting small random samples, 
        the real frequency can be distorted. 
        to avoid this, it is better to run several small 
        samples and then bootstrap them}
    \item{recycling:}{ if the sample size for a tumor type 
        requested by the user is above the available 
        number of cBioportal samples, the samples are recycled. 
        This has the effect of stabilizing the frequencies 
        but y_measure = "absolute" will have no real meaning when 
        the heterogeneity of the samples is lost.}
}

A user balanced design can be also obtained using 
\code{tumor.freqs} parameter. In this case the fraction of altered samples are 
first calculated tumor-wise and then reaggregated using 
the weights provided by \code{tumor.freqs}. 
If the fraction of altered samples are 0.3 and 0.4 for 
breast cancer and lung cancer respectively, 
if you set tumor.freqs = c(brca=0.9 , luad=0.1), 
the full design will have a frequency equal to 0.3*0.9 + 0.4*0.1 = 0.31, 
that is basically equal to the one of breast samples. 
If this parameter is not set, the total amount of 
samples available is used with unpredictable balancing.

Both tumor.freqs and tumor.weights can achieve a balanced design 
according to user specification. 
To have a quick idea of the sample size required, 
it is better to use the former. To get an idea about the possible 
distribution of sample size giving a few samples 
(for example a minimum and a maximum sample size) 
it is better to run the function with 
tumor.weights several times and aggregate the results.

If priority.trial is set, a cascade design is build up. 
Given a set of parameter (power, pCase, pControl , alpha, etc.) 
an Eligible Sample Size (ESS) 
is calculated that is the same across drugs/groups. 
The total Screening Sample Size (SSS) is calculated following this scheme:
\enumerate{
    \item{ Start screening with the first drug/group, 
        reaching the sample size necessary to reach ESS}
    \item{ From the samples not eligible for the first drug/group, 
        test the second drug/group and collects as many 
        samples as possible up to ESS}
    \item{ Continue using the samples not eligible to the 
        end of all drugs/levels. Stop if there are no leftovers.}
    \item{ If all the drugs/groups have reached ESS, stop. 
        Otherwise start a new screening with the first 
        drug/group that has not reached ESS}
    \item{ Repeat from point 2 up to completion}
}
If priority.trial.order is set, the user can decide if 
the drugs/group levels must follow a precise order (as.is) or if the screening 
can start from the rarest drug/group level up to the most common (optimal). 
Following the optimal priority trial guarantees the best possible allocation 
with the minimum screening.
}
\value{
If noPlot = FALSE (default) a scatter plot is returned. 
If noPlot = TRUE, a data.frame is returned. 
In case priority.trial is set, a list of length-5-lists is reported. 
See vignette for details.
}
\references{
Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 
2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 85/89
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link{coveragePlot}} \code{\link{survPowerSampleSize}}
}
\examples{
# Load example CancerPanel object
data(cpObj)
# Show the study design by tumor type:
# 3 pCase - pControl couples and 4 power levels
# The full design is weighted using tumor.freqs
# The final sample size is composed by 90\% luad and 10\% brca
propPowerSampleSize(cpObj 
  , var = "tumor_type"
  , pCase = c(0.7, 0.8 , 0.9)
  , pControl = rep(0.5 , 3) 
  , power=c(0.5 , 0.6 , 0.7 , 0.8 , 0.9)
  , tumor.freqs = c(brca=0.1 , luad=0.9))
# Return power levels giving sample sizes at screening
propPowerSampleSize(cpObj
  , var = NA 
  , pCase = c(0.7, 0.8 , 0.9)
  , pControl = rep(0.5 , 3) 
  , sample.size=c(100 , 300 , 500 , 1000) 
  , noPlot=FALSE)
}

