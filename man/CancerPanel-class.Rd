\name{CancerPanel-class}
\docType{class}
\alias{CancerPanel-class}
\alias{show,CancerPanel-method}
\alias{getAlterations,CancerPanel-method}
\alias{subsetAlterations,CancerPanel-method}
\alias{panelDesigner,CancerPanel-method}
\alias{cpArguments,CancerPanel-method}
\alias{cpData,CancerPanel-method}
\alias{cpDataSubset,CancerPanel-method}
\alias{coveragePlot,CancerPanel-method}
\alias{coverageStackPlot,CancerPanel-method}
\alias{saturationPlot,CancerPanel-method}
\alias{coocMutexPlot,CancerPanel-method}
\alias{survPowerSampleSize,CancerPanel-method}
\alias{survPowerSampleSize1Arm,CancerPanel-method}
\alias{cpFreq,CancerPanel-method}
\alias{panelOptimizer,CancerPanel-method}
\alias{appendRepo,CancerPanel-method}
\alias{dataExtractor,CancerPanel-method}
\alias{propPowerSampleSize,CancerPanel-method}
\alias{filterMutations,CancerPanel-method}
\alias{filterFusions,CancerPanel-method}
\title{Class \code{CancerPanel}}
\value{
Class object
}
\description{
Class \code{CancerPanel} defines an object (S4 class) that contains data and 
custom arguments to run panel simulations and panel designs.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{newCancerPanel(panel, padding_length , utr)}.
}
\section{Constructor}{
newCancerPanel(panel=data.frame, padding_length=numeric , utr=logical)
}
\section{Slots}{

\describe{
\item{\code{arguments}}{Object of class \code{list} containing the all panel related information. 
It is filled during the \code{CancerPanel} object construction.}

\item{\code{dataFull}}{Object of class \code{list} containing all the alterations 
for all the genes in the panel and every tumor type requested. 
It is filled by the \code{getAlterations} method.}

\item{\code{dataSubset}}{Object of class \code{list} containing all the alteration requested. 
It is filled by the \code{subsetAlterations} method. 
Every alteration type is defined by a slot; all slots have the same data.frame layout.}
}
}
\section{Methods}{
  \describe{
    \item{show}{\code{show(object = "CancerPanel")}: ... }
    \item{getAlterations}{\code{getAlterations(object = "CancerPanel")}: ... }
    \item{subsetAlterations}{\code{subsetAlterations(object = "CancerPanel")}: ... }
    \item{panelDesigner}{\code{panelDesigner(object = "CancerPanel")}: ... }
    \item{cpArguments}{\code{cpArguments(object = "CancerPanel")}: ... }
    \item{cpData}{\code{cpData(object = "CancerPanel")}: ... }
    \item{cpDataSubset}{\code{cpDataSubset(object = "CancerPanel")}: ... }
    \item{coveragePlot}{\code{coveragePlot(object = "CancerPanel")}: ... }
    \item{coverageStackPlot}{\code{coverageStackPlot(object = "CancerPanel")}: ... }
    \item{saturationPlot}{\code{saturationPlot(object = "CancerPanel")}: ... }
    \item{coocMutexPlot}{\code{coocMutexPlot(object = "CancerPanel")}: ... }
    \item{survPowerSampleSize}{\code{survPowerSampleSize(object = "CancerPanel")}: ... }
    \item{survPowerSampleSize1Arm}{\code{survPowerSampleSize(object = "CancerPanel")}: ... }
    \item{cpFreq}{\code{cpFreq(object = "CancerPanel")}: ... }
    \item{panelOptimizer}{\code{panelOptimizer(object = "CancerPanel")}: ... }
    \item{appendRepo}{\code{appendRepo(object = "CancerPanel")}: ... }
    \item{dataExtractor}{\code{dataExtractor(object = "CancerPanel")}: ... }
    \item{propPowerSampleSize}{\code{propPowerSampleSize(object = "CancerPanel")}: ... }
    \item{filterMutations}{\code{filterMutations(object = "CancerPanel")}: ... }
    \item{filterFusions}{\code{filterFusions(object = "CancerPanel")}: ... }
	}
}
\examples{
showClass("CancerPanel")
}
\seealso{
\code{\link{newCancerPanel}}
}
\keyword{classes}

