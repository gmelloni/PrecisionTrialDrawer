\name{panelOptimizer}
\alias{panelOptimizer}
\title{
A function that open a shiny app to optimize panel size 
while retaining the majority of mutations
}
\description{
Given a CancerPanel object, it automatically reads the genes 
and panel information and allows a custom 
subsetting of the panel to retain the largest amount 
of mutations while saving genomic space
}
\usage{
panelOptimizer(object)
}
\arguments{
  \item{object}{
A CancerPanel object filled with genomic data.
}
}
\details{
This function reads the panel objects and retrieves 
all the genes requested for SNV and full sequence. 
It also collects all the mutation data and the tumor types 
available and performs an analysis 
on the position of the mutations on 
the canonical protein sequence for every gene. 

At function call, a shiny app with four tabs is opened. 
Select a gene and the tumor types to use and click Run.

On the first and second tabs, the user is guided in the 
choice of most appropriate regions using an in-house Bioconductor 
package called \code{LowMACA}. \code{LowMACA} creates a null model 
where all the mutations on each gene are 
randomly permutated along the sequence. Every position that exceeds 
the threshold of 95\% confidence 
interval is considered not random and represents an hotspot. 
The user can also decide to use a bandwidth in 
this calculation and apply a Gaussian density to the distribution 
of mutations along the sequence. 
Alterations that are closed to each other in the sequence will be aggregated 
to form significant regions.
If the LowMACA analysis succeded, a table will appears under the plot 
in the first tab and on the second tab a list 
of significant positions is also shown. 
Click on Store LowMACA yellow button and all the regions 
identified by the algorithm will be stored in the fourth tab 
(Optimize Panel).

If you want to select your own custom regions, 
go to the third tab (Manual Selection). You can click on the red dots and 
retrieve information on specific mutations. 
If you drag a region with the mouse, a table will appear below, 
with the same information as the LowMACA analysis. 
Click on Store yellow button to keep 
the region selected and move it to the fourth tab.

For a new analysis on a different gene, 
just select the gene and click Run again.

When satisfied, just click on Close and save or simply close the browser page. 
All the regions selected will be merged and 
returned in standard output.
}
\value{
A list of three elements after closing the shiny session.
\enumerate{
  \item{regions}{ a data.frame with the regions selected, 
               the percentage of space occupied and percentage of mutations captured}
  \item{mergedRegions}{ merged protein regions by gene from user section}
  \item{panel}{ if no regions were selected, it reports the original panel, 
               otherwise it substitutes the gene requested 
               in full sequence with the regions in mergedRegions }
}
}
\author{
Giorgio Melloni, Alessandro Guida
}
\seealso{
\code{\link[LowMACA]{entropy}} \code{\link[LowMACA]{lmPlot}}
}
\examples{
## Only run this example in interactive R sessions
if (interactive()) {
  # Load example CancerPanel object
  data(cpObj)
  # Optimize the space on the shiny app. 
  # All changes mad on the app will be saved 
  newpanel <- panelOptimizer(cpObj)
  # If some changes have been made, recreate a new CancerPanel object
  if(!is.null(newpanel$regions)){
    cpObjOptimized <- newCancerPanel(newpanel$panel)
    # Fill the object with the same data of the non optimized panel
    cpObjOptimized <- getAlterations(cpObjOptimized , repos=cpData(cpObj))
    # Subset alterations on the new panel directives
    cpObjOptimized <- subsetAlterations(cpObjOptimized)
  }
}
}
