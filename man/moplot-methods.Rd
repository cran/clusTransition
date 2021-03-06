% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\docType{methods}
\name{moplot}
\alias{moplot}
\alias{moplot,Monic-method}
\alias{moplot-method}
\title{plot Method for output}
\usage{
moplot(object)

\S4method{moplot}{Monic}(object)
}
\arguments{
\item{object}{An object of class Monic}
}
\description{
This method plot 3 barplot and 1 line graph. The first stack barplot shows
SurvivalRatio and AbsorptionRatio, second barplot shows number of newly emerged clusters
at each time stamp, third barplot shows number of disapeared clusters at each time stamp. The line
graph shows passforward Ratio and Survival Ratio.
}
