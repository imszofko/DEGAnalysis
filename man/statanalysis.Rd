% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statanalysis_functions.R
\name{statanalysis}
\alias{statanalysis}
\title{Statistical Analysis}
\usage{
statanalysis(filterCounts, sampleTable)
}
\arguments{
\item{filterCounts, sampleTable}{Variable names for the filter count data and sample table}
}
\value{
dge and desMatrix
}
\description{
Runs Statistical Analysis on the count data
}
\examples{
stats <- statanalysis(countTable, sampleTable);
}
