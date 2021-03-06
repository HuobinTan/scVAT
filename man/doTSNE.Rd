% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{doTSNE}
\alias{doTSNE}
\title{do tSNE using Rtsne library}
\usage{
doTSNE(vat, dims = 2, analysis.key = "PCA", use.col = 50,
  key = "tSNE", seed = NULL, ...)
}
\arguments{
\item{vat}{vat entity}

\item{dims}{tSNE dimension}

\item{analysis.key}{tSNE input data key, default "PCA", if NULL, use vat@data}

\item{use.col}{columns number for tSNE}

\item{key}{set tSNE result's key, default "tSNE"}

\item{seed}{random seed for tSNE, default NULL}

\item{...}{see also Rtsne::Rtsne parameters, For example, set check_duplicates = FALSE for Error "Remove duplicates before running TSNE."}
}
\description{
do tSNE using Rtsne library
}
\examples{
vat <- doTSNE(vat, dims = 2, use.col = 50, seed = 100)

}
