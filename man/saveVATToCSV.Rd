% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dealing.R
\name{saveVATToCSV}
\alias{saveVATToCSV}
\title{Save VAT data or analysis result}
\usage{
saveVATToCSV(vat, filename, analysis.key = NULL, use.genes = NULL,
  use.cells = NULL, ...)
}
\arguments{
\item{vat}{VAT entity}

\item{filename}{CSV filename}

\item{analysis.key}{analysis result's key, default NULL, saving vat@data}

\item{use.genes}{which genes will be saved}

\item{use.cells}{which cells will be saved}
}
\description{
Save VAT data or analysis result
}
