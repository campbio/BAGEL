% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_tables.R
\name{create_sbs96_table}
\alias{create_sbs96_table}
\title{Uses a genome object to find context and generate standard SBS96 tables}
\usage{
create_sbs96_table(bay, g, overwrite = FALSE)
}
\arguments{
\item{bay}{Input samples}

\item{g}{A \linkS4class{BSgenome} object indicating which genome
reference the variants and their coordinates were derived from.}

\item{overwrite}{Overwrite existing count table}
}
\value{
Returns the created SBS96 count table object
}
\description{
Uses a genome object to find context and generate standard SBS96 tables
}
