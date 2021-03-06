% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_class.R
\name{set_sample_annotations}
\alias{set_sample_annotations}
\title{Set sample level annotations for bagel object}
\usage{
set_sample_annotations(bay, annotations)
}
\arguments{
\item{bay}{Bagel object we input sample into}

\item{annotations}{Sample DataFrame}
}
\value{
Sets sample_annotations slot {no return}
}
\description{
Set sample level annotations for bagel object
}
\examples{
bay <- readRDS(system.file("testdata", "bagel_sbs96.rds", package = "BAGEL"))
sample_annotations <- read.table(system.file("testdata",
"sample_annotations.txt", package = "BAGEL"), sep = "\t", header=TRUE)
set_sample_annotations(bay, data.table::as.data.table(sample_annotations))
}
