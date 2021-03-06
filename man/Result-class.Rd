% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_class.R
\docType{class}
\name{Result-class}
\alias{Result-class}
\title{Object containing deconvolved/predicted signatures, sample weights, and
the bagel object the result was generated from}
\description{
Object containing deconvolved/predicted signatures, sample weights, and
the bagel object the result was generated from
}
\section{Slots}{

\describe{
\item{\code{signatures}}{A matrix of signatures by mutational motifs}

\item{\code{exposures}}{A matrix of samples by signature weights}

\item{\code{tables}}{A character vector of table names used to make the result}

\item{\code{type}}{Describes how the signatures/weights were generated}

\item{\code{bagel}}{The bagel object the results were generated from}

\item{\code{log_lik}}{Posterior likelihood of the result (LDA only)}

\item{\code{perplexity}}{Metric of goodness of model fit}

\item{\code{umap}}{List of umap data.frames for plotting and analysis}
}}

