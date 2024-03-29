% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{select_features}
\alias{select_features}
\title{Select features for bulk data}
\usage{
select_features(X, loess_span = 0.3, high_quantile = 0.5)
}
\arguments{
\item{X}{a library-size-normalized matrix or a UMI count matrix or a TPM matrix.}

\item{loess_span}{smooth factor for loess, default = 0.3}

\item{high_quantile}{quantile for selecting the variable genes,}
}
\value{
A list object with two named field: \emph{fit} for loess fit object; \emph{hvg_df} for a dataframe with highly variable genes.
}
\description{
Choose highly variable genes, using 'vst' algorithm adapted from Seurat.
}
