% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_normExpr.TPM}
\alias{calc_normExpr.TPM}
\title{Calculate normalized expresssion for full-length RNA data}
\usage{
calc_normExpr.TPM(x, logarithmetic = T)
}
\arguments{
\item{x}{Matrix or sparse matrix. Each column specify a sample.}

\item{logarithmetic}{Logical. Should the results be log-transformed?}
}
\value{
Matrix
}
\description{
'calc_normExpr.TPM' is a wrapper to calculate the normalized expression for count data derived from full-length protocol.
It would calculate the normalized expresssion as TPM for each sample.
}
