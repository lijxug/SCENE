% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/context.R
\name{detectDiffContext}
\alias{detectDiffContext}
\title{Detect different co-varying pattern in contexts
This function performs multi-variable linear regression on context and programs, then select
program pairs that are significantly influenced by context}
\usage{
detectDiffContext(bulk_pg_dH_mt, context_num_tbl, signif_cutoff = 0.05)
}
\arguments{
\item{bulk_pg_dH_mt}{matrix, with columns as programs, rows as bulk samples, calculated as \code{d} \%*\% \code{H}, where \code{d} and \code{H} are from \code{cNMF}}

\item{context_num_tbl}{tibble, with colnames \code{label}, \code{context}, and \code{num}, indicating the sample name, context, and context number respectively}

\item{signif_cutoff}{Cutoff to determine cProgram pairs with significant context-factor}
}
\value{
A list contains the result of multi-linear regreassion
}
\description{
Detect different co-varying pattern in contexts
This function performs multi-variable linear regression on context and programs, then select
program pairs that are significantly influenced by context
}
