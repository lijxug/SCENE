% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/context.R
\name{calcCorrelative4}
\alias{calcCorrelative4}
\title{Identify correlative pattern.ver4}
\usage{
calcCorrelative4(
  program_mt,
  sources,
  celltypes,
  nonzero_ratio = 0.1,
  rho_cutoff = 0,
  n_permtest = 2000,
  n_cores = parallel::detectCores() - 1,
  qval_cutoff = 0.05,
  contexts = NULL,
  progs = NULL,
  .context = NULL,
  verbose = T
)
}
\arguments{
\item{program_mt}{Matrix, processed}

\item{sources}{Named strings, c\('sample_name' = 'bulk|singlecell'\)}

\item{celltypes}{Named vector, should be a cell-id-named vector of cell types}

\item{nonzero_ratio}{Value between 0-1, c\('sample_name' = 'bulk|singlecell'\)}

\item{rho_cutoff}{The remaining program pairs must have spearman correlation rho larger than this. Default = 0.}

\item{n_permtest}{Number of permutation rounds to perform, default = 2000}

\item{n_cores}{Number of cores to use for permutation, default = the maximum available cores - 1}

\item{qval_cutoff}{The remaining program pairs must have qvalue smaller than this. Default = 0.05. BH method.}

\item{contexts}{Named vector, specify the context information of each bulk sample. If not specified, all bulk samples are considered to be from the same conditions.}

\item{progs}{Character vector, names of the list that derived from \code{detectDiffContext}, indicating program pairs that are significantly influenced by the context. Must be provided when running in context mode.}

\item{.context}{Character, default NULL. Specify one of the contexts to run.}
}
\value{
A named list with the correlative information between celltypes.
}
\description{
Identify correlative pattern.ver4
}
