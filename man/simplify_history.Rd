% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{simplify_history}
\alias{simplify_history}
\title{Simplify the mutational history based on likelihood evidence}
\usage{
simplify_history(G, l_matrix, max_cost = 150, n_cut = 0, verbose = TRUE)
}
\arguments{
\item{G}{igraph Mutation graph}

\item{l_matrix}{matrix Mutation placement likelihood matrix (node by mutation)}
}
\value{
igraph Mutation graph
}
\description{
Simplify the mutational history based on likelihood evidence
}
\keyword{internal}
