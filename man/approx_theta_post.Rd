% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{approx_theta_post}
\alias{approx_theta_post}
\title{Laplace approximation of the posterior of allelic imbalance theta}
\usage{
approx_theta_post(
  pAD,
  DP,
  p_s,
  lower = 0.001,
  upper = 0.499,
  start = 0.25,
  gamma = 20
)
}
\arguments{
\item{pAD}{numeric vector Variant allele depth}

\item{DP}{numeric vector Total allele depth}

\item{p_s}{numeric vector Variant allele frequency}

\item{lower}{numeric Lower bound of theta}

\item{upper}{numeric Upper bound of theta}

\item{start}{numeric Starting value of theta}

\item{gamma}{numeric Gamma parameter of the beta-binomial distribution}
}
\description{
Laplace approximation of the posterior of allelic imbalance theta
}
\keyword{internal}
