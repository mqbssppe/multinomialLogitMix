\name{mala_proposal}
\alias{mala_proposal}
\alias{mala_proposal_cpp}
\alias{log_mix_prior_derivative}
\alias{log_prior_mix}
\title{
Proposal mechanism of the MALA step.
}
\description{
Only the mala_proposal_cpp function is used in the package - which is written as an RCPP function. 
}
\usage{
mala_proposal(y, X, b, z, tau, A = FALSE, pr, nu2)
}
\arguments{
  \item{y}{
count data
}
  \item{X}{
design matrix
}
  \item{b}{
coefficients (array
}
  \item{z}{
allocation vector
}
  \item{tau}{
prior variance
}
  \item{A}{
A
}
  \item{pr}{
mixing proportions
}
  \item{nu2}{
parameter nu2
}
}
\value{
\item{theta}{theta values}
\item{b}{coeeficients}
\item{acceptance}{log-likelihood.}
\item{gradient}{log-likelihood.}
}
\author{
Panagiotis Papastamoulis
}


