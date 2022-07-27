\name{mixLoglikelihood_GLM}
\alias{mixLoglikelihood_GLM}
\title{
Log-likelihood of the multinomial logit. 
}
\description{
Log-likelihood of the multinomial logit. 
}
\usage{
mixLoglikelihood_GLM(y, theta, pi)
}
\arguments{
  \item{y}{
matrix of counts
}
  \item{theta}{
a three-dimensional array containing the multinomial probabilities per cluster, for each observation. 
}
  \item{pi}{
 a numeric vector of length K (the number of mixture components) containing the mixing proportions.
}
}
\value{
Log-likelihood value.
}
\author{
Panagiotis Papastamoulis
}

