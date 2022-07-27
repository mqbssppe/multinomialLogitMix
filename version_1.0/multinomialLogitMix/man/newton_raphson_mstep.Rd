\name{newton_raphson_mstep}
\alias{newton_raphson_mstep}
\title{
M-step of the EM algorithm
}
\description{
Implements the maximization step of the EM algorithm based on a ridge-stabilized version of the Newton-Raphson algorithm, see Goldfeld et al. (1966). 
}
\usage{
newton_raphson_mstep(y, X, b, w, maxNR = 5, R0 = 0.1, method = 5, verbose = FALSE)
}
\arguments{
  \item{y}{
count data matrix
}
  \item{X}{
design matrix (including const).
}
  \item{b}{
coefficients of the multinomial logit mixture
}
  \item{w}{
mixing proportions
}
  \item{maxNR}{
threshold
}
  \item{R0}{
inital value for the parameter that controls the step-size of the update. 
}
  \item{method}{
set to 5. Always.
}
  \item{verbose}{
Boolean.
}
}
\value{
\item{b}{coefficients}
\item{theta}{theta values}
\item{ll}{log-likelihood.}
}
\references{
Goldfeld, S. M., Quandt, R. E., and Trotter, H. F. (1966). Maximization by quadratic hill-climbing. Econometrica: Journal of the Econometric Society, 541-551.
}
\author{
Panagiotis Papastamoulis
}
