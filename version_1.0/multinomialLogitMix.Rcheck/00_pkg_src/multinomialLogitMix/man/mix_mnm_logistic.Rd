\name{mix_mnm_logistic}
\alias{mix_mnm_logistic}
\title{
EM algorithm
}
\description{
Estimation of the multinomial logit mixture using the EM algorithm. The algorithm exploits a careful initialization procedure (Papastamoulis et al., 2016) combined with a ridge-stabilized implementation of the Newton-Raphson
method (Goldfeld et al., 1966) in the M-step.
}
\usage{
mix_mnm_logistic(y, X, Kmax = 10, maxIter = 100, emthreshold = 1e-08, 
	maxNR = 5, nCores, tsplit = 8, msplit = 5, split = TRUE, 
	shake = TRUE, random = TRUE, criterion = "ICL", 
	plotting = FALSE, R0 = 0.1, method = 5)
}
\arguments{
  \item{y}{
matrix of counts
}
  \item{X}{
design matrix (including constant term).
}
  \item{Kmax}{
Maximum number of mixture components. 
}
  \item{maxIter}{
Maximum number of iterations. 
}
  \item{emthreshold}{
Minimum loglikelihood difference between successive iterations in order to terminate. 
}
  \item{maxNR}{
maximum number of Newton Raphson iterations
}
  \item{nCores}{
number of cores for parallel computations. 
}
  \item{tsplit}{
positive integer denoting the number of different runs for each call of the splitting small EM used by split-small EM initialization procedure.
}
  \item{msplit}{
positive integer denoting the number of different runs for each call of the splitting small EM.
}
  \item{split}{
Boolean indicating if the split initialization should be enabled in the small-EM scheme.
}
  \item{shake}{
Boolean indicating if the shake initialization should be enabled in the small-EM scheme.
}
  \item{random}{
Boolean indicating if random initializations should be enabled in the small-EM scheme.
}
  \item{criterion}{
set to "ICL" to select the number of clusters according to the ICL criterion. 
}
  \item{plotting}{
Boolean for displaying intermediate graphical output. 
}
  \item{R0}{
controls the step size of the update: smaller values result to larger step
sizes. See description in paper.
}
  \item{method}{
this should be set to 5. 
}
}
\value{
\item{estimated_K}{selected value of the number of clusters.}
\item{all_runs}{detailed output per run.}
\item{BIC_values}{values of bayesian information criterion.}
\item{ICL_BIC_values}{values of ICL-BIC.}
\item{estimated_clustering}{Single best-clustering of the data, according to the MAP rule.}
}
\references{
Papastamoulis P (2022). Model-based clustering of multinomial count data.  arXiv:2207.13984 [stat.ME]
}
\author{
Panagiotis Papastamoulis
}

\examples{
#	Generate synthetic data

	K <- 2
	p <- 2
	D <- 3
	n <- 2
	set.seed(116)
	simData <- simulate_multinomial_data(K = K, p = p, D = D, n = n, size = 20, prob = 0.025)   

	
	SplitShakeSmallEM <- mix_mnm_logistic(y = simData$count_data, 
		X = simData$design_matrix, Kmax = 2, maxIter = 1, 
		emthreshold = 1e-8, maxNR = 1, nCores = 2, tsplit = 1, 
		msplit = 2, split = TRUE, R0 = 0.1, method = 5, 
		plotting = FALSE)
	#selected number of clusters
	SplitShakeSmallEM$estimated_K
	#estimated single best-clustering, according to MAP rule
	SplitShakeSmallEM$estimated_clustering
	# detailed output for all parameters of the selected number of clusters
	SplitShakeSmallEM$all_runs[[SplitShakeSmallEM$estimated_K]]
	
}
