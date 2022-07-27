\name{gibbs_mala_sampler}
\alias{gibbs_mala_sampler}
\title{
The core of the Hybrid Gibbs/MALA MCMC sampler for the multinomial logit mixture.
}
\description{
This function implements Gibbs sampling to update the mixing proportions and latent allocations variables of the mixture model. The coefficients of the logit model are updated according to Metropolis-Hastings type move, based on a Metropolis adjusted Langevin (MALA) proposal. 
}
\usage{
gibbs_mala_sampler(y, X, tau = 3e-05, nu2, K, mcmc_iter = 100, 
	alpha_prior = NULL, start_values = "EM", em_iter = 10, 
	thin = 10, verbose = FALSE, checkAR = NULL, 
	probsSave = FALSE, ar_low = 0.4, ar_up = 0.6)
}
\arguments{
  \item{y}{
matrix of counts.
}
  \item{X}{
design matrix (including constant term). 
}
  \item{tau}{
the variance of the normal prior distribution of the logit coefficients.
}
  \item{nu2}{
scale of the MALA proposal (positive).
}
  \item{K}{
number of components of the (overfitting) mixture model.
}
  \item{mcmc_iter}{
Number of MCMC iterations.
}
  \item{alpha_prior}{
Parameter of the Dirichlet prior distribution for the mixing proportions.
}
  \item{start_values}{
Optional list of starting values. Random initialization is used if this is not provided.
}
  \item{em_iter}{
Maximum number of iterations if an EM initialization is enabled.
}
  \item{thin}{
optional thinning of the generated MCMC output.
}
  \item{verbose}{
Boolean.
}
  \item{checkAR}{
Number of iterations to adjust the scale of the proposal in MALA mechanism during the initial warm-up phase of the sampler.
}
  \item{probsSave}{
Optional.
}
  \item{ar_low}{
Lowest threshold for the acceptance rate of the MALA proposal (optional) .
}
  \item{ar_up}{
Highest threshold for the acceptance rate of the MALA proposal (optional). 
}
}
\value{
\item{nClusters}{sampled values of the number of clusters (non-empty mixture components).}
\item{allocations}{sampled values of the latent allocation variables.}
\item{logLikelihood}{Log-likelihood values per MCMC iteration.}
\item{mixing_proportions}{sampled values of mixing proportions.}
\item{coefficients}{sapled values of the coefficients of the multinomial logit.}
\item{complete_logLikelihood}{Complete log-likelihood values per MCMC iteration.}
\item{class_probs}{Classification probabilities per iteration (optional).}
\item{AR}{Acceptance rate of the MALA proposal.}
}
\author{
Panagiotis Papastamoulis
}
\note{
This function is used inside the prior tempering scheme, which is the main function. 
}
\seealso{
\code{\link{gibbs_mala_sampler_ppt}}
}
\examples{
#	Generate synthetic data
	K <- 2
	p <- 2
	D <- 2
	n <- 2
	set.seed(116)
	simData <- simulate_multinomial_data(K = K, p = p, D = D, n = n, size = 20, prob = 0.025)   


	gs <- gibbs_mala_sampler(y = simData$count_data, X = simData$design_matrix, 
		tau = 0.00035, nu2 = 100, K = 2, mcmc_iter = 3, 
		alpha_prior = rep(1,K), start_values = "RANDOM", 
		thin = 1, verbose = FALSE, checkAR = 100)

}
