\name{gibbs_mala_sampler_ppt}
\alias{gibbs_mala_sampler_ppt}
\title{
Prior parallel tempering scheme of hybrid Gibbs/MALA MCMC samplers for the multinomial logit mixture.
}
\description{
The main MCMC scheme of the package. Multiple chains are run in parallel and swaps between are proposed. Each chain uses different parameters on the Dirichlet prior of the mixing proportion. The smaller concentration parameter should correspond to the first chain, which is the one that used for inference. Subsequent chains should have larger values of concentration parameter for the Dirichlet prior. 
}
\usage{
gibbs_mala_sampler_ppt(y, X, tau = 3e-05, nu2, K, 
	mcmc_cycles = 100, iter_per_cycle = 10, dirPriorAlphas, 
	start_values = "EM", em_iter = 10, nChains = 4, nCores = 4, 
	warm_up = 100, checkAR = 50, probsSave = FALSE, 
	showGraph = 50, ar_low = 0.4, ar_up = 0.6, withRandom = TRUE)
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
  \item{mcmc_cycles}{
Number of MCMC cycles. At the end of each cycle, a swap between chains is attempted.
}
  \item{iter_per_cycle}{
Number of iterations per cycle. 
}
  \item{dirPriorAlphas}{
Vector of concentration parameters for the Dirichlet priors in increasing order.
}
  \item{start_values}{
Optional list of start values. Randomly generated values are used if this is not provided. 
}
  \item{em_iter}{
Maximum number of iterations if an EM initialization is enabled.
}
  \item{nChains}{
Total number of parallel chains. 
}
  \item{nCores}{
Total number of CPU cores for parallel processing of the \code{nChains}.
}
  \item{warm_up}{
Initial warm-up period of the sampler, in order to adaptively tune the scale of the MALA proposal (optional). 
}
  \item{checkAR}{
Number of iterations to adjust the scale of the proposal in MALA mechanism during the initial warm-up phase of the sampler.
}
  \item{probsSave}{
Optional.
}
  \item{showGraph}{
Optional.
}
  \item{ar_low}{
Lowest threshold for the acceptance rate of the MALA proposal (optional) .
}
  \item{ar_up}{
Highest threshold for the acceptance rate of the MALA proposal (optional). 
}

  \item{withRandom}{
Logical. If TRUE (default) then a random permutation is applied to the supplied starting values. 
}



}
\details{
See the paper for details.
}
\note{The output of the MCMC sampler is not identifiable, due to possible label switching. In order to draw meaningful inferences, the output should be post-processed by \code{\link{dealWithLabelSwitching}}.}
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
\references{
Papastamoulis, P (2022). Model-based clustering of multinomial count data. 
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



# apply mcmc sampler based on random starting values 

Kmax = 2
nChains = 2
dirPriorAlphas  = c(1, 1 + 5*exp((seq(2, 14, length = nChains - 1)))/100)/(200)
nCores <- 2
mcmc_cycles <- 2
iter_per_cycle = 2
warm_up <- 2

mcmc_random1 <-  gibbs_mala_sampler_ppt( y = simData$count_data, X = simData$design_matrix, 
		tau = 0.00035, nu2 = 100,  K = Kmax, dirPriorAlphas = dirPriorAlphas,
		mcmc_cycles = mcmc_cycles, iter_per_cycle = iter_per_cycle, 
		start_values = 'RANDOM', 
		nChains = nChains, nCores = nCores, warm_up = warm_up, showGraph = 1000, 
		checkAR = 1000)

#sampled values for the number of clusters (non-empty mixture components) per chain (columns)
mcmc_random1$nClusters
}
  
