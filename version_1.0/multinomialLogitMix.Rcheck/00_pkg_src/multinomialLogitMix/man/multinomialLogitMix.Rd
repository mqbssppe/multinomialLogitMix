\name{multinomialLogitMix}
\alias{multinomialLogitMix}
\title{
Main function
}
\description{
The main function of the package. 
}
\usage{
multinomialLogitMix(response, design_matrix, method, 
	Kmax = 10, mcmc_parameters = NULL, em_parameters = NULL, 
	nCores, splitSmallEM = TRUE)
}
\arguments{
  \item{response}{
matrix of counts.
}
  \item{design_matrix}{
design matrix (including constant term). 
}
  \item{method}{
character with two possible values: "EM" or "MCMC" indicating the desired method in order to estimate the model. 
}
  \item{Kmax}{
number of components of the (overfitting) mixture model.
}
  \item{nCores}{
Total number of CPU cores for parallel processing.
}
\item{mcmc_parameters}{List with the parameter set-up of the MCMC sampler. See details for changing the defaults.}
\item{em_parameters}{List with the parameter set-up of the EM algorithm. See details  for changing the defaults.}
\item{splitSmallEM}{Boolean value, indicating whether the split-small EM scheme should be used to initialize the \code{method}. Default: true (suggested).}
}
\details{
The details of the parameter setup of the EM algorithm and MCMC sampler. The following specification correspond to the minimal default settings. Larger values of \code{tsplit} will result to better performance. 

em_parameters <- list(maxIter = 100, emthreshold = 1e-08, 
    maxNR = 10, tsplit = 16, msplit = 10, split = TRUE, 
    R0 = 0.1, plotting = TRUE)

mcmc_parameters <- list(tau = 0.00035, nu2 = 100, mcmc_cycles = 2600, 
    iter_per_cycle = 20, nChains = 8, dirPriorAlphas = c(1, 
        1 + 5 * exp((seq(2, 14, length = nChains - 1)))/100)/(200), 
    warm_up = 48000, checkAR = 500, probsSave = FALSE, 
    showGraph = 100, ar_low = 0.15, ar_up = 0.25, burn = 100, 
    thin = 1, withRandom = TRUE)



}
\value{
\item{EM}{List with the results of the EM algorithm.}
\item{MCMC_raw}{List with the raw output of the MCMC sampler - not identifiable MCMC output.}
\item{MCMC_post_processed}{Post-processed MCMC, used for the inference.}
}
\references{
Papastamoulis, P. (2022). Model-based clustering of replicated multinomial data. arXiv:2207.13984 [stat.ME]
}
\author{
Panagiotis Papastamoulis
}
\examples{
#	Generate synthetic data

	K <- 2	#number of clusters
	p <- 2	#number of covariates (constant incl)
	D <- 5	#number of categories
	n <- 20 #generated number of observations
	set.seed(1)
	simData <- simulate_multinomial_data(K = K, p = p, D = D, n = n, size = 20, prob = 0.025)   


	# EM parameters
em_parameters <- list(maxIter = 100, emthreshold = 1e-08, 
    maxNR = 10, tsplit = 16, msplit = 10, split = TRUE, 
    R0 = 0.1, plotting = TRUE)

	#  MCMC parameters - just for illustration
	#	typically, set `mcmc_cycles` and `warm_up`to a larger values
	#	such as` mcmc_cycles = 2500` or more 
	#	and `warm_up = 40000` or more.
	nChains <- 2 #(set this to a larger value, such as 8 or more)
	mcmc_parameters <- list(tau = 0.00035, nu2 = 100, mcmc_cycles = 260, 
	    iter_per_cycle = 20, nChains = nChains, dirPriorAlphas = c(1, 
		1 + 5 * exp((seq(2, 14, length = nChains - 1)))/100)/(200), 
	    warm_up = 4800, checkAR = 500, probsSave = FALSE, 
	    showGraph = 100, ar_low = 0.15, ar_up = 0.25, burn = 100, 
	    thin = 1, withRandom = TRUE)

	# run EM with split-small-EM initialization, and then use the output to 
	#	initialize MCMC algorithm for an overfitting mixture with 
	#	Kmax = 5 components (max number of clusters - usually this is 
	#	set to a larger value, e.g. 10 or 20).
	#	Note: 
	#		1. the MCMC output is based on the non-empty components
	#		2. the EM algorithm clustering corresponds to the selected 
	#			number of clusters according to ICL.
	#		3. `nCores` should by adjusted according to your available cores.
	\donttest{
	mlm <- multinomialLogitMix(response = simData$count_data, 
		design_matrix = simData$design_matrix, method = "MCMC", 
             Kmax = 5, nCores = 2, splitSmallEM = TRUE, 
             mcmc_parameters = mcmc_parameters, em_parameters = em_parameters)
	# retrieve clustering according to EM
	mlm$EM$estimated_clustering
	# retrieve clustering according to MCMC
	mlm$MCMC_post_processed$cluster	
	}
}
  
