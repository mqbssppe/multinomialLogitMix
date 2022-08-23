pkgname <- "multinomialLogitMix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "multinomialLogitMix-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('multinomialLogitMix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("gibbs_mala_sampler")
### * gibbs_mala_sampler

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gibbs_mala_sampler
### Title: The core of the Hybrid Gibbs/MALA MCMC sampler for the
###   multinomial logit mixture.
### Aliases: gibbs_mala_sampler

### ** Examples

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gibbs_mala_sampler", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gibbs_mala_sampler_ppt")
### * gibbs_mala_sampler_ppt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gibbs_mala_sampler_ppt
### Title: Prior parallel tempering scheme of hybrid Gibbs/MALA MCMC
###   samplers for the multinomial logit mixture.
### Aliases: gibbs_mala_sampler_ppt

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gibbs_mala_sampler_ppt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mix_mnm_logistic")
### * mix_mnm_logistic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mix_mnm_logistic
### Title: EM algorithm
### Aliases: mix_mnm_logistic

### ** Examples

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
	



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mix_mnm_logistic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multinomialLogitMix")
### * multinomialLogitMix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multinomialLogitMix
### Title: Main function
### Aliases: multinomialLogitMix

### ** Examples

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
	## No test: 
	mlm <- multinomialLogitMix(response = simData$count_data, 
		design_matrix = simData$design_matrix, method = "MCMC", 
             Kmax = 5, nCores = 2, splitSmallEM = TRUE, 
             mcmc_parameters = mcmc_parameters, em_parameters = em_parameters)
	# retrieve clustering according to EM
	mlm$EM$estimated_clustering
	# retrieve clustering according to MCMC
	mlm$MCMC_post_processed$cluster	
	
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multinomialLogitMix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
