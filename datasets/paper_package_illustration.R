library(multinomialLogitMix)
set.seed(727)

# generate a synthetic dataset 
n <- 250
p <- 3  # number of covariates (constant term included)
D <- 6  # number of multinomial categories
K <- 2  #number of clusters

simData <- simulate_multinomial_data(K = K, p= p, D = D, n = n , size = 20, prob = 0.025)   

# observed data
y <- simData$count_data	#multinomial count data matrix (response)
X <- simData$design_matrix	#desing matrix

# run the proposed methods (EM+MCMC) using 
#	EM algorithm with split-shake-random initialization)
#	and then use the selected model (according to ICL)
#	in order to initialize the MCMC sampler using an 
#	overfitting mixture model with Kmax = 10 components
#	(maximum number of clusters)
#	We will use nCores = 8 cores in parallel. 
#	All other parameters are set to their default values.

Kmax <- 10
y <- simData$count_data
X <- simData$design_matrix
nCores <- 8
#	Run the algorithm
#	NOTE: this will first run the EM algorithm for models with 1, 2, ..., 10 components
#		and then an MCMC sampler using an overfitting mixture of 10 components,
#		using a prior papalell tempering scheme of 8 chains in total. 
#		consisting of 100000 iterations, where the first half is discarded as burn-in
#		and the second half is retained for inference (after thinning the chain)
#		and post-processing for correcting label-switching.

mlm_split <- multinomialLogitMix(response = y, 
          design_matrix = X, method = "MCMC", 
          Kmax = Kmax, nCores = nCores, splitSmallEM = TRUE)


mlm_split$EM$estimated_K

mlm_split$EM$estimated_clustering

mlm_split$EM$all_runs[[2]]$weights

mlm_split$EM$all_runs[[2]]$beta

# mcmc
mlm_split$MCMC_post_processed$nClusters_posterior

mlm_split$MCMC_post_processed$cluster


# post-processed simulated allocation vectors
mlm_split$MCMC_post_processed$mcmc[[1]]

# post-processed beta_{k=2, j = 4, p = 2}
beta_123 <- mlm_split$MCMC_post_processed$mcmc[[2]][[3]][,1,2]
beta_223 <- mlm_split$MCMC_post_processed$mcmc[[2]][[3]][,2,2]
betas <- as.mcmc(cbind(beta_123, beta_223))
summary(betas)
par(mar = c(4,4,2,1))
plot(betas)

round(mlm_split$EM$all_runs[[mlm_split$EM$estimated_K]]$posteriorProbabilities, 2)

round(mlm_split$MCMC_post_processed$posteriorProbabilities, 2)
library(RColorBrewer)
myCol <- brewer.pal(12, "Set3")[-2]
par(mar = c(4,5,1,4))
matplot(mlm_split$MCMC_raw$coefficients[,2,3,], type = "p", pch = 16, col = myCol, ylab = bquote("MCMC values of "*beta[k*","*2*","*3]), xlab = "iteration")
abline(h = mlm_split$EM$all_runs[[mlm_split$EM$estimated_K ]]$beta[2,3,], lty = 2)
axis(4, at = mlm_split$EM$all_runs[[mlm_split$EM$estimated_K ]]$beta[2,3,], labels = c(expression(hat(beta)[1*","*2*","*3]^"(EM)"), expression(hat(beta)[2*","*2*","*3]^"(EM)")), las = 2)
legend("bottomleft", paste0("k = ", 1:10), col = myCol[1:10], pch = 16, ncol = 5, title = "mixture component")


