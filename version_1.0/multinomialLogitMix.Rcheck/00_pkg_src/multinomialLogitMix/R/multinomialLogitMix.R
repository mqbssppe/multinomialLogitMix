# auto den xrisimopoieitai kapou
expected_complete_LL <- function(y, X, b, w, pr){
	n <- dim(w)[1]
	p <- dim(b)[2]
	K <- dim(b)[3]
	D <- dim(b)[1]+1
	
	G <- array(data = NA, dim = c(n, D - 1, K))
	s <- rowSums(y)
	sum1 <- matrix(0, nrow = n, ncol=D)
	cll <- 0
	for(k in 1:K){
		for(d in 1:(D-1)){
			G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
			sum1[,d+1] <- G[,d,k]
		}
		lastCol <- rowLogSumExps(sum1)
		G[,,k] <- G[,,k] - lastCol
		#                       G <- G/sum1
		theta[ , 1:(D-1), k] <- G[,,k]
		theta[ , D, k] <- -lastCol
		cll <- cll + w[,k] * ( log(pr[k]) + rowSums(y * theta[,,k]) )

	}       
	cll <- sum(cll)

	return(cll)
}


newton_raphson_mstep <- function(y, X, b, w, maxNR = 5, R0 = 0.1, method = 5, verbose = FALSE){

		if(method == 5){
#		ridge + adjust step size + complete log-ll evaluation
	#	Input
	#---------------------------------------------------------------------------
	#	b:	Regression coefficients on the previous iteration of the EM
	#		format:	(D-1) x p x K array		
	#	w:	posterior probabilities computed at the E-step
	#		format:	n x K matrix (row-sums=1)
	#	maxNR:	Maximum number of NR loops
	#		format:	positive integer
	#---------------------------------------------------------------------------
	#
	#	Output
	#
	#---------------------------------------------------------------------------
	#	b:	updated values of b
	#	theta:	corresponding multinomial prob per subject (n x D x K array)
	#---------------------------------------------------------------------------
		n <- dim(w)[1]
		p <- dim(b)[2]
		K <- dim(b)[3]
		D <- dim(b)[1]+1
		G <- array(data = NA, dim = c(n, D - 1, K))
		s <- rowSums(y)
		normConst <- lgamma(s + 1) - rowSums(lgamma(y + 1))
		FirstDerivative <- array(data = NA, dim = c(p, D-1) )
		SecondDerivative <- array(data = NA, dim = c((D-1)*p, (D-1)*p) )
		theta <- array(data = NA, dim = c(n, D, K))
		llValues <- numeric(maxNR+1)
		terminateNR <- FALSE
		nr_iterations <- 0
		sc <- numeric(maxNR)
		sum1 <- matrix(0, nrow = n, ncol=D)
		bStart <- b

		pr <- apply(w, 2, mean)
		pr <- pr/sum(pr)
		unityMatrix <- diag(p*(D-1))
		# initial state
		bLongOld <- c()
		cll <- 0
		for(k in 1:K){
			for(d in 1:(D-1)){
				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				sum1[,d+1] <- G[,d,k]
			}
			lastCol <- rowLogSumExps(sum1)
			G[,,k] <- G[,,k] - lastCol
	#			G <- G/sum1
			theta[ , 1:(D-1), k] <- G[,,k]
			theta[ , D, k] <- -lastCol
			bLongOld <- c(bLongOld, b[,,k])
			cll <- cll + w[,k] * ( log(pr[k]) + rowSums(y * theta[,,k]) + normConst)
		}	
		logLikelihoodOld <- sum(cll)
		llValues[1] <- logLikelihoodOld			
		thetaOld <- theta
		Gold <- G
	#	theta1 <- theta
		bOld <- b
		if(verbose){print(logLikelihoodOld)}
		R <- R0
		llDIFF <- 100
		nr_iterations <- 0
		SD_long <- matrix(0, K*(D-1)*p,  K*(D-1)*p)
		while((nr_iterations < maxNR) & (abs(llDIFF) > 1e-3)){
			nr_iterations	<- nr_iterations + 1
			FD_long <- c()
			bLongNew <- c()
			cll <- 0
			for(k in 1:K){
	#			for(d in 1:(D-1)){
	#				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
	#				sum1[,d+1] <- G[,d,k]
	#			}
	#			lastCol <- rowLogSumExps(sum1)
	#			G[,,k] <- G[,,k] - lastCol
	#			theta[ , 1:(D-1), k] <- G[,,k]
	#			theta[ , D, k] <- -lastCol
	#			theta2 <- theta
	#			print(sum(abs(theta2-theta1)))
				G[,,k] <- exp(G[,,k])
				motorhead <- matrix(rep(w[,k],D-1),ncol=D-1) * (y[,-D] - s*G[,,k])
				protector <- s * w[,k] * G[,,k]
				for(r in 1:p){
					FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
					for(h in 1:p){
						nocturnus <- matrix(protector*X[,r]*X[,h], n, D - 1)
						for(d in 1:(D-1)){
							i <- (r-1)*(D-1) + d
							for(d2 in 1:(D-1)){
								j <- (h-1)*(D-1) + d2
								if(j >= i){
									kron_delta <- 0
									if(d2 == d){kron_delta = 1}
									SecondDerivative[i, j] <- sum(
										nocturnus[,d]*(kron_delta - G[,d2,k]) )
								}else{
									SecondDerivative[i, j] <- SecondDerivative[j, i]
								}
							}
						}
					}
				}
				SecondDerivative <- -SecondDerivative
				SecondDerivative[is.na(SecondDerivative)] <- 0
				grammes <- ((k - 1) * (D-1) * p + 1):(k * (D-1) * p) 
				stiles  <- grammes
				SD_long[grammes, stiles] <- SecondDerivative
				# compute alpha 
				lambda.1 <- max(eigen(SecondDerivative)$values)
				norm.SS <-sqrt(sum(FirstDerivative^2))
				alpha <- lambda.1 + R * norm.SS		
				# compute modified Hessian
				BB <- SecondDerivative - (alpha > 0) * alpha *unityMatrix			
				FD_long <- c(FD_long, c(t(FirstDerivative)))
				# N-R update
				sepultura <- tryCatch(
					qr.coef(qr(BB, tol = 1e-300), c(t(FirstDerivative))),
					error=function(e)return(0)
				)
				bVector <- c(b[,,k]) - sepultura
				bVectorOld <- as.vector(b[,,k])
				b[,,k] <- matrix(bVector, nrow = D-1, ncol = p)
				bLongNew <- c(bLongNew, b[,,k])
				for(d in 1:(D-1)){
					G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
					sum1[,d+1] <- G[,d,k]
				}
				lastCol <- rowLogSumExps(sum1)
				G[,,k] <- G[,,k] - lastCol
				theta[ , 1:(D-1), k] <- G[,,k]
				theta[ , D, k] <- -lastCol
	#			theta1 <- theta
				cll <- cll + w[,k] * ( log(pr[k]) + rowSums(y * theta[,,k]) + normConst)
				
			}
			logLikelihoodNew <- sum(cll)
			llValues[nr_iterations + 1] <- logLikelihoodNew
			llDIFF <- logLikelihoodNew - logLikelihoodOld 
#			if (llDIFF < 0 ){
#				
#			}
			# adjust R
			dTheta <- bLongNew - bLongOld
			ll.LQA <- logLikelihoodOld + FD_long %*% dTheta + 
				(1/2)*t(dTheta) %*% SD_long %*% dTheta
			delta.ll <- llDIFF
			delta.ll.LQA <- ll.LQA - logLikelihoodOld
			z <- delta.ll / delta.ll.LQA

	# na valw an logdiff == 0 na stamataei 
			if(delta.ll.LQA != 0){
				if( z < 0){
					nr_iterations	<- nr_iterations - 1
					R <- R * 2
					G <- Gold
					theta <- thetaOld
					b <- bOld
				}else{
		#			if(z < 0.1){
		#				R <- R * 0.5
		#
		#			} 

					if((.7 < z) & (z < 1.3)){
						R <- R * 0.5

					} else if(z > 2){
						R <- R * 2
					}
					logLikelihoodOld <- logLikelihoodNew
					bLongOld <- bLongNew
					bOld <- b
					thetaOld <- theta
					Gold <- G
					if(verbose){cat(paste0('i = ', nr_iterations, ', logL = ', round(logLikelihoodNew,2), ', R = ', R, 
						', llDIFF = ', round(llDIFF,2), ', z = ', round(z,3)), '\n')}			

				} 
			}
		
		}

		
		
		result <- vector('list', length=3)
		result[[1]] <- b
		result[[2]] <- theta
		result[[3]] <- llValues[1:(nr_iterations+1)]
		names(result) <- c('b','theta','ll')
		return(result)

	}else if(method == 4){
#		ridge + adjust step size + observed log-ll evaluation
	#	Input
	#---------------------------------------------------------------------------
	#	b:	Regression coefficients on the previous iteration of the EM
	#		format:	(D-1) x p x K array		
	#	w:	posterior probabilities computed at the E-step
	#		format:	n x K matrix (row-sums=1)
	#	maxNR:	Maximum number of NR loops
	#		format:	positive integer
	#---------------------------------------------------------------------------
	#
	#	Output
	#
	#---------------------------------------------------------------------------
	#	b:	updated values of b
	#	theta:	corresponding multinomial prob per subject (n x D x K array)
	#---------------------------------------------------------------------------
		n <- dim(w)[1]
		p <- dim(b)[2]
		K <- dim(b)[3]
		D <- dim(b)[1]+1
		G <- array(data = NA, dim = c(n, D - 1, K))
		s <- rowSums(y)
		FirstDerivative <- array(data = NA, dim = c(p, D-1) )
		SecondDerivative <- array(data = NA, dim = c((D-1)*p, (D-1)*p) )
		theta <- array(data = NA, dim = c(n, D, K))
		terminateNR <- FALSE
		nr_iterations <- 0
		sc <- numeric(maxNR)
		sum1 <- matrix(0, nrow = n, ncol=D)
		bStart <- b

		pr <- apply(w, 2, mean)
		pr <- pr/sum(pr)
		unityMatrix <- diag(p*(D-1))
		# initial state
		bLongOld <- c()
		for(k in 1:K){
			for(d in 1:(D-1)){
				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				sum1[,d+1] <- G[,d,k]
			}
			lastCol <- rowLogSumExps(sum1)
			G[,,k] <- G[,,k] - lastCol
	#			G <- G/sum1
			theta[ , 1:(D-1), k] <- G[,,k]
			theta[ , D, k] <- -lastCol
			bLongOld <- c(bLongOld, b[,,k])
		}	
		logLikelihoodOld <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
		thetaOld <- theta
		Gold <- G
	#	theta1 <- theta
		bOld <- b
		if(verbose){print(logLikelihoodOld)}
		R <- R0
		llDIFF <- 100
		nr_iterations <- 0
		SD_long <- matrix(0, K*(D-1)*p,  K*(D-1)*p)
		while((nr_iterations < maxNR) & (abs(llDIFF) > 1e-3)){
			nr_iterations	<- nr_iterations + 1
			FD_long <- c()
			bLongNew <- c()
			for(k in 1:K){
	#			for(d in 1:(D-1)){
	#				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
	#				sum1[,d+1] <- G[,d,k]
	#			}
	#			lastCol <- rowLogSumExps(sum1)
	#			G[,,k] <- G[,,k] - lastCol
	#			theta[ , 1:(D-1), k] <- G[,,k]
	#			theta[ , D, k] <- -lastCol
	#			theta2 <- theta
	#			print(sum(abs(theta2-theta1)))
				G[,,k] <- exp(G[,,k])
				motorhead <- matrix(rep(w[,k],D-1),ncol=D-1) * (y[,-D] - s*G[,,k])
				protector <- s * w[,k] * G[,,k]
				for(r in 1:p){
					FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
					for(h in 1:p){
						nocturnus <- matrix(protector*X[,r]*X[,h], n, D - 1)
						for(d in 1:(D-1)){
							i <- (r-1)*(D-1) + d
							for(d2 in 1:(D-1)){
								j <- (h-1)*(D-1) + d2
								if(j >= i){
									kron_delta <- 0
									if(d2 == d){kron_delta = 1}
									SecondDerivative[i, j] <- sum(
										nocturnus[,d]*(kron_delta - G[,d2,k]) )
								}else{
									SecondDerivative[i, j] <- SecondDerivative[j, i]
								}
							}
						}
					}
				}
				SecondDerivative <- -SecondDerivative
				SecondDerivative[is.na(SecondDerivative)] <- 0
				grammes <- ((k - 1) * (D-1) * p + 1):(k * (D-1) * p) 
				stiles  <- grammes
				SD_long[grammes, stiles] <- SecondDerivative
				# compute alpha 
				lambda.1 <- max(eigen(SecondDerivative)$values)
				norm.SS <-sqrt(sum(FirstDerivative^2))
				alpha <- lambda.1 + R * norm.SS		
				# compute modified Hessian
				BB <- SecondDerivative - (alpha > 0) * alpha *unityMatrix			
				FD_long <- c(FD_long, c(t(FirstDerivative)))
				# N-R update
				sepultura <- tryCatch(
					qr.coef(qr(BB, tol = 1e-300), c(t(FirstDerivative))),
					error=function(e)return(0)
				)
				bVector <- c(b[,,k]) - sepultura
				bVectorOld <- as.vector(b[,,k])
				b[,,k] <- matrix(bVector, nrow = D-1, ncol = p)
				bLongNew <- c(bLongNew, b[,,k])
				for(d in 1:(D-1)){
					G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
					sum1[,d+1] <- G[,d,k]
				}
				lastCol <- rowLogSumExps(sum1)
				G[,,k] <- G[,,k] - lastCol
				theta[ , 1:(D-1), k] <- G[,,k]
				theta[ , D, k] <- -lastCol
	#			theta1 <- theta
				
			}
			logLikelihoodNew <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
			llDIFF <- logLikelihoodNew - logLikelihoodOld 
			if (llDIFF < 0 ){
				
			}
			# adjust R
			dTheta <- bLongNew - bLongOld
			ll.LQA <- logLikelihoodOld + FD_long %*% dTheta + 
				(1/2)*t(dTheta) %*% SD_long %*% dTheta
			delta.ll <- llDIFF
			delta.ll.LQA <- ll.LQA - logLikelihoodOld
			z <- delta.ll / delta.ll.LQA

	# na valw an logdiff == 0 na stamataei 
			if(delta.ll.LQA != 0){
				if( z < 0){
					nr_iterations	<- nr_iterations - 1
					R <- R * 2
					G <- Gold
					theta <- thetaOld
					b <- bOld
				}else{
		#			if(z < 0.1){
		#				R <- R * 0.5
		#
		#			} 

					if((.7 < z) & (z < 1.3)){
						R <- R * 0.5

					} else if(z > 2){
						R <- R * 2
					}
					logLikelihoodOld <- logLikelihoodNew
					bLongOld <- bLongNew
					bOld <- b
					thetaOld <- theta
					Gold <- G
					if(verbose){cat(paste0('i = ', nr_iterations, ', logL = ', round(logLikelihoodNew,2), ', R = ', R, 
						', llDIFF = ', round(llDIFF,2), ', z = ', round(z,3)), '\n')}			

				} 
			}
		
		}

		
		
		result <- vector('list', length=2)
		result[[1]] <- b
		result[[2]] <- theta
		names(result) <- c('b','theta')
		return(result)

	}else if(method == 3){
#		ridge + no adjust step size
	#	Input
	#---------------------------------------------------------------------------
	#	b:	Regression coefficients on the previous iteration of the EM
	#		format:	(D-1) x p x K array		
	#	w:	posterior probabilities computed at the E-step
	#		format:	n x K matrix (row-sums=1)
	#	maxNR:	Maximum number of NR loops
	#		format:	positive integer
	#---------------------------------------------------------------------------
	#
	#	Output
	#
	#---------------------------------------------------------------------------
	#	b:	updated values of b
	#	theta:	corresponding multinomial prob per subject (n x D x K array)
	#---------------------------------------------------------------------------
		n <- dim(w)[1]
		p <- dim(b)[2]
		K <- dim(b)[3]
		D <- dim(b)[1]+1
		G <- array(data = NA, dim = c(n, D - 1, K))
		s <- rowSums(y)
		FirstDerivative <- array(data = NA, dim = c(p, D-1) )
		SecondDerivative <- array(data = NA, dim = c((D-1)*p, (D-1)*p) )
		theta <- array(data = NA, dim = c(n, D, K))
		terminateNR <- FALSE
		nr_iterations <- 0
		sc <- numeric(maxNR)
		sum1 <- matrix(0, nrow = n, ncol=D)
		bStart <- b

		pr <- apply(w, 2, mean)
		pr <- pr/sum(pr)
		unityMatrix <- diag(p*(D-1))
		# initial state
		
		for(k in 1:K){
			for(d in 1:(D-1)){
				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				sum1[,d+1] <- G[,d,k]
			}
			lastCol <- rowLogSumExps(sum1)
			G[,,k] <- G[,,k] - lastCol
	#			G <- G/sum1
			theta[ , 1:(D-1), k] <- G[,,k]
			theta[ , D, k] <- -lastCol
		}	
		logLikelihoodOld <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
		thetaOld <- theta
		bOld <- b
		if(verbose){print(logLikelihoodOld)}
		R <- R0
		llDIFF <- 100
		nr_iterations <- 0
		while((nr_iterations < maxNR) & (llDIFF > 1e-3)){
			nr_iterations	<- nr_iterations + 1
			for(k in 1:K){
				#for(d in 1:(D-1)){
				#	G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				#	sum1[,d+1] <- G[,d,k]
				#}
				#lastCol <- rowLogSumExps(sum1)
				#G[,,k] <- G[,,k] - lastCol
				#theta[ , 1:(D-1), k] <- G[,,k]
				#theta[ , D, k] <- -lastCol
				G[,,k] <- exp(G[,,k])
				motorhead <- matrix(rep(w[,k],D-1),ncol=D-1) * (y[,-D] - s*G[,,k])
				protector <- s * w[,k] * G[,,k]
				for(r in 1:p){
					FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
					for(h in 1:p){
						nocturnus <- matrix(protector*X[,r]*X[,h], n, D - 1)
						for(d in 1:(D-1)){
							i <- (r-1)*(D-1) + d
							for(d2 in 1:(D-1)){
								j <- (h-1)*(D-1) + d2
								if(j >= i){
									kron_delta <- 0
									if(d2 == d){kron_delta = 1}
									SecondDerivative[i, j] <- sum(
										nocturnus[,d]*(kron_delta - G[,d2,k]) )
								}else{
									SecondDerivative[i, j] <- SecondDerivative[j, i]
								}
							}
						}
					}
				}
				SecondDerivative <- -SecondDerivative
				SecondDerivative[is.na(SecondDerivative)] <- 0
				
				# compute alpha 
				lambda.1 <- max(eigen(SecondDerivative)$values)
				norm.SS <-sqrt(sum(FirstDerivative^2))
				alpha <- lambda.1 + R * norm.SS		
				# compute modified Hessian
				BB <- SecondDerivative - (alpha > 0) * alpha *unityMatrix			
				
				# N-R update
				sepultura <- tryCatch(
					qr.coef(qr(BB, tol = 1e-300), c(t(FirstDerivative))),
					error=function(e)return(0)
				)
				bVector <- c(b[,,k]) - sepultura
				bVectorOld <- as.vector(b[,,k])
				b[,,k] <- matrix(bVector, nrow = D-1, ncol = p)
				
				for(d in 1:(D-1)){
					G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
					sum1[,d+1] <- G[,d,k]
				}
				lastCol <- rowLogSumExps(sum1)
				G[,,k] <- G[,,k] - lastCol
				theta[ , 1:(D-1), k] <- G[,,k]
				theta[ , D, k] <- -lastCol
				
			}
			logLikelihoodNew <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
			llDIFF <- logLikelihoodNew - logLikelihoodOld 
			# adjust R
			#dTheta <- bVector - bVectorOld
			#negll.LQA <- - (-logLikelihoodOld + c(t(FirstDerivative)) %*% dTheta + 
			#	(1/2)*t(dTheta) %*% SecondDerivative %*% dTheta)
			#delta.negll <- - llDIFF
			#delta.negll.LQA <- negll.LQA + logLikelihoodOld
			#z <- delta.negll / delta.negll.LQA
			#if( z < 0){R <- R * 4} else if((.7 < z) & (z < 1.3)){R <- R * 0.4} else if(z > 2){R <- R * 4}
			if(verbose){cat(paste0('logL = ', round(logLikelihoodNew,2), ', R = ', R, 
				', llDIFF = ', round(llDIFF,2)), '\n')}			
			logLikelihoodOld <- logLikelihoodNew

		
		}

		
		
		result <- vector('list', length=2)
		result[[1]] <- b
		result[[2]] <- theta
		names(result) <- c('b','theta')
		return(result)

	}else if(method == 2){
#		nr with adjust step size
	#	Input
	#---------------------------------------------------------------------------
	#	b:	Regression coefficients on the previous iteration of the EM
	#		format:	(D-1) x p x K array		
	#	w:	posterior probabilities computed at the E-step
	#		format:	n x K matrix (row-sums=1)
	#	maxNR:	Maximum number of NR loops
	#		format:	positive integer
	#---------------------------------------------------------------------------
	#
	#	Output
	#
	#---------------------------------------------------------------------------
	#	b:	updated values of b
	#	theta:	corresponding multinomial prob per subject (n x D x K array)
	#---------------------------------------------------------------------------
		n <- dim(w)[1]
		p <- dim(b)[2]
		K <- dim(b)[3]
		D <- dim(b)[1]+1
		G <- array(data = NA, dim = c(n, D - 1, K))
		s <- rowSums(y)
		FirstDerivative <- array(data = NA, dim = c(p, D-1) )
		SecondDerivative <- array(data = NA, dim = c((D-1)*p, (D-1)*p) )
		theta <- array(data = NA, dim = c(n, D, K))
		terminateNR <- FALSE
		nr_iterations <- 0
		sc <- numeric(maxNR)
		sum1 <- matrix(0, nrow = n, ncol=D)
		bStart <- b

		pr <- apply(w, 2, mean)
		pr <- pr/sum(pr)
		
		# initial state
		stepSize <- 1
		for(k in 1:K){
			for(d in 1:(D-1)){
				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				sum1[,d+1] <- G[,d,k]
			}
			lastCol <- rowLogSumExps(sum1)
			G[,,k] <- G[,,k] - lastCol
	#			G <- G/sum1
			theta[ , 1:(D-1), k] <- G[,,k]
			theta[ , D, k] <- -lastCol
			G[,,k] <- exp(G[,,k])
		}

		Gold <- G	
		logLikelihoodOld <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
		thetaOld <- theta
		bOld <- b
		#print(logLikelihoodOld)	
		reDO <- FALSE

		while(nr_iterations < maxNR & (terminateNR == FALSE)){
			nr_iterations <- nr_iterations + 1
			for(k in 1:K){
				#for(d in 1:(D-1)){
				#	G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				#	sum1[,d+1] <- G[,d,k]
				#}
				#lastCol <- rowLogSumExps(sum1)
				#G[,,k] <- G[,,k] - lastCol
				#theta[ , 1:(D-1), k] <- G[,,k]
				#theta[ , D, k] <- -lastCol
				#G[,,k] <- exp(G[,,k])
				if(reDO == FALSE){
					motorhead <- matrix(rep(w[,k],D-1),ncol=D-1) * (y[,-D] - s*G[,,k])
					protector <- s * w[,k] * G[,,k]
					for(r in 1:p){
						FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
						for(h in 1:p){
							nocturnus <- matrix(protector*X[,r]*X[,h], n, D - 1)
							for(d in 1:(D-1)){
								i <- (r-1)*(D-1) + d
								for(d2 in 1:(D-1)){
									j <- (h-1)*(D-1) + d2
									if(j >= i){
										kron_delta <- 0
										if(d2 == d){kron_delta = 1}
										SecondDerivative[i, j] <- sum(
											nocturnus[,d]*(kron_delta - G[,d2,k]) )
									}else{
										SecondDerivative[i, j] <- SecondDerivative[j, i]
									}
								}
							}
						}
					}
					SecondDerivative <- -SecondDerivative
					SecondDerivative[is.na(SecondDerivative)] <- 0
					# N-R update
					sepultura <- tryCatch(
						qr.coef(qr(SecondDerivative, tol = 1e-300), c(t(FirstDerivative))),
						error=function(e)return(0)
					)
				}
				bVector <- c(b[,,k]) - sepultura*stepSize
				b[,,k] <- matrix(bVector, nrow = D-1, ncol = p)
				#
				
			}
			if(is.na(min(range(b)))){
				#print('NA in NR')
				terminateNR = TRUE
			}
			b[is.na(b)] = 0
			theta[is.na(theta)] <- -745
			b[abs(b) > 10000] <- 0
			
			for(k in 1:K){
				for(d in 1:(D-1)){
					G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
					sum1[,d+1] <- G[,d,k]
				}
				lastCol <- rowLogSumExps(sum1)
				G[,,k] <- G[,,k] - lastCol
				theta[ , 1:(D-1), k] <- G[,,k]
				theta[ , D, k] <- -lastCol
				G[,,k] <- exp(G[,,k])
			}	


			logLikelihoodNew <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
			derivNormDiff <- logLikelihoodOld - logLikelihoodNew

			if( (derivNormDiff > 0) | (logLikelihoodNew == 0)  ){
			#	cat('[WARNING] Newton-Raphson adjusted because logL decreased.','\n')
				b <- bOld
				theta <- thetaOld
				G <- Gold
				stepSize <- stepSize/2
				nr_iterations <- nr_iterations - 1
				reDO <- TRUE
			}else{
				if((abs(derivNormDiff) < 1e-6) | (terminateNR == TRUE)){
					terminateNR <- TRUE
				}else{
					bOld <- b
					thetaOld <- theta
					Gold <- G
					logLikelihoodOld <- logLikelihoodNew
					reDO <- FALSE
				}
			}
			#print(c(nr_iterations, logLikelihoodNew, stepSize))
			#print(b)
		}
		
		
		result <- vector('list', length=2)
		result[[1]] <- b
		result[[2]] <- theta
		names(result) <- c('b','theta')
		return(result)

	}else if(method == 1){
#		nr no adjust step size
	#	Input
	#---------------------------------------------------------------------------
	#	b:	Regression coefficients on the previous iteration of the EM
	#		format:	(D-1) x p x K array		
	#	w:	posterior probabilities computed at the E-step
	#		format:	n x K matrix (row-sums=1)
	#	maxNR:	Maximum number of NR loops
	#		format:	positive integer
	#---------------------------------------------------------------------------
	#
	#	Output
	#
	#---------------------------------------------------------------------------
	#	b:	updated values of b
	#	theta:	corresponding multinomial prob per subject (n x D x K array)
	#---------------------------------------------------------------------------
		n <- dim(w)[1]
		p <- dim(b)[2]
		K <- dim(b)[3]
		D <- dim(b)[1]+1
		G <- array(data = NA, dim = c(n, D - 1, K))
		s <- rowSums(y)
		FirstDerivative <- array(data = NA, dim = c(p, D-1) )
		SecondDerivative <- array(data = NA, dim = c((D-1)*p, (D-1)*p) )
		theta <- array(data = NA, dim = c(n, D, K))
		llValues <- numeric(maxNR+1)
		terminateNR <- FALSE
		nr_iterations <- 0
		sc <- numeric(maxNR)
		sum1 <- matrix(0, nrow = n, ncol=D)
		bStart <- b

		pr <- apply(w, 2, mean)
		pr <- pr/sum(pr)
		
		# initial state
		
		for(k in 1:K){
			for(d in 1:(D-1)){
				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
				sum1[,d+1] <- G[,d,k]
			}
			lastCol <- rowLogSumExps(sum1)
			G[,,k] <- G[,,k] - lastCol
	#			G <- G/sum1
			theta[ , 1:(D-1), k] <- G[,,k]
			theta[ , D, k] <- -lastCol
		}	
		logLikelihoodOld <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
		llValues[1] <- logLikelihoodOld
		if(verbose){cat(paste0('iter = ', 0, ', ll = ', round(logLikelihoodOld,2)),'\n')}
		thetaOld <- theta
		bOld <- b
	#	print(logLikelihoodOld)	
		while(nr_iterations < maxNR & (terminateNR == FALSE)){
			nr_iterations <- nr_iterations + 1
			for(k in 1:K){
				for(d in 1:(D-1)){
					G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
					sum1[,d+1] <- G[,d,k]
				}
				lastCol <- rowLogSumExps(sum1)
				G[,,k] <- G[,,k] - lastCol
				theta[ , 1:(D-1), k] <- G[,,k]
				theta[ , D, k] <- -lastCol
				G[,,k] <- exp(G[,,k])
				motorhead <- matrix(rep(w[,k],D-1),ncol=D-1) * (y[,-D] - s*G[,,k])
				protector <- s * w[,k] * G[,,k]
				for(r in 1:p){
					FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
					for(h in 1:p){
						nocturnus <- matrix(protector*X[,r]*X[,h], n, D - 1)
						for(d in 1:(D-1)){
							i <- (r-1)*(D-1) + d
							for(d2 in 1:(D-1)){
								j <- (h-1)*(D-1) + d2
								if(j >= i){
									kron_delta <- 0
									if(d2 == d){kron_delta = 1}
									SecondDerivative[i, j] <- sum(
										nocturnus[,d]*(kron_delta - G[,d2,k]) )
								}else{
									SecondDerivative[i, j] <- SecondDerivative[j, i]
								}
							}
						}
					}
				}
				SecondDerivative <- -SecondDerivative
				# N-R update
				sepultura <- tryCatch(
					qr.coef(qr(SecondDerivative, tol = 1e-10), c(t(FirstDerivative))),
					error=function(e)return(0)
				)
				bVector <- c(b[,,k]) - sepultura
				b[,,k] <- matrix(bVector, nrow = D-1, ncol = p)
				#
				
			}
			
			b[is.na(b)] <- 0
			theta[is.na(theta)] <- -745
			b[abs(b) > 10000] <- 0
			
			for(k in 1:K){
				for(d in 1:(D-1)){
					G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
					sum1[,d+1] <- G[,d,k]
				}
				lastCol <- rowLogSumExps(sum1)
				G[,,k] <- G[,,k] - lastCol
				theta[ , 1:(D-1), k] <- G[,,k]
				theta[ , D, k] <- -lastCol
			}	


			logLikelihoodNew <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
			llValues[nr_iterations+1] <- logLikelihoodNew
			if(verbose){cat(paste0('iter = ', nr_iterations, ', ll = ', 
					round(logLikelihoodNew, 2)),'\n')}
			derivNormDiff <- logLikelihoodOld - logLikelihoodNew
			if( is.na(derivNormDiff) ){
	#			cat('[WARNING] Newton-Raphson terminated early due to underflows.','\n')
				terminateNR <- TRUE
				b <- bOld
				theta <- thetaOld
			}else{
				if( derivNormDiff > 0 ){
	#				cat('[WARNING] Newton-Raphson terminated because logL decreased.','\n')
					terminateNR <- TRUE
					b <- bOld
					theta <- thetaOld
				}
				if(abs(derivNormDiff) < 1e-6){
					terminateNR <- TRUE
				}else{
					bOld <- b
					thetaOld <- theta
					logLikelihoodOld <- logLikelihoodNew
				}
			}
			#print(logLikelihoodNew)
		}
		
		
		result <- vector('list', length=3)
		result[[1]] <- b
		result[[2]] <- theta
		result[[3]] <- llValues[1:(nr_iterations+1)]
		names(result) <- c('b','theta','ll')
		return(result)

	}
	
	
}




mixLoglikelihood_GLM <- function (y, theta, pi){
	# theta is in log-scale!
#    if (is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
#        stop(paste(sQuote("y"), "must be a matrix"))
#    if (min(y) < 0 | sum(round(y)) != sum(y)) 
#        stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
#    if (length(mean) != length(pi)) 
#        stop(paste(sQuote("mean"), "must be a list of the same length as", 
#            sQuote("pi")))
    s <- rowSums(y)
    log_multinom_coefficient <-  lgamma(s+1) - rowSums(lgamma(y + 1))
    g <- length(pi)
    n <- dim(y)[1]
    cols <- dim(y)[2]
    nas <- 0
    y <- matrix(y, nrow = n, ncol = cols)
    logLike <- rep(0, n)
    index <- 1:g
    epsilon <- exp(-720)
    thresh <- -745
    nn <- 0
    logpi <- log(pi)
    ef <- matrix(logpi, nrow = n, ncol = g, byrow = T)
    for (k in 1:g) {
        ef[, k] <- ef[, k] + rowSums( y * theta[,,k] ) + log_multinom_coefficient
    }
    efmax <- apply(ef, 1, max)
    wf_k <- ef
    ef <- ef - efmax
    logLike <- efmax + log(rowSums(exp(ef)))
    return(
	list(
		ll = sum(logLike, na.rm = TRUE),
		log_weighted_density_k = wf_k
		)
	)
}


multinomial_logistic_EM <- function(y, x, K, w_start, b_start, maxIter = 1000, emthreshold = 1e-8, maxNR = 5, nCores = NULL, verbose = FALSE, R0, method){
#	y:	multinomial response matrix
#	x:	design matrix
	n <- dim(y)[1]
	D <- dim(y)[2]
	p <- dim(x)[2]
	epsilon <- 1e-10
	thresh <- -744
	loglik <- numeric(maxIter)

	b <-  array( data = rnorm((D-1)*p*K, mean = 0, sd = 0.1), dim = c(D-1, p, K) ) #array( data = 0, dim = c(D-1, p, K) )
	theta <- array(data = NA, dim = c(n, D, K))

	if(K == 1){
		w <- t(t(rep(1,n)))
		mStep <- newton_raphson_mstep(y = y, X = x, b = b, w = w, maxNR = maxNR+5, R0 = R0, method = method)
		b <- mStep$b
		theta[,,1] <- mStep$theta
		loglik <- mixLoglikelihood_GLM(y, theta, pi = 1)$ll
		if(verbose){
		cat(paste0('     Log-likelihood: ', round(loglik,2)),'            ',' \n')
		}
		result <- vector('list',length = 4)
		result[[1]] <- 1
		result[[2]] <- b
		result[[3]] <- matrix(1,nrow = n, ncol = K)
		result[[4]] <- loglik
		names(result) <- c('weights', 'beta', 'posteriorProbabilities', 'logLikelihood')
		#cat('\n')
		return(result)
	}

#	INITIALIZATION
	if( missing(w_start) ){
		w_start <- matrix(runif(n*K), nrow = n, ncol = K)
		w_start <- t(apply(w_start, 1, function(z)z/sum(z)))
		# this is for k-means scheme
		#km <- kmeans(y, centers = K)
		#for(k in 1:K){
		#	ind <- which(km$cluster==k)
		#	w_start[ind,] <- 0
		#	w_start[ind,k] <- 1
		#}
	}
	w <- w_start
#	MIX WEIGHTS
	pr <- apply(w, 2, mean)
	pr <- pr/sum(pr)
#	REGRESSION COEFFICIENTS (EDW PREPEI NA EXW b=beta_start)
	if(missing(b_start)){
		b_start <- array(data=-0.2*runif((D-1)*p*K)+0.1,dim=c(D-1,p,K))		
	}
	b <- b_start
#	if(!is.null(nCores)){
#		registerDoParallel(nCores)
#		mStep <- newton_raphson_mstep2(y = y, X = x, b = b, w = w, maxNR = maxNR)
#	}else{
		mStep <- newton_raphson_mstep(y = y, X = x, b = b, w = w, maxNR = maxNR, R0 = R0, method = method)
#	}

	b <- mStep$b
	theta <- theta_start <- mStep$theta
	i <- 1
	LL <- mixLoglikelihood_GLM(y, theta, pr)
	loglik[i] <- LL$ll
	if(verbose){
		cat(paste0('     Log-likelihood at starting value: ', round(loglik[i],2)), '\n')
	}
	criterion <- 10^8

	while (criterion > emthreshold & i < maxIter) {
		i <- i + 1
	#	EXPECTATION STEP
		w <- matrix(data = 0, nrow = n, ncol = K, byrow = T)
		for (j in 1:K) {
		    w[, j] <- w[, j] + LL$log_weighted_density[,j]
		}
		if (K > 1) {
	            #print(head(w))
		    w[is.na(w)] <- thresh
		    v1 <- which(apply(w, 1, max) < thresh)
		    v3 <- 1:n
		    len <- length(v1)
		    if (len > 0) {
			v2 <- apply(array(w[v1, ], dim = c(len, K)), 1, order)[K, ]
			ddd <- cbind(v1, v2)
			w[v1, ] <- 0
			w[ddd] <- 1
			v3 <- -v1
		    }
		    w[v3, ] <- exp(w[v3, ])
		}
		w <- w/rowSums(w)
		#print(head(w))
		sl <- length(w[w < epsilon])
		bl <- length(w[w > 1 - epsilon])
		w[w < epsilon] <- rep(epsilon, sl)
		w[w > 1 - epsilon] <- rep(1 - epsilon, bl)
		w <- w/rowSums(w)
	#	MIX WEIGHTS
		pr <- apply(w, 2, mean)
	#	REGRESSION  COEFFICIENTS
#		if(is.null(nCores)){
			mStep <- newton_raphson_mstep(y = y, X = x, b = b, w = w, maxNR = maxNR, R0 = R0, method = method)
#		}else{
#			mStep <- newton_raphson_mstep2(y = y, X = x, b = b, w = w, maxNR = maxNR)
#		}
		b <- mStep$b
		theta <- mStep$theta
		LL <- mixLoglikelihood_GLM(y, theta, pr)
		loglik[i] <- LL$ll
		criterion <- abs((loglik[i - 1] - loglik[i])/loglik[i - 1])
	#	if(loglik[i - 1] - loglik[i] > 0){
	#		cat('[WARNING] Log-likelihood decreased. Quitting.','\n')
	#		b <- b_start
	#		theta <- theta_start
	#		loglik[i] <- loglik[1]
	#		criterion <- emthreshold - 1
	#	}
#		if((i %% 10) == 0){
#			plot(loglik[1:i], xlab = 'iteration', ylab = 'Log-likelihood',type='l',col='blue')
		if(verbose){
			cat(paste0('     Iteration: ', i, ', Log-likelihood: ', round(loglik[i],2)),'            ',' \n')
		}
#		}
	}
#	if(!is.null(nCores)){stopImplicitCluster()}
	result <- vector('list',length=4)
	result[[1]] <- pr
	result[[2]] <- b
	result[[3]] <- w
	result[[4]] <- loglik[1:i]
	names(result) <- c('weights', 'beta', 'posteriorProbabilities', 'logLikelihood')
	#cat('\n')
	return(result)

}



# this is the parallel version
splitEM_GLM <- function(y, x, K, smallerModel = NULL, tsplit = 10, maxIter = 20, emthreshold = 1e-8, maxNR=5, nCores, split = TRUE, R0, method){
	registerDoParallel(cores = nCores)
	n <- dim(y)[1]
	D <- dim(y)[2]
	p <- dim(x)[2]
	pr <- numeric(K)
	b <- array( data = 0, dim = c(D-1, p, K) )
	theta <- array(data = NA, dim = c(n, D, K))
	epsilon <- 1e-10
	thresh <- -744
	loglik <- numeric(maxIter)
#	if(is.null(smallerModel)){
#		smallerModel = 0
#	}
	if(split){
		previousclust <- apply(smallerModel$posteriorProbabilities, 1, which.max)
		clust <- previousclust
		index1 <- 1:n
		cc <- numeric(K-1)
		for(k in 1:(K-1)) cc[k] <- length(which(clust==k))
		full.components <- which(cc>1)
		m <- length(full.components)
	}

#	for(iter in 2:tsplit){
	parLoop <- foreach(iter = 1:tsplit, .export = ls(envir = globalenv())) %dopar% 
	{
		if(split){
			w <- matrix(0, nrow = n, ncol = K)
			split.cluster <- sample(full.components, 1)
			#cat(paste0('***   splitting component ', split.cluster,'   ***'),'\n')
			index1 <- which(clust == split.cluster)
			w[ ,1:(K-1)] <- smallerModel$posteriorProbabilities
			u.numbers <- runif(length(index1))
#			u.numbers <- rbeta(length(index1), shape1 = 0.5, shape2=0.5)
			w[index1, K] <- w[index1, split.cluster]*u.numbers
			w[index1, split.cluster] <- w[index1, split.cluster]*(1-u.numbers)
#			proposing values for starting the NR iterates
			# OPTION A: keep the same as the split component
			b[,,-K] = smallerModel$b
			b[,,K] = smallerModel$b[,,split.cluster] 

			# OPTION B: run NR inside each one of them
	#		wm <- apply(w[index1,c(split.cluster,K)],1,which.max)
	#		group1 <- index1[which(wm == 1)]
	#		group2 <- index1[-group1]
	#		if(length(group1)>1){
	#			b[,,split.cluster] <- multinomial_logistic_EM(y = y[group1,], x = x[group1,], K = 1, 
	#				maxIter = maxIter, 
	#				emthreshold = emthreshold, maxNR=maxNR, verbose = F)$beta
	#		}
	#		if(length(group2)>1){
	#			b[,,K] <- multinomial_logistic_EM(y = y[group2,], x = x[group2,], K = 1, maxIter = maxIter, 
	#				emthreshold = emthreshold, maxNR=maxNR, verbose = F)$beta
	#		}
			smallEM <- multinomial_logistic_EM(y=y, x=x, K=K, w_start = w, b_start = b, 
					maxIter = maxIter, emthreshold = emthreshold, maxNR=maxNR, R0 = R0, method = method)
		}else{
			smallEM <- multinomial_logistic_EM(y=y, x=x, K=K, 
					maxIter = maxIter, emthreshold = emthreshold, maxNR=maxNR, R0 = R0, method = method)			
		}
		smallEM
	}
	stopImplicitCluster()
	maxIndex <- which.max(unlist(lapply(parLoop, function(x)tail(x$logLikelihood,1))))
	bestRun <- parLoop[[maxIndex]]
	return(bestRun)
}


shakeEM_GLM <- function(y, x, K, equalModel, tsplit = 10, maxIter = 20, emthreshold = 1e-8, maxNR=5, nCores, split = TRUE, R0, method){
	registerDoParallel(cores = nCores)
	n <- dim(y)[1]
	D <- dim(y)[2]
	p <- dim(x)[2]
	pr <- numeric(K)
	b <- array( data = 0, dim = c(D-1, p, K) )
	theta <- array(data = NA, dim = c(n, D, K))
	epsilon <- 1e-10
	thresh <- -744
	loglik <- numeric(maxIter)
	if(split){
		previousclust <- apply(equalModel$posteriorProbabilities, 1, which.max)
		clust <- previousclust
		index1 <- 1:n
		cc <- numeric(K)
		for(k in 1:K) cc[k] <- length(which(clust==k))
		full.components <- which(cc>1)
		m <- length(full.components)
	}

#	for(iter in 2:tsplit){
	parLoop <- foreach(iter = 1:tsplit, .export = ls(envir = globalenv())) %dopar% 
	{
		if(length(full.components) > 1){
			w <- matrix(0, nrow = n, ncol = K)
			clusters <- sample(full.components, 2)
			#cat(paste0('***   splitting component ', split.cluster,'   ***'),'\n')
			index1 <- which(clust == clusters[1])
			index2 <- which(clust == clusters[2])
			index12 <- sort(c(index1, index2))
			index12_EM <- multinomial_logistic_EM(y = y[index12, ], x = x[index12, ], K = 2, 
					maxIter = maxIter, 
					emthreshold = emthreshold,
					maxNR = maxNR, nCores = 1,
					verbose = F, R0 = R0, method = method)
			w[ ,1:K] <- equalModel$posteriorProbabilities
			u.numbers <- index12_EM$posteriorProbabilities[,1]

			z <- rowSums(w[index12, clusters])
			w[index12, clusters[1]] <- u.numbers * z
			w[index12, clusters[2]] <- (1-u.numbers) * z

			
#			w[index12, clusters[1]] <- w[index12, clusters[1]] * u.numbers
#			w[index12, clusters[2]] <- w[index12, clusters[2]] * (1-u.numbers)

			b = equalModel$b
			b[,,clusters[1]] = index12_EM$b[,,1]
			b[,,clusters[2]] = index12_EM$b[,,2]  
			smallEM <- multinomial_logistic_EM(y = y, x = x, K = K, w_start = w, b_start = b, 
					maxIter = maxIter, emthreshold = emthreshold, maxNR=maxNR, R0 = R0, method = method)
			#cat(paste0('shake: (', clusters[1], ', ', clusters[2],') LogL = ', 
			#	tail(round(smallEM$logLikelihood,2),1)),'\n')
		}else{
			smallEM <- multinomial_logistic_EM(y=y, x=x, K=K, 
					maxIter = maxIter, emthreshold = emthreshold, maxNR=maxNR, R0 = R0, method = method)			
		}
		smallEM
	}
	stopImplicitCluster()
	maxIndex <- which.max(unlist(lapply(parLoop, function(x)tail(x$logLikelihood,1))))
	bestRun <- parLoop[[maxIndex]]
	return(bestRun)
}









mix_mnm_logistic <- function(y, X, Kmax = 10, maxIter = 100, emthreshold = 1e-8, maxNR = 5, nCores, tsplit=8, msplit=5, split = TRUE, shake = TRUE, random = TRUE, criterion = 'ICL', plotting = FALSE, R0 = 0.1, method = 5){
	s <- rowSums(y)
	p <- dim(X)[2]
	n <- dim(y)[1]
	D <- dim(y)[2]
	thresh <- -744
#	set.seed(1)
	runEM <- vector('list',length = Kmax)
	cat(paste0('Running EM for K = ',1) ,'\n')
	runEM[[1]] <- multinomial_logistic_EM(y=y, x = X, K = 1, maxIter = maxIter, emthreshold = emthreshold, maxNR=3*maxNR, verbose = TRUE, R0 = R0, method = method)
	bic <- iclbic <- numeric(Kmax)
	bic[1] <- iclbic[1] <- -2*runEM[[1]]$logLikelihood + p*(D-1)*log(n) 
	
       oldpar <- par(no.readonly = TRUE)
       on.exit(par(oldpar)) 
	
	
	for( K in 2:Kmax){
		if(split){
			cat(paste0('SPLIT Small EM initialization for K = ',K, '...') )
			findStartValues <- splitEM_GLM(y=y, x=X, K = K, smallerModel = runEM[[K-1]], 
					tsplit = tsplit, maxIter = msplit, emthreshold = emthreshold, 
					maxNR = maxNR, nCores = nCores, R0 = R0, method = method)
			if(shake){
				cat(paste0(' log-likelihood = ', round(tail(findStartValues$logLikelihood,1),2)), '\n' )
				cat(paste0('SHAKE Small EM initialization for K = ',K, '...') )
				shakeEM <- shakeEM_GLM(y = y, x = X, K = K, equalModel = findStartValues, 
				tsplit = tsplit, emthreshold=emthreshold, nCores = nCores, 
				maxIter = msplit, maxNR = maxNR, R0 = R0, method = method)
				splitL <- tail(findStartValues$logLikelihood,1)
				shakeL <- tail(shakeEM$logLikelihood,1)
				if(is.numeric(shakeL) == F){shakeL <- splitL - 1}
				if( shakeL > splitL ){
					findStartValues <- shakeEM
				}else{
					cat(paste0('shake aborted'),'\n')
				}
			}
			if(random){
				cat(paste0('RANDOM Small EM initialization for K = ',K, '...') )
				findStartValues2 <- splitEM_GLM(y=y, x=X, K = K,  
						tsplit = tsplit, maxIter = msplit, emthreshold = emthreshold, 
						maxNR = maxNR, nCores = nCores, split = FALSE, R0 = R0, method = method)
				l1 <- tail(findStartValues$logLikelihood,1)
				l2 <- tail(findStartValues2$logLikelihood,1)		
				if(l2 > l1){
					findStartValues <- findStartValues2
					cat(paste0('random starts accepted'),'\n')
				}
			
			}
		}else{
			cat(paste0('RANDOM Small EM initialization for K = ',K, '...') )
			findStartValues <- splitEM_GLM(y=y, x=X, K = K,  
					tsplit = tsplit, maxIter = msplit, emthreshold = emthreshold, 
					maxNR = maxNR, nCores = nCores, split = FALSE, R0 = R0, method = method)			
		}
		cat(paste0(' log-likelihood = ', round(tail(findStartValues$logLikelihood,1),2)), '\n' )
		cat(paste0('Running EM for K = ',K) ,'\n')
		runEM[[K]] <- multinomial_logistic_EM(y = y, x = X, K = K, 
				w_start = findStartValues$posteriorProbabilities,
				b_start =  findStartValues$beta,
				maxIter = maxIter, 
				emthreshold = emthreshold,
				maxNR = maxNR, nCores = nCores,
				verbose = TRUE, R0 = R0, method = method)

		if(FALSE){
			cat(paste0('     SHAKE Small EM initialization for K = ',K, '...') )
			shakeEM <- shakeEM_GLM(y = y, x = X, K = K, equalModel = runEM[[K]], 
				tsplit = tsplit, emthreshold=emthreshold, nCores = nCores, 
				maxIter = msplit, maxNR = maxNR, R0 = R0, method = method)
			findStartValues <- shakeEM
			cat(paste0(' log-likelihood = ', round(tail(findStartValues$logLikelihood,1),2)), '\n' )
			cat(paste0('     RE-Running EM for K = ',K) ,'\n')
			runEM[[K]] <- multinomial_logistic_EM(y = y, x = X, K = K, 
				w_start = findStartValues$posteriorProbabilities,
				b_start =  findStartValues$beta,
				maxIter = maxIter, 
				emthreshold = emthreshold,
				maxNR = maxNR, nCores = nCores,
				verbose = TRUE, R0 = R0, method = method)

		}

		cat('\n')
		bic[K] <- -2*tail(runEM[[K]]$logLikelihood, 1) + (K - 1 + K * (D - 1) * p)*log(n) 
		nz <- runEM[[K]]$posteriorProbabilities
		ind <- 1:K
		for (i in 1:n) {
			index <- ind[nz[i, ] < exp(thresh)]
			nz[i, index] <- rep(exp(thresh), length(index))
			nz[i, ] <- nz[i, ]/sum(nz[i, ])
		}
		iclbic[K] <- bic[K] - 2 * sum(nz * log(nz))
		if(plotting){
			matplot( cbind(bic[1:K], iclbic[1:K]), type='b', xlab = 'K', 
				ylab='criterion',lwd=1.5, col=c('blue', 'darkorange'),pch = c(0,1) )
			points(which.min(bic[1:K]),min(bic[1:K]),col='blue',pch=16)
			points(which.min(iclbic[1:K]),min(iclbic[1:K]),col='darkorange',pch=c(15, 16))
			legend('topright',c('BIC','ICL-BIC'),col = c('blue','darkorange'),lty = 1)
		}
	}

	if(criterion == 'ICL'){	
		selectedK <- which.min(iclbic[1:K])
	}else{
		selectedK <- which.min(bic[1:K])
	}
	zEstimate <- apply(runEM[[selectedK]]$posteriorProbabilities, 1, which.max)
	

	results <- vector('list', length = 5)
	results[[1]] <- selectedK
	results[[2]] <- runEM
	results[[3]] <- bic
	results[[4]] <- iclbic
	results[[5]] <- zEstimate
	names(results) <- c('estimated_K', 'all_runs', 'BIC_values', 'ICL_BIC_values', "estimated_clustering")
	cat('\n')

	if(criterion == 'ICL'){	
		cat(paste0('According to ICL-BIC, the estimated number of components is K = ', selectedK,'.'),'\n' )
		cat(paste0('The ICL-BIC for the selected model is ', round(iclbic[results[[1]]],2),'.'),'\n' )
	}else{
		cat(paste0('According to BIC, the estimated number of components is K = ', selectedK,'.'),'\n' )
		cat(paste0('The BIC for the selected model is ', round(bic[results[[1]]],2),'.'),'\n' )
	}
	
	cat('\n')
	return(results)
}




#######################################################################################################################
#####################################################################################################################
## functions for mcmc
######################################################################################################################
####################################################################################################################



myDirichlet <- function (alpha){
	k <- length(alpha)
	theta <- rgamma(k, shape = alpha, rate = 1)
	return(theta/sum(theta))
}

log_dirichlet_pdf <- function (alpha, weights){
    min_weight <- 1e-300
    weights[weights < min_weight] <- min_weight	
    weights <- weights/sum(weights)
    normConstant <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    pdf <- sum((alpha - 1) * log(weights)) - normConstant
    return(pdf)
}

# b	current value of regression coef
# z	sim allocations
# tau	step >0
# A	optional scaling positive definite matrix (inverse Hessian)
# pr	mixing proportions
# nu2	prior variance for betas

# note that mala_proposal() is not used
# the mala_proposal_cpp() is used instead.
mala_proposal <- function(y, X, b, z, tau, A = FALSE, pr, nu2){
        n <- length(z)
        p <- dim(b)[2]
        K <- dim(b)[3]
        D <- dim(b)[1]+1
        G <- array(data = 0, dim = c(n, D - 1, K))
        s <- rowSums(y)
        normConst <- lgamma(s + 1) - rowSums(lgamma(y + 1))
        FirstDerivative <- array(data = NA, dim = c(p, D-1) )
        log_posterior_gradient <-  array(data = NA, dim = c(D-1, p, K) )
        SecondDerivative <- array(data = NA, dim = c((D-1)*p, (D-1)*p) )
        theta <- theta_prop <- array(data = 0, dim = c(n, D, K))
	sum1 <- matrix(0, nrow = n, ncol=D)

	w <- matrix(0, n, K)
	for(i in 1:n){
		w[i, z[i]] = 1
	}
	cll <- 0
	for(k in 1:K){
		for(d in 1:(D-1)){
		        G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
		        sum1[,d+1] <- G[,d,k]
		}
		lastCol <- rowLogSumExps(sum1)
#		print(lastCol[1])
		G[,,k] <- G[,,k] - lastCol
	#                       G <- G/sum1
		theta[ , 1:(D-1), k] <- G[,,k]
		theta[ , D, k] <- -lastCol
		#bLongOld <- c(bLongOld, b[,,k])
#		cll <- cll + w[,k] * ( log(pr[k]) + rowSums(y * theta[,,k]) + normConst)
		cll <- cll + w[,k] * ( log(pr[k]) + rowSums(y * theta[,,k]))
		G[,,k] <- exp(G[,,k])
		motorhead <- matrix(rep(w[,k], D-1), ncol = D-1) * (y[,-D] - s*G[,,k])
#		print(paste0('k = ', k))
#		print(head(motorhead))
		#protector <- w[ ,k] * G[ , , k]
		for(r in 1:p){
			FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
		}
		log_posterior_gradient[,,k] <- t(FirstDerivative) - 	# log-likelihood term
						b[,,k]/nu2	# prior term
	}
	results <- vector('list', length = 4)
	results[[4]] <- log_posterior_gradient	
	
#	print(G[1,,])
#	print(theta[1,,])
#	print(paste0('cll = ', sum(cll)))
	npars <- (D-1)*p*K
	if(A == FALSE){
		A = diag(npars)
	}else{
#		inv Hessian
	}

	norm_rv <- array(rnorm(npars), dim = c(D - 1, p, K))


#	prop_mean <- b + tau * A %*% log_posterior_gradient
	prop_mean <- b + tau * log_posterior_gradient
#	prop_var <- 2*tau*A
	prop_var <- 2*tau
	b_prop <- prop_mean + sqrt(prop_var) * norm_rv


	cll_prop <- 0
	sum1 <- matrix(0, nrow = n, ncol=D)
	for(k in 1:K){
		for(d in 1:(D-1)){
		        G[,d,k] <- apply(X,1,function(x)t(b_prop[d,,k]) %*% x)
		        sum1[,d+1] <- G[,d,k]
		}
		lastCol <- rowLogSumExps(sum1)
		G[,,k] <- G[,,k] - lastCol
	#                       G <- G/sum1
		theta_prop[ , 1:(D-1), k] <- G[,,k]
		theta_prop[ , D, k] <- -lastCol
		#bLongOld <- c(bLongOld, b[,,k])
#		cll_prop <- cll_prop + w[,k] * ( log(pr[k]) + rowSums(y * theta_prop[,,k]) + normConst)
		cll_prop <- cll_prop + w[,k] * ( log(pr[k]) + rowSums(y * theta_prop[,,k]))
		G[,,k] <- exp(G[,,k])
		motorhead <- matrix(rep(w[,k], D-1), ncol = D-1) * (y[,-D] - s*G[,,k])
		#protector <- w[ ,k] * G[ , , k]
		for(r in 1:p){
			FirstDerivative[r, ] <- t(motorhead) %*% X[,r]
		}
		log_posterior_gradient[,,k] <- t(FirstDerivative) - 	# log-likelihood term
						b_prop[,,k]/nu2	# prior term
	}
	old_mean <- b_prop + tau * log_posterior_gradient
#	prop_var <- 2*tau*A
	old_var <- 2*tau

#	print(paste0('cll_prop = ', sum(cll_prop)))
	
	log_priorRatio <- sum(dnorm(b_prop, mean = 0, sd = sqrt(nu2), log = TRUE)) - 
		sum(dnorm(b, mean = 0, sd = sqrt(nu2), log = TRUE))
	log_proposalRatio <- sum(dnorm(b, old_mean, sqrt(old_var), log = TRUE)) - 
				sum(dnorm(b_prop, prop_mean, sqrt(prop_var), log = TRUE))
	lar <- sum(cll_prop) - sum(cll) + log_priorRatio + log_proposalRatio
#	print(paste0('log_priorRatio = ', log_priorRatio))
#	print(paste0('log_proposalRatio = ', log_proposalRatio))
#	print(paste0('lar = ', lar))

	acceptance = 0
	if(log(runif(1)) < lar){
		theta = theta_prop
		b = b_prop
		acceptance = 1
	}

	results[[1]] <- theta
	results[[2]] <- b
	results[[3]] <- acceptance

	names(results) <- c('theta', 'b', 'acceptance', 'gradient')
	return(results)
}



gibbs_mala_sampler <- function(y, X, tau = 0.00003, nu2,  K, mcmc_iter = 100, alpha_prior=NULL, start_values = 'EM', em_iter = 10, thin = 10, verbose = FALSE, checkAR = NULL, probsSave = FALSE, ar_low = 0.4, ar_up = 0.6){
	n <- dim(y)[1]
	D <- dim(y)[2]
	p <- dim(X)[2]
	b <-  array( data = rnorm((D-1)*p*K, mean = 0, sd = 0.1), dim = c(D-1, p, K) ) #array( data = 0, dim = c(D-1, p, K) )
	theta <- array(data = NA, dim = c(n, D, K))
	s <- rowSums(y)
	log_multinom_coefficient <-  lgamma(s+1) - rowSums(lgamma(y + 1))
	G <- array(data = 0, dim = c(n, D - 1, K))
	sum1 <- matrix(0, nrow = n, ncol = D)
	cluster_size <- numeric(K)
	nSaved <- floor(mcmc_iter/thin)
	kValues <- numeric(nSaved)
	zValues <- matrix(NA, nSaved, n)
	prValues <- matrix(NA, nSaved, K)
	bValues <- array(NA, dim = c(nSaved, D-1, p, K))
	if(is.null(alpha_prior)){
		#alpha_prior <- rep(1,K)
		alpha_prior <- rep(1, K)/K
	}
	mh_acceptance_rate <- window_accept <- window_accept_rate <- 0
	llValues <- completeLL <- numeric(nSaved)
	if(probsSave){
	class_probs <- array(data = 0, dim = c(nSaved, n, K))
	}
	
       oldpar <- par(no.readonly = TRUE)
       on.exit(par(oldpar)) 

	
	if(start_values[1] == 'EM'){	
#		random EM initialization
		cat(paste0('    EM initialization...'))
		mnm_em <- multinomial_logistic_EM(y = y, x = X, K = K, maxIter = em_iter, 
			emthreshold = 1e-8, maxNR = 5, nCores = NULL, verbose = FALSE, R0 = 0.1, method = 5)
		cat(paste0('    OK'), '\n')
#		mnm_em <- splitEM_GLM(y=y, x=X, K = K,  
#				tsplit = em_starts, maxIter = em_iter, emthreshold = 1e-8, 
#				maxNR = 5, nCores = 8, split = FALSE, R0 = 0.1, method = 5)                

		pr <- mnm_em$weights
		b <- mnm_em$beta
		z <- apply(mnm_em$posteriorProbabilities, 1, which.max)

	}else if(start_values[1] == 'RANDOM'){
#		random start	
		pr <- myDirichlet(rep(1,K))
		z <- sample(1:K, n, prob = pr, replace = TRUE)
		b <- array(data = -0.2*runif((D-1)*p*K) + 0.1, dim = c(D-1,p,K))
	}else{
		if (sum(names(start_values)  == c('b', 'z','pr')) == 3){
			pr <- start_values$pr
			z <- start_values$z
			b <- start_values$b
		}else{
			stop('start_values should be a list with entries `b`, `z` and `pr` corresponding to regression coefficients, allocations and mixing proportions, respectively')		
		}
	}


	for(k in 1:K){
		for(d in 1:(D-1)){
		        G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
		        sum1[,d+1] <- G[,d,k]
		}
		lastCol <- rowLogSumExps(sum1)
		G[,,k] <- G[,,k] - lastCol
	#                       G <- G/sum1
		theta[ , 1:(D-1), k] <- G[,,k]
		theta[ , D, k] <- -lastCol
	}       
	
	llValues[1] <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
	kValues[1] <- length(table(z))
	zValues[1,] <- z
	probs <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$log_weighted_density_k
	probs_log <- probs
	probs <- array(t(apply(probs, 1, function(tmp) {return(exp(tmp - max(tmp)))})), dim = c(n, K))
	completeLL[1] <- 0
	for (k in 1:K) {
	    index <- which(z == k)
	    completeLL[1] <- completeLL[1] + sum(probs_log[index,k])
	}
	if(probsSave){
		class_probs[1,,] <- probs 
	}
#------------
	min_weight <- 1e-300
	for(iter in 2:mcmc_iter){

#		metropolis_step <- mala_proposal(y = y, X = X, b = b, z = z, tau = tau, A = FALSE, pr = pr, nu2 = nu2)
		metropolis_step <- mala_proposal_cpp(y = y, X = X, b = b, z = z, tau = tau, pr = pr, nu2 = nu2)
		theta <- metropolis_step$theta
		if(metropolis_step$acceptance == 1){
			b <- metropolis_step$b
			mh_acceptance_rate <- mh_acceptance_rate + 1
			window_accept <- window_accept + 1
		}

#		for(k in 1:K){
#			for(d in 1:(D-1)){
#				G[,d,k] <- apply(X,1,function(x)t(b[d,,k]) %*% x)
#				sum1[,d+1] <- G[,d,k]
#			}
#			lastCol <- rowLogSumExps(sum1)
#			G[,,k] <- G[,,k] - lastCol
#			theta[ , 1:(D-1), k] <- G[,,k]
#			theta[ , D, k] <- -lastCol
#		}       

	#	ta probs einai ta log p + log f_k
		ll <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)
		probs <- ll$log_weighted_density_k
		probs_log <- probs
		probs <- array(t(apply(probs, 1, function(tmp) {return(exp(tmp - max(tmp)))})), dim = c(n, K))

		z <- apply(probs, 1, function(tmp) {
				if (anyNA(tmp)) {
				    tmp <- rep(1, K)
				}
				return(sample(K, 1, prob = tmp))
			}
		)
		for (k in 1:K) {
		    index <- which(z == k)
		    cluster_size[k] <- length(index)
		}


		pr <- myDirichlet(alpha_prior + cluster_size)


		pr[pr < min_weight] <- min_weight	
		pr <- pr/sum(pr)



		if(iter%%thin == 0){
			kValues[iter/thin] <- length(table(z))
			zValues[iter/thin,] <- z

			ll <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)
			probs <- ll$log_weighted_density_k
			probs_log <- probs
			probs <- array(t(apply(probs, 1, function(tmp) {return(exp(tmp - max(tmp)))})), dim = c(n, K))
			completeLL[iter/thin] <- 0
			for (k in 1:K) {
			    index <- which(z == k)
			    completeLL[iter/thin] <- completeLL[iter/thin] + sum(probs_log[index,k])
			}

			llValues[iter/thin] <- ll$ll
			#llValues[iter/thin] <- mixLoglikelihood_GLM(y = y, theta = theta, pi = pr)$ll
			prValues[iter/thin,] <- pr		
			bValues[iter/thin,,,] <- b
			if(probsSave){
				class_probs[iter/thin,,] <- probs 
			}

		}
		
		if(is.null(checkAR)==FALSE){
			if(iter %% checkAR == 0){
				window_accept_rate <- window_accept/checkAR
				if(window_accept_rate < ar_low){
					tau <- 0.9 * tau
				}
				if(window_accept_rate > ar_up){
					tau <- tau/0.9
				}
				
				window_accept <- 0
			}
		}

		if(verbose){
			if(iter %% checkAR == 0){
				par(mfrow = c(1,2))
				plot(llValues[1:(iter/thin)], type = 'l')
				plot(kValues[1:(iter/thin)], type = 'l')
				cat(paste0('iter = ', iter, ', logL = ', round(llValues[iter/thin],2), 
					', ar = ', 100*mh_acceptance_rate/iter, ', war = ', window_accept_rate, 
						', tau = ', tau), '\n')
			}
		}



	}
	results <- vector('list', length = 9)
	results[[1]] <- kValues
	results[[2]] <- zValues
	results[[3]] <- llValues
	results[[4]] <- prValues
	results[[5]] <- bValues
	results[[6]] <- tau
	results[[7]] <- completeLL	
	results[[8]] <- NA
	results[[9]] <- 100*mh_acceptance_rate/iter
	if(probsSave){
		results[[8]] <- class_probs
	}
	names(results) <- c('nClusters', 'allocations', 'logLikelihood', 'mixing_proportions', 'coefficients', 'tau', 'complete_logLikelihood', 'class_probs', "AR")
	return(results)
}


# prior parallel tempering

gibbs_mala_sampler_ppt <- function(
	y, X, tau = 0.00003, nu2,  K, 
	mcmc_cycles = 100, 
	iter_per_cycle = 10,
	dirPriorAlphas,
	start_values = 'EM', 
	em_iter = 10,  nChains = 4, 
	nCores = 4, warm_up = 100, checkAR = 50, probsSave = FALSE, showGraph = 50, ar_low = 0.4, ar_up = 0.6, withRandom = TRUE){
	
       oldpar <- par(no.readonly = TRUE)
       on.exit(par(oldpar)) 
 
	
	weights <- array(data = NA, dim = c(2, K))
	mh_acceptance_rate <- 0
	n <- dim(y)[1]
	D <- dim(y)[2]
	p <- dim(X)[2]
	b <-  array( data = rnorm((D-1)*p*K, mean = 0, sd = 0.1), dim = c(D-1, p, K) ) #array( data = 0, dim = c(D-1, p, K) )
	kValues <- matrix(NA, mcmc_cycles, nChains)
	zValues <- matrix(NA, mcmc_cycles, n)
	prValues <- matrix(NA, mcmc_cycles, K)
	bValues <- array(NA, dim = c(mcmc_cycles, D-1, p, K))
	llValues <- completeLL <- numeric(mcmc_cycles)
	if(probsSave){
		class_probs <- array(data = 0, dim = c(mcmc_cycles, n, K))
	}

	

	if (missing(dirPriorAlphas)) {
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN * (2:nChains - 1))/(1*K)
	}

	registerDoParallel(nCores)

	if(is.list(start_values) == FALSE){withRandom = FALSE}
	cat(paste0('Warm-up period...'), '\n')
	parLoop <- foreach(iter = 1:nChains, .export = ls(envir = globalenv())) %dopar% {
		#alpha_prior = rep(p*(D-1), K)
		alpha_prior = rep(dirPriorAlphas[iter], K)


		sv = start_values
		if( withRandom ){
			#if(iter %% 2 == 0){
			#	sv = 'RANDOM'
			#}
			if(iter  > 1){
				# this draws a random permutation of the starting values
				myPerm <- sample(K, K, replace = FALSE)
				#print(myPerm)
				sv$b <- start_values$b[,,myPerm] 
				#print(sv$b)
				sv$pr <- start_values$pr[myPerm]
				#print(sv$pr)
				invPerm <- order(myPerm)
				#print(invPerm)
				sv$z <- invPerm[start_values$z]
				#print(sv$z)
			}

		}
		
		mcmc_chain <- gibbs_mala_sampler(y = y, X = X, alpha_prior = alpha_prior, 
				tau = tau, nu2 = nu2,  K = K, mcmc_iter = warm_up, 
				start_values = sv, em_iter = em_iter, thin = warm_up, checkAR = checkAR, probsSave = probsSave, ar_low = ar_low, ar_up = ar_up)


	}
	stopImplicitCluster()
	cat(paste0('OK'), '\n')
	tauM <- numeric(nChains)
	for(iter in 1:nChains){	
		tauM[iter] <- parLoop[[iter]]$tau
		cat(paste0('tau = ', tauM[iter]),"\n")
	}
#	tauF <- median(tau)
#	for(iter in 1:nChains){	
#		tau[iter] <- tauF
#	}
		
	
	for(mcmc_iter in 1:mcmc_cycles){

		registerDoParallel(nCores)
		parLoop2 <- foreach(iter = 1:nChains, .export = ls(envir = globalenv())) %dopar% {
			start_values <- vector('list', length = 3)
			start_values$b <- parLoop[[iter]]$coefficients[1,,,]
			start_values$z <- parLoop[[iter]]$allocations[1,]
			start_values$pr <- parLoop[[iter]]$mixing_proportions[1,]
		
			alpha_prior = rep(dirPriorAlphas[iter], K)
			mcmc_chain <- gibbs_mala_sampler(y = y, X = X, alpha_prior = alpha_prior, 
#					tau = parLoop[[iter]]$tau, 
					tau = tauM[iter],
					nu2 = nu2,  K = K, mcmc_iter = iter_per_cycle, 
					start_values = start_values, thin = iter_per_cycle, checkAR = NULL, probsSave = probsSave)
		}
		stopImplicitCluster()

		if (nChains > 1) {
			chains <- sample(nChains - 1, 1)
			chains <- c(chains, chains + 1)
			weights[1, ] <- parLoop2[[chains[1]]]$mixing_proportions[1,]
			weights[2, ] <- parLoop2[[chains[2]]]$mixing_proportions[1,]
			mh_denom <- log_dirichlet_pdf(rep(dirPriorAlphas[chains[1]], K), weights[1, ]) + 
				log_dirichlet_pdf(rep(dirPriorAlphas[chains[2]], K), weights[2, ])
		    	mh_nom <- log_dirichlet_pdf(rep(dirPriorAlphas[chains[2]], K), weights[1, ]) +
		    			log_dirichlet_pdf(rep(dirPriorAlphas[chains[1]], K), weights[2, ])
		    	mh_ratio <- mh_nom - mh_denom
			if (log(runif(1)) < mh_ratio) {
				tmp <- parLoop2[[chains[2]]]
				parLoop2[[chains[2]]] <- parLoop2[[chains[1]]]
				parLoop2[[chains[1]]] <- tmp
				mh_acceptance_rate <- mh_acceptance_rate + 1
			}
		}



		parLoop <- parLoop2
		kValues[mcmc_iter,] <- unlist(lapply( parLoop, function(x)x$nClusters ))
		zValues[mcmc_iter,] <- parLoop[[1]]$allocations[1,]
		prValues[mcmc_iter,] <- parLoop[[1]]$mixing_proportions[1,]
		bValues[mcmc_iter,,,] <- parLoop[[1]]$coefficients[1,,,]
		llValues[mcmc_iter] <- parLoop[[1]]$logLikelihood
		completeLL[mcmc_iter] <- parLoop[[1]]$complete_logLikelihood
		if(probsSave){
			class_probs[mcmc_iter, , ] <- parLoop[[1]]$class_probs
		}
		if(mcmc_iter %% showGraph == 0){
			par(mfrow = c(1,3))
			matplot(kValues[1:mcmc_iter,], type = 'l', col = 1:nChains, lwd = c(3,rep(1,nChains -1)))
			plot(llValues[1:mcmc_iter], type = 'l')
			theta <- by(y, zValues[mcmc_iter,], FUN = function(x)colMeans(x/rowSums(x) ), simplify=FALSE)
			plot(theta[[1]], type = 'l', col = 1, lwd = 2, ylim = c(0,1))
			for(k in 1:length(table(zValues[mcmc_iter,]))){
				points(theta[[k]], type = 'l', col = k, lwd = 2)
			}
			#plot(completeLL[1:mcmc_iter], type = 'l')

			print(table(zValues[mcmc_iter,]))
			cat(paste0('tau = ', parLoop2[[1]]$tau, ", proposal AR = ", parLoop[[1]]$AR,'%, mh_ar = ', round(100*mh_acceptance_rate/mcmc_iter,2),'%.'),'\n')
			
		}
	}	

	results <- vector('list', length = 8)
	results[[1]] <- kValues
	results[[2]] <- zValues
	results[[3]] <- llValues
	results[[4]] <- prValues
	results[[5]] <- bValues
	results[[6]] <- mh_acceptance_rate
	results[[7]] <- completeLL
	results[[8]] <- NA
	if(probsSave){
		results[[8]] <- class_probs
	}
	names(results) <- c('nClusters', 'allocations', 
			'logLikelihood', 'mixing_proportions', 'coefficients', 'mh_acceptance_rate', 'complete_logLikelihood', 'class_probs')
	return(results)


}


dealWithLabelSwitching <- function(gs, burn, thin = 10, zPivot = NULL, returnRaw = FALSE, maxM = NULL){

	if(is.null(maxM)){
        	m <- dim(gs$nClusters)[1]
        }else{
        	m <- maxM
        }
        gs$nClusters <- gs$nClusters[seq(1,m, by = thin),]
        gs$allocations <- gs$allocations[seq(1,m, by = thin),]
        gs$logLikelihood <- gs$logLikelihood[seq(1,m, by = thin)]
        gs$mixing_proportions <- gs$mixing_proportions[seq(1,m, by = thin),]
        gs$coefficients <- gs$coefficients[seq(1,m, by = thin),,,]
        gs$complete_logLikelihood <-    gs$complete_logLikelihood[seq(1,m, by = thin)]
        
        
        kTab <- table(gs$nClusters[-(1:burn),1])
        k_posterior <- kTab/(length(gs$nClusters[,1])-burn)
        K_est = as.numeric(names(kTab)[which.max(k_posterior)])
        n <- dim(gs$allocations)[2]
        index <- which(gs$nClusters[,1] == K_est)
        index <- index[index>burn]
        Kindex <- index
        logl <- gs$logLikelihood[Kindex]
        z <- gs$allocations[Kindex, ]

#        zPivot <-z[which.max(logl),]

        K_max <- max(z)
        K_max2 <- dim(gs$mixing_proportions)[2]
        p <- dim(gs$coefficient)[3]
        theta <- vector('list', length = p)
        thetaRetained <- vector('list', length = p)
        
        
	newDim <- dim(gs$coefficients)
	newDim[c(1,4)] <- c(length(Kindex), K_est)
	retainedCoefficients <- array(data = NA, dim = newDim)

	m <- 0
	for(i in Kindex){
		m <- m + 1
		alive <- sort(unique(z[m,]))
		for(j in 1:p){
			retainedCoefficients[m, ,j, ] <-  gs$coefficients[i,,j,alive]
		}
	}

	z <- t(apply(z, 1, function(y)as.numeric(as.factor(y))))
	K_max <- max(z)
	if(is.null(zPivot)==TRUE){
	zPivot <-z[which.max(logl),]
	}else{
		zPivot <- as.numeric(as.factor(zPivot))
	}
        




    if(K_est > 1){
                ls <- label.switching(method = 'ECR', K = K_max, zpivot = zPivot, z = z)

                cluster_estimate <- ls$clusters[1,]

                m <- length(Kindex)
                allocationsECR <- matrix(NA, m, n)
                for (i in 1:m) {
                    myPerm <- order(ls$permutations[['ECR']][i, ])
                    allocationsECR[i, ] <- myPerm[z[i, ]]
                }
                for(j in 1:p){
                        tmp <- array(data = NA, dim = c(length(Kindex), K_max, dim(gs$coefficients)[2]))
                        for(k in 1:K_max){
                                tmp[,k,] <- retainedCoefficients[,,j,k]
                        }
                        theta[[j]] <- permute.mcmc( tmp, ls$permutations[['ECR']] )$output
                }
                
        }else{
                cluster_estimate  <- zPivot
                allocationsECR <- zPivot
        }
        results <- vector('list', length = 5)
        results[[1]] <- cluster_estimate
        results[[2]] <- k_posterior

        mcmc <- vector('list', length = 2)
        mcmc[[1]] <- allocationsECR
        mcmc[[2]] <- theta
        results[[3]] <- mcmc
        results[[4]] <- matrix(0, n, K_est)
	 if(K_est > 1){
		for(i in 1:n){
			for(k in 1:K_est){
				results[[4]][i, k] <- length( which(allocationsECR[, i] == k) )/m
			}
		}
        }else{
        	results[[4]] <- matrix(1, n, K_est)
        }
        if(returnRaw){
	        results[[5]] <- retainedCoefficients
        }

        names(results) <- c('cluster', 'nClusters_posterior', 'mcmc', 'posteriorProbabilities', 'raw')
        return(results)
}



simulate_multinomial_data <- function(K, p, D, n, size = 20, prob = 0.025, betaTrue = NULL){
	nSize <- n 
	s <- rnbinom(nSize, size = size, prob = prob) 
	sd <- 1 + 4*runif(1)
	if(is.null(betaTrue)){
	betaTrue <- array(data = 0, dim = c(D-1,p,K))
	for(k in 1:K){
		for(d in 1:(D-1)){
			for(r in 1:p){
				u <- runif(1)
				if(u < 0.5){
					betaTrue[d,r,k] <- rnorm(1,mean=0,sd=sd)
	#				betaTrue[j,r,k] <- rt(1,1,0)
				}
			}
		}
	}
	}
	#	design matrix (w const term)
	X <- matrix(data = 0, nrow = n, ncol = p)
	if(p > 1){
		X[,1] <- 1
		Sigma <- matrix(data = 1, nrow = p-1,ncol=p-1)
		diag(Sigma) <- rep(2,p-1)
		X[,-1] <- mvrnorm(n = n, rep(1, p-1), Sigma)
	}
	if(p==1){
	#	design matrix (no const term)
	X <- matrix(data = 0, nrow = n, ncol = p)
	X[,1] <- rnorm(n)
	}

	#	multinomial probabilities (per subject)
	theta <- array(data = 0, dim = c(n, D, K))
	for(i in 1:n){
		for(k in 1:K){
			denom <- 0
			for(d in 1:(D-1)){
				theta[i,d,k] <-  exp(sum(betaTrue[d,,k]*X[i,]))
				denom <- denom + theta[i,d,k]
			}
			denom <- denom + 1
			theta[i,1:(D-1),k] <- theta[i,1:(D-1),k]/denom
			theta[i,D,k] <- 1/denom
		}
	}
	#	data simulation
	y <- matrix(data = NA, nrow = n, ncol = D)
	zTrue <- numeric(n)
	pTrue <- numeric(K)
	pTrue <- 1:K
	pTrue <- pTrue/sum(pTrue)
	zTrue <- sample(K, n, replace = TRUE, prob = pTrue)
	for(k in 1:K){
		ind <- which(zTrue == k)
		for(i in ind){
			y[i, ] <- rmultinom(1,s[i],theta[i,,k])
		}
	}
	ll <- mixLoglikelihood_GLM(y = y, theta = log(theta), pi = pTrue)
	postProb <- round((exp(ll$log_weighted_density_k)/rowSums(exp(ll$log_weighted_density_k))),2)
	mProb <- apply(postProb,1,max)
	zz <- apply(postProb,1,which.max)

	result <- vector("list", length = 4)
	result[[1]] <- y
	result[[2]] <- X
	result[[3]] <- zTrue
	result[[4]] <- zz	
	names(result) <- c("count_data", "design_matrix", "clustering", "z_MAP_true")
	return(result)
}


multinomialLogitMix <- function(response, 
				design_matrix, 
				method, 
				Kmax = 10, 
				mcmc_parameters = NULL, 
				em_parameters = NULL, 
				nCores,
				splitSmallEM = TRUE
				){

	if(missing(nCores)){
		stop("Please supply the number of cores in `nCores`.")
	}
	if(is.integer(response) == FALSE){
		stop("The `response` matrix should consist of integers.")
	}


#	default values for MCMC
		if(is.null(mcmc_parameters)){
			mcmc_parameters <- list()
			if(is.null(mcmc_parameters$tau)){mcmc_parameters$tau <- 0.00035}
			if(is.null(mcmc_parameters$nu2)){mcmc_parameters$nu2 <- 100}
			if(is.null(mcmc_parameters$mcmc_cycles)){mcmc_parameters$mcmc_cycles <- 2600}
			if(is.null(mcmc_parameters$iter_per_cycle)){mcmc_parameters$iter_per_cycle <- 20}	
			if(is.null(mcmc_parameters$nChains)){mcmc_parameters$nChains <- 8}	
			if(is.null(mcmc_parameters$dirPriorAlphas)){
				mcmc_parameters$dirPriorAlphas = c(1,  1 + 
					5 * exp((seq(2, 14, length = mcmc_parameters$nChains - 1)))/100)/(200)}	
			if(is.null(mcmc_parameters$warm_up)){mcmc_parameters$warm_up <- 48000}	
			if(is.null(mcmc_parameters$checkAR)){mcmc_parameters$checkAR <- 500}	
			if(is.null(mcmc_parameters$probsSave )){mcmc_parameters$probsSave <- FALSE}	
			if(is.null(mcmc_parameters$showGraph)){mcmc_parameters$showGraph <- 100}	
			if(is.null(mcmc_parameters$ar_low)){mcmc_parameters$ar_low <- 0.15}	
			if(is.null(mcmc_parameters$ar_up)){mcmc_parameters$ar_up <- 0.25}	
			if(is.null(mcmc_parameters$burn)){mcmc_parameters$burn <- 100}			
			if(is.null(mcmc_parameters$thin)){mcmc_parameters$thin <- 1}				
			if(is.null(mcmc_parameters$withRandom)){mcmc_parameters$withRandom = TRUE}


		}else{
			if(is.null(mcmc_parameters$tau)){mcmc_parameters$tau <- 0.00035}
			if(is.null(mcmc_parameters$nu2)){mcmc_parameters$nu2 <- 100}
			if(is.null(mcmc_parameters$mcmc_cycles)){mcmc_parameters$mcmc_cycles <- 2600}
			if(is.null(mcmc_parameters$iter_per_cycle)){mcmc_parameters$iter_per_cycle <- 20}	
			if(is.null(mcmc_parameters$nChains)){mcmc_parameters$nChains <- 8}	
			if(is.null(mcmc_parameters$dirPriorAlphas)){
				mcmc_parameters$dirPriorAlphas = c(1,  1 + 
					5 * exp((seq(2, 14, length = mcmc_parameters$nChains - 1)))/100)/(200)}	
			if(is.null(mcmc_parameters$warm_up)){mcmc_parameters$warm_up <- 48000}	
			if(is.null(mcmc_parameters$checkAR)){mcmc_parameters$checkAR <- 500}	
			if(is.null(mcmc_parameters$probsSave )){mcmc_parameters$probsSave <- FALSE}	
			if(is.null(mcmc_parameters$showGraph)){mcmc_parameters$showGraph <- 100}	
			if(is.null(mcmc_parameters$ar_low)){mcmc_parameters$ar_low <- 0.15}	
			if(is.null(mcmc_parameters$ar_up)){mcmc_parameters$ar_up <- 0.25}	
			if(is.null(mcmc_parameters$burn)){mcmc_parameters$burn <- 100}			
			if(is.null(mcmc_parameters$thin)){mcmc_parameters$thin <- 1}				
			if(is.null(mcmc_parameters$withRandom)){mcmc_parameters$withRandom = TRUE}
		}
#	default values for EM
	if(is.null(em_parameters)){
		em_parameters <- list(maxIter = 100, emthreshold = 1e-08, 
		    maxNR = 10, tsplit = 16, msplit = 10, split = TRUE, 
		    R0 = 0.1, plotting = TRUE)
		
	}

	if(splitSmallEM){
		em_run <- mix_mnm_logistic(y = response, X = design_matrix, Kmax = Kmax, maxIter = em_parameters$maxIter, 
			emthreshold = em_parameters$emthreshold, maxNR = em_parameters$maxNR, nCores = nCores, tsplit = em_parameters$tsplit, 
			msplit = em_parameters$msplit, split = TRUE, R0 = em_parameters$R0, method = 5, plotting = em_parameters$plotting)
	}else{
		if(method == "EM"){
			em_run <- mix_mnm_logistic(y = response, X = design_matrix, Kmax = Kmax, maxIter = em_parameters$maxIter, 
				emthreshold = em_parameters$emthreshold, maxNR = em_parameters$maxNR, nCores = nCores, tsplit = em_parameters$tsplit, 
				msplit = em_parameters$msplit, split = FALSE, R0 = em_parameters$R0, method = 5, plotting = em_parameters$plotting)
		}else{
		em_run = NA}
	}
	mcmc_run <- postProcessMCMC <- NA
	if(method == "MCMC"){

		D <- dim(response)[2]
		p <- dim(design_matrix)[2]

		if(splitSmallEM){
			start_values <- vector('list', length = 3)
			start_values$b <- array(data = 0, dim = c(D-1, p, Kmax))
			start_values$z <- numeric(dim(response)[1])
			start_values$pr <- numeric(Kmax)

			start_values$b[,,1:em_run$estimated_K] <- em_run$all_runs[[em_run$estimated_K]]$beta
			start_values$z <- em_run$estimated_clustering
			start_values$pr[1:em_run$estimated_K] <- em_run$all_runs[[em_run$estimated_K]]$weights
		
			mcmc_run <-  gibbs_mala_sampler_ppt( y = response, X = design_matrix, tau = mcmc_parameters$tau, 
					nu2 = mcmc_parameters$nu2,  K = Kmax, dirPriorAlphas = mcmc_parameters$dirPriorAlphas, 
					mcmc_cycles = mcmc_parameters$mcmc_cycles, iter_per_cycle = mcmc_parameters$iter_per_cycle, 
					start_values = start_values, nChains = mcmc_parameters$nChains, nCores = nCores, 
					warm_up = mcmc_parameters$warm_up, showGraph = mcmc_parameters$showGraph, checkAR = mcmc_parameters$checkAR, 
					withRandom = mcmc_parameters$withRandom)

		
		
		}else{
			mcmc_run <- gibbs_mala_sampler_ppt(y = response, X = design_matrix, tau = mcmc_parameters$tau, 
					nu2 = mcmc_parameters$nu2,  K = Kmax, dirPriorAlphas = mcmc_parameters$dirPriorAlphas, 
					mcmc_cycles = mcmc_parameters$mcmc_cycles, iter_per_cycle = mcmc_parameters$iter_per_cycle, 
					start_values = "RANDOM", nChains = mcmc_parameters$nChains, nCores = nCores, 
					warm_up = mcmc_parameters$warm_up, showGraph = mcmc_parameters$showGraph, checkAR = mcmc_parameters$checkAR)
		
		}

		postProcessMCMC <- dealWithLabelSwitching(mcmc_run, burn = mcmc_parameters$burn, thin = mcmc_parameters$thin)
	}

	result <- list(
		em_run,
		mcmc_run,
		postProcessMCMC
	)
	names(result) <- c("EM", "MCMC_raw", "MCMC_post_processed")
	return(result)

}








