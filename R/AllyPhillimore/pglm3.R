
library(ape)
library(mvtnorm) # need to download from CRAN
library(MASS)

# ----------------------------------------------------------------------
# This generates random data
#
rand.dat <- function(lambda, Vmat) {
	imat <- matrix(0, nrow = dim(Vmat)[1], ncol = dim(Vmat)[2] )
	diag(imat) <- 1	
	Vmat <- lambda * Vmat + (1 - lambda) * imat
	dat <-  rmvnorm(1, matrix( 0, nrow = dim(Vmat)[1] ), Vmat) 
	dat <- as.matrix(t(dat))
	return(dat)
	}

# ----------------------------------------------------------------------
# This first section of code fits lambda to data on single traits
#

lik.lambda <- function(x, V, lambda) {
	
	lamTrans <- function(Vmat, lambda) {
		V1 <- Vmat
		diag(V1) <- 0
		V2 <- diag( diag(Vmat), ncol = length(Vmat[1,]), nrow = length(Vmat[,1]))
		Vmat <- V1 * lambda + V2
		return(Vmat)
		}
	
	est.mean <- function(y, V) {
		iV <- solve(V, tol = .Machine$double.eps)
		xdum <- matrix(1, nrow = length(y))
		xVix <- crossprod(xdum, iV %*% xdum)
		xViy <- crossprod(xdum, iV %*% y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy # This is a bad thing!
		return(mu[1])
		}

	est.var <- function(x, V) {
		mu <- est.mean(x, V)
		iV <- solve(V, tol = .Machine$double.eps)
		e <- x - mu
		s2 <- crossprod(e, iV%*%e)
		n <- length(x) 
		return(s2 / (n - 1) )
		}

	mv.lik <- function(x, Vmat) {
		mu <- est.mean(x, Vmat)
		s2 <- est.var(x, Vmat)
		n <- length(x)
		logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n-1) / 2.0
		return( list(ll = ll, mu = mu, s2 = s2) )
		}

	Vmat <- lamTrans(V, lambda)
	return(mv.lik(x, Vmat))
	}


max.lik.lambda <- function(x, V) {
		ll <- function(lambda) return(-1* lik.lambda(x, V, lambda)$ll )
		optimize(ll, c(0,1)) 
		}

profile.lambda <- function(x, V) {
		lb <- seq(0, 1, by = 0.01)
		ll <- sapply(lb, function(l) return(lik.lambda(x, V, l )$ll ) )
		plot(lb, ll, type = "l", col =  "red")
		}

# This fits a GLM, correcting for phylogeny, with the option
# of setting the value of lambda, the index of phylogenetic 
# dependence (default is 1.0)
pglm <- function(formula, data, phylomat, lambda = 1.0, ...) {

	prune <- function(dat, Vmat) {
		complete <- complete.cases(dat)
		idx <- which(complete == TRUE)
		Vmat <- Vmat[idx, idx]
		dat <- dat[idx,]

# ------ Delete from here if you want to take the risk! ---------------------	
		# Makes sure data and matrix are in the same order
		nms <- row.names(dat)
		if(length(nms) == 0) stop("Need to supply row names for the data")
		idx <- sort(nms, index.return = TRUE)$ix
		dat <- dat[idx,]
		Vnms <- row.names(Vmat)
		if(length(Vnms) == 0) stop("Need to supply row names for the Variance matrix")
		idx <- sort(Vnms, index.return = TRUE)$ix
		Vmat <- Vmat[idx, idx]
		
		nms <- sort(nms)
		Vn <- sort(row.names(Vmat))
		idx <- which(nms == Vn)
		idx <- which(idx == FALSE)
		
		if(length(idx) > 0) stop("Error, taxon names do not match")
# ------ Delete to here if you want to take the risk! ---------------------


		return(list(dat = dat , Vmat = Vmat))
		}
	
	Dfun <- function(Cmat) {
		iCmat <- solve(Cmat,  tol = .Machine$double.eps)
		svdCmat <- La.svd(iCmat)
		D <- svdCmat$u %*% diag(sqrt( svdCmat$d )) %*% t(svdCmat$v)
		return( t(D) )
		}
			
	lamTrans <- function(Vmat, lambda) {
		V1 <- Vmat
		diag(V1) <- 0
		V2 <- diag( diag(Vmat), ncol = length(Vmat[1,]), nrow = length(Vmat[,1]))
		Vmat <- V1 * lambda + V2
		return(Vmat)
		}
			
	resVar <- function(y, Vmat, p) {
		iV <- solve(Vmat, tol = .Machine$double.eps)
		e <- y - p
		s2 <- crossprod(e, iV %*% e)
		if( s2 < 0) { cat("ERROR -- negative variance\n")
			s2 <- s2 * -1}
		n <- length(y) 
		return( s2 / (n- nx) )
		}
		
	# Estimates the GLS parameters for given data
	get.coeffs <- function(Y, V, X) {
		iV <- solve(V, tol = .Machine$double.eps)
		xVix <- crossprod(X, iV %*% X)
		xViy <- crossprod(X, iV %*% Y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy 	#This is  a bad thing to do!!!!
		return(mu)
		}

	# Estimates the variance of a given trait (accounting for phylogeny)
	est.var <- function(y, V, x, mu ) {
		iV <- solve(V, tol = .Machine$double.eps)
		e <- y - x %*% mu
		s2 <- crossprod(e, iV %*% e)
		n <- length(y) 
		k <- length(x[1,])
		return( s2 / (n- k) )
		}
	
	# Full ML estimation for given x and V
	log.likelihood <- function(y, x, V) {
		mu <- get.coeffs(y, V, x)
		s2 <- est.var(y, V, x, mu)
		n <- length(x[,1])
		logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - 1)/2.0
		ypred <- x%*%mu	
		return( list(ll = ll, mu = mu, s2 = s2) )
		}
		
		
	null.var <- function(y, V) {
		X <- matrix(1, nrow = length(y))
		mu <- get.coeffs(y, V, X)
		return(est.var(y, V, X, mu))
		}
	

	Vmat <- as.matrix(phylomat)
	
	Vmat <- lamTrans(phylomat, lambda)
	
	prune.dat <- prune(data, Vmat)
	Vmat <- prune.dat$Vmat
	data <- prune.dat$dat
	nm <- names(data)

	
	
	n <- length(data[,1])
	
	
	
	# Get the design matrix
	m <- model.frame(formula, data)
	y <- m[,1]
	x <- model.matrix(formula, m)
	k <- length(x[1,])
	
	namey <- names(m)[1]
	
	ll <- log.likelihood(y, x, Vmat)
	
	log.lik <- ll$ll	

	aic <- -2 * log.lik + 2 * k
	aicc <- -2 * log.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
	
	coeffs <- ll$mu
	coeffs <- data.frame(t(coeffs))
	names(coeffs) <- colnames(x)
	varNames = names(m)

	
	pred <- x %*% ll$mu 
	
	res <- y - pred
	D <- Dfun(Vmat)
	pres <- D %*% res
	
	fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
	
	logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
 	
	logLikY <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(ll$s2) - logDetV / 2.0  - (n - 1 )/ 2.0
	
	RMS <- ll$s2
	RSSQ <- ll$s2 * (n - k)
	NMS <- RMS
	NSSQ <- RSSQ
	
	if(k > 0) {
		NMS <- null.var(y, Vmat)
		NSSQ <- NMS * (n - 1)
		}

	# Bits for parameter errors	
	errMat <- t(x)%*% solve(Vmat) %*% x  
	errMat <- solve(errMat) * RMS[1] 
	sterr <- diag(errMat)
	sterr <- sqrt(sterr)
	
	
	ret <- list(model = fm, formula = formula, logLikY = logLikY, RMS = RMS, NMS = NMS, NSSQ = NSSQ[1], RSSQ = RSSQ[1], 
	aic = aic, aicc = aicc, n = n, k = k, sterr = sterr, vcv = errMat, fitted = pred, residuals = res, phyres = pres, x = x, data = data,  varNames = varNames, y = y, V = Vmat, lambda = lambda, L0 = NULL, L1 = NULL, LamOptimised = FALSE, namey = namey)
	class(ret) <- "pglm"
	return(ret)
	
	}

# This returns the coefficients from the model
coef.pglm <- function(obj) { ret <- obj$model$coef 
							return(ret) }

# This returns the residuals from the model
residuals.pglm <- function(obj, phylo = FALSE) { ret <- NULL
				if(phylo == FALSE){ret <- obj$res} else {ret <- obj$phyres}
								return(ret) }

								
# This returns the fitted values							
fitted.pglm <- function(obj){ ret <- obj$fitted
								return(ret) }
# This predicts for given x
predict.pglm <- function(obj, x) { mu <- as.matrix(coef(obj) )
									ret <- cbind(1, x) %*% t(mu)
									return(ret) }
									
# This returns the AIC
AIC.pglm <- function(obj) { ret <- obj$aic
							return(ret[1]) }
							
# This returns the AICc
AICc.pglm <- function(obj) { ret <- obj$aicc
							return(ret[1]) }
							
# This returns the value of lambda at which the pglm was evaluated
lambda.pglm <- function(obj) { ret <- obj$lambda
								return(ret[1])}
							
# Very rough function for displaying a pglm output							
summary.pglm <- function(obj) {
		
		testLambda <- function(pobj) {
			
			lrt0 <- 2 * (pobj$logLikY - pobj$L0)
			lrt1 <- 2 * (pobj$logLikY - pobj$L1)
			
			p0 <- 1 - pchisq(lrt0, 1)
			p1 <- 1 - pchisq(lrt1, 1)
			
			cat("     Test of Lambda = 0: chisq = ", lrt0, " P = ", p0, "\n")
			cat("     Test of Lambda = 1: chisq = ", lrt1, " P = ", p1, "\n")
			}
			
	
	cat("\n\n--------------------------------------------------------\n")
	cat("Summary of Generalised Least Squares, correcting for \n")
	cat("Phylogeny:\n\n")
	cat("Number of parameters = ", obj$k,"\n")
	cat("Number of data points = ", obj$n,"\n\n")
	cat("Lambda statistic = ", obj$lambda, "\n")
	if(obj$LamOptimised == TRUE) { testLambda(obj)}
	cat("Maximised log-likelihood = ", obj$logLikY,"\n\n")
	cat("Model AIC = ", obj$aic, "\n")
	cat("Model AICc = ", obj$aicc, "\n\n")

	cat("Null Mean Square = ", obj$NSSQ,"\n")
	cat("Residual Mean Square = ", obj$RSSQ, "\n\n")
	cat("Raw R^2 = ", (obj$NSSQ - obj$RSSQ) / obj$NSSQ, "\n")
	cat("Adjusted R^2 = ", (obj$NMS - obj$RMS) / obj$NMS, "\n")
	
	Fstat <- ((obj$NSSQ - obj$RSSQ) / obj$RMS) / (obj$k - 1)
	
	cat("F statistic = ",  ((obj$NSSQ - obj$RSSQ) / obj$RMS) / (obj$k - 1), " ")
	cat("P model = ", pf(Fstat, obj$k - 1, obj$n - obj$k,  ncp=0, lower.tail = FALSE, log.p = FALSE), "\n\n")
	cat("Summary of coefficients:\n\n")
	coeffs <- coef(obj)
	errs <- obj$sterr
	cat("Term\t\tEstimate\t\tStd Err\t\tT-value\tP\n")
	storet<-c()
	for(i in 1:length(coeffs) ) {
		est <- coeffs[1,i]
		nm <- names(coeffs)[i]
		se <- errs[i]
		Tstat <- est / se
		storet<-c(storet,Tstat)	
		Pval <- 2 * ( 1 - pt( abs(Tstat), obj$n - obj$k) )
		cat(nm,"\t")
		cat(est, "\t", se, "\t", Tstat, "\t", Pval, "\n")
		}
	
	cat("\n\n--------------------------------------------------------\n")
	}

# Simple print function
print.pglm <- function(obj) {
	return(summary(obj))
	}


# Diagnostic plots
plot.pglm <- function(obj) {
	layout(matrix(c(1,2,3,4), 2, 2, byrow = FALSE))
	res <- residuals(obj, phylo = TRUE)
	res <- res / sqrt(var(res))[1]
	truehist(res, xlab = "Residual value (corrected for phylogeny)")
	qqnorm(res)
	abline(0, 1)
	plot(fitted(obj), res, xlab = "Fitted value", ylab = "Residual value (corrected for phylogeny)"  )
	plot(obj$y, fitted(obj), xlab = "Observed value", ylab = "Fitted value")
	}

# ----------------------------------------------------------------------
# This optimises the value of lambda for a given dataset
# and phylogeny. The return is the model at the ML value of
# lambda together with the ML estimate of lambda and the maximised
# log-likelihood. Also returned are the logLikelihoods at Lambda = 0
# and Lambda = 1.0
pglmEstLambda <- function(formula, data, phylomat, plotit = FALSE,  ...) {
	
	ll.fun <- function(lam) {
		pg <- pglm(formula, data, phylomat, lam)
		ll <- pg$logLikY
		return( ll )
		}
	
	oL <- optimize( ll.fun, interval = c(0,1), maximum = TRUE )
	L1 <- ll.fun(1)
	L0 <- ll.fun(0)
	fm <- pglm(formula, data, phylomat, oL$maximum)
	
	if(plotit == TRUE) {
		lambda <- seq(0, 1, by = 0.01)
		logLikelihood <- sapply(lambda, ll.fun)
		plot(lambda, logLikelihood, type = "l")
		maxlikelambda<-which(logLikelihood==max(logLikelihood))[1]
		sigdiff<-logLikelihood[maxlikelambda]-(.5*3.8414587)
		if(lambda[maxlikelambda]<0.01){lowerlamb<-0}
		if(lambda[maxlikelambda]>0.99){upperlamb<-1}
		if(lambda[maxlikelambda]>=0.01){
		lowerlamb<-lambda[which(abs(logLikelihood[1:maxlikelambda]-sigdiff)==min(abs(logLikelihood[1:maxlikelambda]-sigdiff)))]}
		if(lambda[maxlikelambda]<=0.99){
		upperlamb<-lambda[maxlikelambda+which(abs(logLikelihood[maxlikelambda:length(lambda)]-sigdiff)==min(abs(logLikelihood[maxlikelambda:length(lambda)]-sigdiff)))]}
		}
		
	
	fm$L1 <- L1
	fm$L0 <- L0
	fm$LamOptimised <- TRUE
	ret <- fm
	return( ret )
	}

# ----------------------------------------------------------------------
# Performs an ANOVA on the pglm Object using Sequential Sums of Squares
#

anova.pglm <- function(pglmObj) {
				V <- pglmObj$V
				dat <- pglmObj$data
				lambda <- pglmObj$lambda
				x <- data.frame( pglmObj$x )
				nm.x <- attr( terms(pglmObj$formula), "term.labels")
				k <- length(nm.x)
				nm.x <- nm.x[1:k]
				SSQ <- matrix(0, nrow = k + 1)
				remDF<- matrix(0, nrow = k + 1)
				SSQ[1] <- pglmObj$NSSQ
				remDF[1] <- pglmObj$n - 1
				for( i in 1: k ) {
						xv <- paste(nm.x[1:i], sep = "")
						fmla <- as.formula(paste( pglmObj$namey, " ~ ", paste(xv, collapse = "+") ))
						plm <- pglm(fmla, dat, V, lambda)
						SSQ[i + 1] <- plm$RSSQ
						remDF[i + 1] <- remDF[1] - plm$k + 1
				}
				errorDF <- remDF[k]
				errorMS <- SSQ[k] / errorDF
				termDF <- remDF[1:k] - remDF[2:(k+1)]
				termSSQ <- SSQ[1:k] - SSQ[2:(k+1)]
				MS <- termSSQ / termDF
				F <- MS / errorMS
				pF <-  df(F, termDF, errorDF)
				cat("Source\tDF\t\tSSQ\t\tMS\t\tF\n")
						for(i in 1:k) {
						cat(nm.x[i], "\t\t", termDF[i], "\t", termSSQ[i], "\t",
						MS[i], "\t", F[i])
						cat("   { P = ", pF[i], " }\n")
				}
				cat("Error\t", errorDF, "\t", SSQ[k], "\t", errorMS, "\n")
}


# ------------------------------------------------------------------------


