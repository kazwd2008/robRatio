#################################################################################################
#################################################################################################
#  Robust estimators for the generalised ratio model
#         by iteratively re-weighted least squares (IRLS) algorithm
#    Weight function:  Tukey's biweight function
#    Scale: Average absolute deviation (AAD) or median absolute deviation (MAD)
#------------------------------------------------------------------------------------------------#
#      Impremented by K. Wada (NSTAC, Japan)
#------------------------------------------------------------------------------------------------#
#       Created from Tirls.r  Ver. 4.0 [2014/09/05] for regression model
#       Ver.0     [2017/03/27]  Functions with AAD scale are implemented.
#       Ver.1     [2017/04/10]  - Functions with MAD scale are added.
#               		- An error flag added, which become 1 and abort the function
#	                          when all weights (w1) become zero.
#------------------------------------------------------------------------------------------------#
#  Functions
# 	RrTa.aad: Gamma = 1,   AAD scale
# 	RrTb.aad: Gamma = 1/2, AAD scale (conventional ratio model)
# 	RrTc.aad: Gamma = 0,   AAD scale (single regression without intercept)
# 	RrTa.mad: Gamma = 1,   MAD scale*
# 	RrTb.mad: Gamma = 1/2, MAD scale* (conventional ratio model)
# 	RrTc.mad: Gamma = 0,   MAD scale* (single regression without intercept)
#           * Since the mad function in R returns values corresponding to the standard deviation
#              (SD), tuning constants for SD are used instead of those for MAD.
#------------------------------------------------------------------------------------------------#
#  Parameters
#   x1      single explanatory variable
#   y1      objective variable
#   c1      tuning parameter for Tukey's biweight function
#   rp.max  maximum number of iteration (default setting : 100)
#------------------------------------------------------------------------------------------------#
#   Recommended range of tuning parameter for Tukey's biweight function
#
#              More robust <----> Less robust | default setting
#      ---------+--------+--------+-----------+-------------------
#	 SD     |  5.01  |  7.52  |  10.03    |  10.03
#	 AAD    |  4     |  6     |   8       |   8
#  MAD    |  3.38  |  5.07  |   6.76    | * Use the value for SD for mad function in R
#------------------------------------------------------------------------------------------------#
#  Returned values
#   par     robustly estimated rate y1/x1
#   wt      robust weights
#   rp      total number of iteration
#   s1      changes of the scale (AAD or MAD)
#   efg	    error flag, 1: acalculia (all weights become zero)  0: successful termination
#################################################################################################
#' @name  RrTa.aad
#' @title RrTa.aad: Robust estimator for a generalized ratio model (gamma=1)
#'        with Tukey biweight function and AAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(average absolute deviation)
#'           and from 5.01 to 10.03 for the scale parameter of MAD(median absolute deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   par   robustly estimated ratio of y1 to x1
#'   res   homoscedastic quasi-residuals
#'   wt    robust weights
#'   rp    total number of iteration
#'   s1    changes of the scale (AAD or MAD)
#'   efg   error flag. 1: acalculia (all weights become zero)  0: successful termination
#'
#' @export
#################################################################################################
# 	RrTa :   gamma = 1
#------------------------------------------------------------------------------------------------#
RrTa.aad <- function(x1, y1, c1=8, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- mean(y1 / x1)	                     	 # initial estimation
  res <- y1 / x1 - par			         # homoscedastic quasi-residuals
  rp1 <- 1				         # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))       # AAD scale

  #### calculating weights
  	u1 <-res/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0
	if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1 * y1 / x1) / sum(w1)     	# robust estimation with weights
      res <- y1 / x1 - par			# homoscedastic quasi-residuals
      rp1 <- rp1 + 1				# number of iteration
      s1 <- s1.cg[rp1] <- mean(abs(res))	# AAD scale
      u1 <-res/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0
      if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break           # convergence condition
      s0 <- s1
    }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
#' @name  RrTa.mad
#' @title RrTa.mad: Robust estimator for a generalized ratio model (gamma=1)
#'        with Tukey biweight function and MAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(average absolute deviation)
#'           and from 5.01 to 10.03 for the scale parameter of MAD(median absolute deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   par   robustly estimated ratio of y1 to x1
#'   res   homoscedastic quasi-residuals
#'   wt    robust weights
#'   rp    total number of iteration
#'   s1    changes of the scale (AAD or MAD)
#'   efg   error flag. 1: acalculia (all weights become zero)  0: successful termination
#'
#' @importFrom stats mad
#' @export
#------------------------------------------------------------------------------------------------#
RrTa.mad <- function(x1, y1, c1=10.03, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)               	 # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- mean(y1 / x1)	                	 # initial estimation
  res <- y1 / x1 - par				 # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)    	 # MAD scale

  #### calculating weights
  	u1 <-res/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0
	if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1 * y1 / x1) / sum(w1)     	# robust estimation with weights
      res <- y1 / x1 - par			# homoscedastic quasi-residuals
      rp1 <- rp1 + 1				# number of iteration
      s1 <- s1.cg[rp1] <- mad(res)	        # MAD scale
      u1 <-res/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0
      if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break           # convergence condition
      s0 <- s1
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#################################################################################################
# 	RetTb :   gamma = 1/2 (conventional ratio estimator)
#------------------------------------------------------------------------------------------------#
#' @name  RrTb.aad
#' @title Robust estimator for conventional ratio model (gamma=1/2)
#'        with Tukey biweight function and AAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(average absolute deviation)
#'           and from 5.01 to 10.03 for the scale parameter of MAD(median absolute deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   par   robustly estimated ratio of y1 to x1
#'   res   homoscedastic quasi-residuals
#'   wt    robust weights
#'   rp    total number of iteration
#'   s1    changes of the scale (AAD or MAD)
#'   efg   error flag. 1: acalculia (all weights become zero)  0: successful termination
#'
#' @export
#------------------------------------------------------------------------------------------------#
RrTb.aad <- function(x1, y1, c1=8, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1) # prevent overflow

  s1.cg <- rep(0, rp.max)                	# preserve changes in s1 (scale)
  efg <- 0					# error flag
  par <- sum(y1) / sum(x1)     	      		# initial estimation
  res <- y1/sqrt(x1) - par*sqrt(x1)	 	# homoscedastic quasi-residuals
  rp1 <- 1					# number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))      # AAD scale

  #### calculating weights
   u1 <-res/(c1*s1)
   w1 <- (1-u1**2)**2
   w1[which(abs(u1)>=1)] <- 0
   if (sum(w1)==0)  return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * w1) / sum(x1 * w1 )     # robust estimation with weights
        res   <-  y1 / sqrt(x1) - par * sqrt(x1)  # homoscedastic quasi-residuals
        rp1 <- rp1 + 1				  # number of iteration
        s1 <- s1.cg[rp1] <- mean(abs(res))        # AAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break           # convergence condition
        s0 <- s1
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
#' @name  RrTb.mad
#' @title RrTb.mad: Robust estimator for conventional ratio model (gamma=1/2)
#'        with Tukey biweight function and MAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(average absolute deviation)
#'           and from 5.01 to 10.03 for the scale parameter of MAD(median absolute deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   par   robustly estimated ratio of y1 to x1
#'   res   homoscedastic quasi-residuals
#'   wt    robust weights
#'   rp    total number of iteration
#'   s1    changes of the scale (AAD or MAD)
#'   efg   error flag. 1: acalculia (all weights become zero)  0: successful termination
#'
#' @importFrom stats mad
#' @export
#------------------------------------------------------------------------------------------------#
RrTb.mad <- function(x1, y1, c1=10.03, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                # preserve changes in s1 (scale)
  efg <- 0				 # error flag
  par <- sum(y1) / sum(x1)     	         # initial estimation
  res <- y1/sqrt(x1) - par*sqrt(x1)	 # homoscedastic quasi-residuals
  rp1 <- 1				 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)     # MAD scale

  #### calculating weights
   u1 <-res/(c1*s1)
   w1 <- (1-u1**2)**2
   w1[which(abs(u1)>=1)] <- 0
   if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * w1) / sum(x1 * w1 )     	# robust estimation with weights
        res   <-  y1 / sqrt(x1) - par * sqrt(x1)	# homoscedastic quasi-residuals
        rp1 <- rp1 + 1					# number of iteration
        s1 <- s1.cg[rp1] <- mad(res)                    # MAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break                 # convergence condition
        s0 <- s1
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#################################################################################################
# 	RrTc :   gamma = 1 (a single regression without intercept)
#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
#' @name  RrTc.aad
#' @title Robust estimator for a single regression model without intercept (gamma=0)
#'        with Tukey biweight function and AAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(average absolute deviation)
#'           and from 5.01 to 10.03 for the scale parameter of MAD(median absolute deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   par   robustly estimated ratio of y1 to x1
#'   res   homoscedastic quasi-residuals
#'   wt    robust weights
#'   rp    total number of iteration
#'   s1    changes of the scale (AAD or MAD)
#'   efg   error flag. 1: acalculia (all weights become zero)  0: successful termination
#'
#' @export
#------------------------------------------------------------------------------------------------#
RrTc.aad <- function(x1, y1, c1=8, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1 * x1) / sum(x1**2)               # initial estimation
  res <- y1 - par*x1	                         # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))       # AAD scale

  #### calculating weights
   u1 <-res/(c1*s1)
   w1 <- (1-u1**2)**2
   w1[which(abs(u1)>=1)] <- 0
   if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * x1 * w1) / sum(x1 **2 * w1 )    # robust estimation with weights
        res   <-  y1  - par * x1		          # homoscedastic quasi-residuals
        rp1 <- rp1 + 1					  # number of iteration
        s1 <- s1.cg[rp1] <- mean(abs(res))          	  # AAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break                   # convergence condition
        s0 <- s1
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
#' @name  RrTc.mad
#' @title Robust estimator for a single regression model without intercept (gamma=0)
#'        with Tukey biweight function and MAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(average absolute deviation)
#'           and from 5.01 to 10.03 for the scale parameter of MAD(median absolute deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   par   robustly estimated ratio of y1 to x1
#'   res   homoscedastic quasi-residuals
#'   wt    robust weights
#'   rp    total number of iteration
#'   s1    changes of the scale (AAD or MAD)
#'   efg   error flag. 1: acalculia (all weights become zero)  0: successful termination
#'
#' @importFrom stats mad
#' @export
#------------------------------------------------------------------------------------------------#
RrTc.mad <- function(x1, y1, c1=10.03, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                   	 # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1 * x1) / sum(x1**2)     	    	 # initial estimation
  res <- y1 - par*x1	                    	 # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)        	 # MAD scale

  #### calculating weights
  u1 <-res/(c1*s1)
  w1 <- (1-u1**2)**2
  w1[which(abs(u1)>=1)] <- 0
  if (sum(w1)==0)  return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * x1 * w1) / sum(x1 **2 * w1 )  # robust estimation with weights
        res   <-  y1  - par * x1		        # homoscedastic quasi-residuals
        rp1 <- rp1 + 1					# number of iteration
        s1 <- s1.cg[rp1] <- mad(res)                    # MAD scale
        u1 <-res/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
	    if (sum(w1)==0) return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break                 # convergence condition
        s0 <- s1
    }
	return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }
#################################################################################################
#################################################################################################
