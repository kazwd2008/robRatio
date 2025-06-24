#################################################################################################
#################################################################################################
#  Robust estimators for the generalised ratio model
#         by iteratively re-weighted least squares (IRLS) algorithm
#    Weight function:  Huber's weight function
#    Scale: Average absolute deviation (AAD) or median absolute deviation (MAD)
#------------------------------------------------------------------------------------------------#
#      Impremented by K. Wada (NSTAC, Japan)
#------------------------------------------------------------------------------------------------#
#       Created from RtT.r Ver. 1  [2017/04/10]
#       Ver.0     [2019/06/19]  Functions with AAD and MAD scale are implemented.
#------------------------------------------------------------------------------------------------#
#    Functions
# 	RrHa.aad:     Gamma = 1,   AAD scale
# 	RrHb.aad:     Gamma = 1/2, AAD scale (conventional ratio model)
# 	RrHc.aad:     Gamma = 0,   AAD scale (single regression without intercept)
# 	RrHa.mad:     Gamma = 1,   MAD scale*
# 	RrHb.mad:     Gamma = 1/2, MAD scale* (conventional ratio model)
# 	RrHc.mad:     Gamma = 0,   AAD scale* (single regression without intercept)
#           * Since the mad function in R returns values corresponding to the standard deviation
#              (SD), tuning constants for SD are used instead of those for MAD.
#------------------------------------------------------------------------------------------------#
#  Parameters
#   x1      single explanatory variable
#   y1      objective variable
#   c1      tuning parameter for Huber's weight function (Huber's k)
#   rp.max  maximum number of iteration (default setting : 100)
#------------------------------------------------------------------------------------------------#
#   Recommended range of tuning parameter c1 for Huber's weight function
#
#              More robust <----> Less robust | default setting
#      ---------+--------+--------+-----------+------------------
#	       SD     |  1.44  |  2.16  |   2.88    |   2.88
#	       AAD    |  1.15  |  1.72  |   2.30    |   2.30
#        MAD    |  0.97  |  1.46  |   1.94    | * Use the value for SD for mad function in R
#------------------------------------------------------------------------------------------------#
#  Returned values
#   par     robustly estimated rate y1/x1
#   wt      robust weights
#   rp      total number of iteration
#   s1      changes of the scale (AAD or MAD)
#   efg	    error flag, 1: acalculia (all weights become zero)  0: successful termination
#################################################################################################
#' @name  RrHa.aad
#' @title RrHa.aad: Robust estimator for a generalized ratio model (gamma=1)
#'               with Huber's weight function and AAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable 
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(average absolute deviation)
#'           and from 1.44 to 2.88 for the scale parameter of MAD(median absolute deviation)
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
# 	RrHa :   gamma = 1
#------------------------------------------------------------------------------------------------#
RrHa.aad <- function(x1, y1, c1=2.3, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1) # prevent overflow

  s1.cg <- rep(0, rp.max)               	# preserve changes in s1 (scale)
  efg <- 0					# error flag
  par <- mean(y1 / x1)	                        # initial estimation
  res <- y1 / x1 - par			        # homoscedastic quasi-residuals
  rp1 <- 1					# number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))      # AAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
   if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1 * y1 / x1) / sum(w1)	# robust estimation with weights
      res <- y1 / x1 - par			# homoscedastic quasi-residuals
      rp1 <- rp1 + 1				# number of iteration
      s1 <- s1.cg[rp1] <- mean(abs(res))	# AAD scale

      w1 <- s1*c1 / abs(res)
      w1[which(abs(res) <= s1*c1)] <- 1
      if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break           # convergence condition
      s0 <- s1
   }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
#' @name  RrHa.mad
#' @title Robust estimator for a generalized ratio model (gamma=1)
#'               with Huber's weight function and MAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable 
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(average absolute deviation)
#'           and from 1.44 to 2.88 for the scale parameter of MAD(median absolute deviation)
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
RrHa.mad <- function(x1, y1, c1=2.88, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1) # prevent overflow

  s1.cg <- rep(0, rp.max)               	# preserve changes in s1 (scale)
  efg <- 0					# error flag
  par <- mean(y1 / x1)	                        # initial estimation
  res <- y1 / x1 - par			        # homoscedastic quasi-residuals
  rp1 <- 1					# number of iteration
  s0 <- s1 <- s1.cg[rp1] <-  mad(res)           # MAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
   if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1 * y1 / x1) / sum(w1)	# robust estimation with weights
      res <- y1 / x1 - par			# homoscedastic quasi-residuals
      rp1 <- rp1 + 1				# number of iteration
      s1 <- s1.cg[rp1] <-  mad(res)             # MAD scale

      w1 <- s1*c1 / abs(res)
      w1[which(abs(res) <= s1*c1)] <- 1
      if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break           # convergence condition
      s0 <- s1
   }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#################################################################################################
# 	RrHb :   gamma = 1/2 (conventional ratio model)
#------------------------------------------------------------------------------------------------#
#' @name  RrHb.aad
#' @title Robust estimator for a cconventional ratio model (gamma=1/2)
#'               with Huber's weight function and AAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(average absolute deviation)
#'           and from 1.44 to 2.88 for the scale parameter of MAD(median absolute deviation)
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
RrHb.aad <- function(x1, y1, c1=2.3, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1) / sum(x1)     	                 # initial estimation
  res <- y1/sqrt(x1) - par*sqrt(x1)	         # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))       # AAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
   if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par   <- sum(y1 * w1) / sum(x1 * w1 )      # robust estimation with weights
      res   <-  y1 / sqrt(x1) - par * sqrt(x1)	 # homoscedastic quasi-residuals
      rp1 <- rp1 + 1				 # number of iteration
      s1 <- s1.cg[rp1] <- mean(abs(res))         # AAD scale

      w1 <- s1*c1 / abs(res)
      w1[which(abs(res) <= s1*c1)] <- 1
      if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break            # convergence condition
      s0 <- s1
    }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
#' @name  RrHb.mad
#' @title Robust estimator for a cconventional ratio model (gamma=1/2)
#'               with Huber's weight function and MAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(average absolute deviation)
#'           and from 1.44 to 2.88 for the scale parameter of MAD(median absolute deviation)
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
RrHb.mad <- function(x1, y1, c1=2.88, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1) / sum(x1)     	                 # initial estimation
  res <- y1/sqrt(x1) - par*sqrt(x1)	         # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <-  mad(res)            # MAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
   if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par   <- sum(y1 * w1) / sum(x1 * w1 )      # robust estimation with weights
      res   <-  y1 / sqrt(x1) - par * sqrt(x1)	 # homoscedastic quasi-residuals
      rp1 <- rp1 + 1				 # number of iteration
      s1 <- s1.cg[rp1] <-  mad(res)              # MAD scale

      w1 <- s1*c1 / abs(res)
      w1[which(abs(res) <= s1*c1)] <- 1
      if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
      if (abs(1-s1/s0) < cg.rt) break            # convergence condition
      s0 <- s1
    }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }

#################################################################################################
# 	RrHc :   gamma = 1 (a single regression without intercept)
#------------------------------------------------------------------------------------------------#
#' @name  RrHc.aad
#' @title Robust estimator for a single regression model without intercept (gamma=0)
#'               with Huber's weight function and AAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(average absolute deviation)
#'           and from 1.44 to 2.88 for the scale parameter of MAD(median absolute deviation)
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
RrHc.aad <- function(x1, y1, c1=2.3, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1 * x1) / sum(x1**2)          	 	 # initial estimation
  res <- y1 - par * x1	                         # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))       # AAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
   if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * x1 * w1) / sum(x1 **2 * w1 )     # robust estimation with weights
        res   <-  y1  - par * x1		    # homoscedastic quasi-residuals
        rp1 <- rp1 + 1				    # number of iteration
        s1 <- s1.cg[rp1] <- mean(abs(res))          # AAD scale
        w1 <- s1*c1 / abs(res)
        w1[which(abs(res) <= s1*c1)] <- 1
        if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break             # convergence condition
        s0 <- s1
    }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }
#------------------------------------------------------------------------------------------------#
#' @name  RrHc.mad
#' @title Robust estimator for a single regression model without intercept (gamma=0)
#'               with Huber's weight function and MAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(average absolute deviation)
#'           and from 1.44 to 2.88 for the scale parameter of MAD(median absolute deviation)
#' @param dat name of the dataframe containing x1 and y1 if any
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
RrHc.mad <- function(x1, y1, c1=2.88, dat="", rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1 * x1) / sum(x1**2)     	 	 # initial estimation
  res <- y1 - par*x1	                         # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)             # MAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
   if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))

  #### iteration
  for (i in 2:rp.max) {
        par.bak <- par
        res.bak <- res
        par   <- sum(y1 * x1 * w1) / sum(x1 **2 * w1 )     # robust estimation with weights
        res   <-  y1  - par * x1		    # homoscedastic quasi-residuals
        rp1 <- rp1 + 1				    # number of iteration
        s1 <- s1.cg[rp1] <-  mad(res)               # MAD scale
        w1 <- s1*c1 / abs(res)
        w1[which(abs(res) <= s1*c1)] <- 1
        if (sum(w1)==0)   return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=1))
        if (abs(1-s1/s0) < cg.rt) break             # convergence condition
        s0 <- s1
    }
    return(list(par = par, res=res, wt=w1, rp=rp1, s1=s1.cg, efg=efg))
  }
#################################################################################################
#################################################################################################
