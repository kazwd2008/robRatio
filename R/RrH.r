#################################################################################################
#################################################################################################
#  Robust estimators for the generalised ratio model
#         by iteratively re-weighted least squares (IRLS) algorithm
#    Weight function:  Huber's weight function
#    Scale: Average absolute deviation (AAD) or median absolute deviation (MAD)
#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
#    Functions
# 	RrH.aad:     AAD scale (former RrHa.aad, RrHb.aad and RrHc.aad) 
# 	RrH.mad:     sandardized MAD scale (former RrHa.mad, RrHb.mad and RrHc.mad)
#
#################################################################################################
#' @name  RrH.aad
#' @title RrH.aad: Robust estimator for a generalized ratio model 
#'               with Huber's weight function and AAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable 
#' @param g1 power (default: g1=0.5(conventional ratio model))
#' @param c1 tuning parameter usually from 1.15 to 2.30 (smaller figure is more robust)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   \item{\code{par}}{robustly estimated ratio of y1 to x1}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#'   \item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
#'
#' @export
#################################################################################################
# 	RrHa :   gamma = 1
#------------------------------------------------------------------------------------------------#
RrH.aad <- function(x1, y1, g1=0.5, c1=2.3, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1) # prevent overflow

  s1.cg <- rep(0, rp.max)               	# preserve changes in s1 (scale)
  efg <- 0					# error flag
  par <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))   # initial estimation
  res <- y1 / x1^g1 - par * x1^(1-g1)	        # homoscedastic quasi-residuals
  rp1 <- 1					# number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mean(abs(res))      # AAD scale

  #### calculating weights
   w1 <- s1*c1 / abs(res)
   w1[which(abs(res) <= s1*c1)] <- 1
  ix1 <- which(((w1*x1)!=0) & (res!=0))  # remove observations that make par and res NaN
  if (length(ix1)==0) {           # reset w1 as all 1 when all the weights become zero 
     w1 <- rep(1, length(x1))
     ix1 <- which((w1*x1) !=0)
  }

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / 
             sum((w1[ix1] * x1[ix1])^(2*(1-g1)))  # robust estimation with weights
      res <- y1 / x1^g1 - par * x1^(1-g1) 	# homoscedastic quasi-residuals
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
#' @name  RrH.mad
#' @title Robust estimator for a generalized ratio model 
#'               with Huber's weight function and MAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable 
#' @param g1 power (default: g1=0.5(conventional ratio model))
#' @param c1 tuning parameter usually from 1.44 to 2.88 (equivalent to those for AAD scale)  
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   \item{\code{par}}{robustly estimated ratio of y1 to x1}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#'   \item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
#'
#' @importFrom stats mad
#' @export
#------------------------------------------------------------------------------------------------#
RrH.mad <- function(x1, y1, g1=0.5, c1=2.88, rp.max=100, cg.rt=0.01) {

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
  ix1 <- which(((w1*x1)!=0) & (res!=0))  # remove observations that make par and res NaN
  if (length(ix1)==0) {             # reset w1 as all 1 when all the weights become zero 
     w1 <- rep(1, length(x1))
     ix1 <- which((w1*x1) !=0)
  }

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / 
             sum((w1[ix1] * x1[ix1])^(2*(1-g1)))  # robust estimation with weights
      res <- y1 / x1^g1 - par * x1^(1-g1)	  # homoscedastic quasi-residuals
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
#################################################################################################
