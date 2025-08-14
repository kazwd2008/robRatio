#################################################################################################
#################################################################################################
#  Robust estimators for the generalised ratio model
#         by iteratively re-weighted least squares (IRLS) algorithm
#    Weight function:  Tukey's biweight function
#    Scale: Average absolute deviation (AAD) or median absolute deviation (MAD)
#------------------------------------------------------------------------------------------------#
#  Functions
# 	RrT.aad: AAD scale (former RrTa.aad, RrTb.aad and RrTc.aad) 
# 	RrT.mad: standardized MAD scale (former RrTa.mad, RrTb.mad and RrTc.mad)
#
#################################################################################################
#' @name  RrT.aad
#' @title RrT.aad: Robust estimator for a generalized ratio model
#'        with Tukey biweight function and AAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param g1 power (default: g1=0.5(conventional ratio model))
#' @param c1 tuning parameter usually from 4 to 8 (smaller figure is more robust) 
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   \item{\code{par}}{robustly estimated ratio of `y1` to `x1`}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#'   \item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
#'
#' @export
#################################################################################################
RrT.aad <- function(x1, y1, g1=0.5, c1=8, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)                        # preserve changes in s1 (scale)
  efg   <- 0					 # error flag
  par   <- sum(y1*x1^(1-2*g1))/sum(x1^(2*(1-g1)))  # initial estimation
  res   <- y1/x1^g1 - par * x1^(1-g1)	         # homoscedastic quasi-residuals
  # res.x <- y1 * x1^g1 - par * x1^(1+g1)        # heteroscedastic residuals (2025.07.26)
  rp1   <- 1				         # number of iteration
  s0    <- s1 <- s1.cg[rp1] <- mean(abs(res))    # AAD scale

  #### calculating weights
  u1 <-res/(c1*s1)
  w1 <- (1-u1**2)**2
  w1[which(abs(u1)>=1)] <- 0

  ix1 <- which(((w1*x1)!=0) & (res!=0))  # remove observations that make par and res NaN
  if (length(ix1)==0) {             # reset w1 as all 1 when all the weights become zero 
     w1 <- rep(1, length(x1))
     ix1 <- which((w1*x1) !=0)
  }

  #### iteration
  for (i in 2:rp.max) {
      par.bak <- par
      res.bak <- res
      par <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
      rs1    <- y1 / x1^g1 - par * x1^(1-g1)	# homoscedastic quasi-residuals
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
#' @name  RrT.mad
#' @title RrT.mad: Robust estimator for a generalized ratio model 
#'        with Tukey biweight function and MAD scale
#'        by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable
#' @param g1 power (default: g1=0.5(conventional ratio model))
#' @param c1 tuning parameter usually from 5.01 to 10.03 (equivalent to those for AAD scale)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#'   \item{\code{par}}{robustly estimated ratio of `y1` to `x1`}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes of the scale (AAD or MAD)}
#'   \item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
#'
#' @importFrom stats mad
#' @export
#------------------------------------------------------------------------------------------------#
RrT.mad <- function(x1, y1, g1=0.5, c1=10.03, rp.max=100, cg.rt=0.01) {

  x1 <- as.numeric(x1);    y1 <- as.numeric(y1)  # prevent overflow

  s1.cg <- rep(0, rp.max)               	 # preserve changes in s1 (scale)
  efg <- 0					 # error flag
  par <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))  # initial ratio estimation
  res <- y1 / x1^g1 - par * x1^(1-g1)		 # homoscedastic quasi-residuals
  rp1 <- 1					 # number of iteration
  s0 <- s1 <- s1.cg[rp1] <- mad(res)    	 # standardized MAD scale

  #### calculating weights
  u1 <- res/(c1*s1)
  w1 <- (1-u1**2)**2
  w1[which(abs(u1)>=1)] <- 0
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
      rp1 <- rp1 + 1				  # number of iteration
      s1 <- s1.cg[rp1] <- mad(res)	          # MAD scale
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
#################################################################################################
