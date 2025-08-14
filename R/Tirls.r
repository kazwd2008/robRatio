#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#' @name  Tirls.aad
#' @title Tirls.aad: Robust estimator for the linear regression model
#'               with Tukey's biweight function and AAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm 
#'
#' @param x1 explanatory variable(s)
#' @param y1 objective variable 
#' @param rt sample weights
#' @param c1 tuning parameter from 4 to 8 for the scale parameter of AAD(Average Absolute Deviation) 
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#' \describe{
#'   \item{\code{TK}}{results of robust regression}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#' }
#' @export
#-------------------------------------------------------------------------------
Tirls.aad <- function(x1, y1, rt=rep(1, length(y1)), c1=8, rp.max=150, cg.rt=0.01) {

  s1.cg <- rep(0, rp.max)               # save changes of scale (2012.12.18)
  R0 <- stats::lm(y1~x1, weights=rt)           # initial estimation by OLS
  
  Tk2 <- R0
  rp1 <- 1				# iteration counter
  s0 <- 0			        # initial value for the scale of residuals
  s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))       # Average Absolute Deviation(AAD)

  #### calculate robust weights
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### iteration
  for (i in 2:rp.max) {                                 # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                   # 2012.12.19
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- stats::lm(y1~x1, weights=w1*rt) 
        rp1 <- rp1 + 1			                # counting iteration
        s1 <- s1.cg[rp1] <- mean(abs(Tk2$residuals))	# AAD
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#' @name  Tirls.mad
#' @title Tirls.mad: Robust estimator for the linear regression model
#'               with Tukey's biweight function and MAD scale
#'               by iteratively re-weighted least squares (IRLS) algorithm 
#'
#' @param x1 explanatory variable(s)
#' @param y1 objective variable 
#' @param c1 tuning parameter from 5.01 to 10.03 for the scale parameter of MAD(Median Absolute Deviation)
#' @param rt sample weights
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#' \describe{
#'   \item{\code{TK}}{results of robust regression}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#' }
#' @export
#-------------------------------------------------------------------------------
Tirls.mad <- function(x1, y1, rt=rep(1, length(y1)), c1=10.03, rp.max=150, cg.rt=0.01) {
  s1.cg <- rep(0, rp.max)               # save changes of scale (2011.03.07)
  R0 <- stats::lm(y1~x1, weights=rt)    # initial estimation by OLS
  
  Tk2 <- R0
  rp1 <- 1				# iteration counter
  s0 <- 0			        # initial value for the scale of residuals
  s1 <- s1.cg[rp1] <- stats::mad(Tk2$residuals) 	    # standardised MAD

  #### calculate robust weights
  	u1 <-Tk2$residuals/(c1*s1)
  	w1 <- (1-u1**2)**2
  	w1[which(abs(u1)>=1)] <- 0

  #### iteration
  for (i in 2:rp.max) {                             # 2011.03.07
      # while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.18
        Tk1 <-Tk2
        s0 <- s1	
        Tk2 <- stats::lm(y1~x1, weights=w1*rt) 
        rp1 <- rp1 + 1			        # counting iteration
        s1 <- s1.cg[rp1] <- stats::mad(Tk2$residuals)	# MAD
        u1 <-Tk2$residuals/(c1*s1)
        w1 <- (1-u1**2)**2
        w1[which(abs(u1)>=1)] <- 0
    }
	return(list(TK=Tk2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-------------------------------------------------------------------------------
