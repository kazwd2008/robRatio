#-------------------------------------------------------------------------------
#' @name  Hirls.aad
#' @title Hirls.aad: Robust estimator for the linear regression model
#'               with Huber's weight function and AAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm 
#'
#' @param x1 explanatory variable(s)
#' @param y1 objective variable 
#' @param rt sample weights
#' @param c1 tuning parameter from 1.15 to 2.30 for the scale parameter of AAD(Average Absolute Deviation)
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#' \describe{
#'   \item{\code{HB}}{results of robust regression}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#' }
#' @export
#-------------------------------------------------------------------------------
Hirls.aad <- function(x1, y1, rt=rep(1, length(y1)), c1=1.15, rp.max=150, cg.rt=0.01) {

  s1.cg <- rep(0, rp.max)               # save changes of scale (2011.03.07)
  R0 <- stats::lm(y1~x1, weights=rt)    # initial estimation by OLS
  
  Hb2 <- R0
  rp1 <- 1                              # iteration counter
  s0 <- 0			        # initial value for the scale of residuals
  s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals))        # 平均絶対残差

  #### calculate robust weights
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### iteration
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- stats::lm(y1~x1, weights=w1*rt) 
      rp1 <- rp1 + 1			
      s1 <- s1.cg[rp1] <- mean(abs(Hb2$residuals))	# AAD
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#' @name  Hirls.mad
#' @title Hirls.mad: Robust estimator for the linear regression model
#'               with Huber's weight function and MAD scal
#'               by iteratively re-weighted least squares (IRLS) algorithm 
#'
#' @param x1 explanatory variable(s)
#' @param y1 objective variable 
#' @param c1 tuning parameter from 1.44 to 2.88 for the scale parameter of MAD(Median Absolute Deviation)
#' @param rt sample weights
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#'
#' @return a list with the following elements
#' \describe{
#'   \item{\code{HB}}{results of robust regression}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes in scale through iterative calculation}
#' }
#' @export
#-------------------------------------------------------------------------------
Hirls.mad <- function(x1, y1, rt=rep(1, length(y1)), c1=1.44, rp.max=150, cg.rt=0.01) {

  s1.cg <- rep(0, rp.max)               # save changes of scale (2011.03.07)
  R0 <- stats::lm(y1~x1, weights=rt)    # initial estimation by OLS
  
  Hb2 <- R0
  rp1 <- 1                              # iteration counter
  s0 <- 0			        # initial value for the scale of residuals
  s1 <- s1.cg[rp1] <- stats::mad(Hb2$residuals)        # standardised MAD

  #### calculate robust weights
    w1 <- s1*c1 / abs(Hb2$residuals)
    w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
  
  #### iteration
  for (i in 2:rp.max) {                             # 2011.03.07
    #  while (abs(1-s1/s0) >= 0.01) {
      if (abs(1-s1/s0) < cg.rt) break                # 2012.12.19
      Hb1 <-Hb2
      s0 <- s1	
      Hb2 <- stats::lm(y1~x1, weights=w1*rt) 
      rp1 <- rp1 + 1			
      s1 <- s1.cg[rp1] <- stats::mad(Hb2$residuals)	# MAD
      w1 <- s1*c1 / abs(Hb2$residuals)
      w1[which(abs(Hb2$residuals) <= s1*c1)] <- 1
    }
	return(list(HB=Hb2, wt=w1, rp=rp1, s1=s1.cg))
  }
#-------------------------------------------------------------------------------
