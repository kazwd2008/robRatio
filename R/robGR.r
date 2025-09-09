####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####
#`
#' @title Robust estimator for a generalized ratio model
#'
#' @description This function simultaneously estimates two parameters of 
#' the generalized ratio model <doi:10.17713/ajs.v50i1.994>. 
#' It uses Tukey's biweight function and AAD for scale of quasi residuals.
#'
####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|
#' @name robGR
#' @param x1 single explanatory variable (a vector)
#' @param y1 objective variable to be imputed (a vector)
#' @param g1 initial gamma value (default g1=0.5)
#' @param c1 tuning constant for Tukey's biweight function. Supposed to choose 4 to 8. 
#'           Smaller figure is more robust (default tp=8).
#' @param rp.max maximum number of iteration (default: rp.max=50)
#' @param cg.rt convergence condition to stop iteration (default: cg.rt=0.001)
#' @return a list with the following elements
#' \describe{
#'   \item{\code{par}}{robustly estimated ratio of y1 to x1 (beta)}
#'   \item{\code{g1}}{robustly estimated power (gamma)}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{efg}}{error flag. 1: calculation not coverged,  0: successful termination}
#'   \item{\code{rt.cg}}{change of par(beta)}
#'   \item{\code{g1.cg}}{changes of g1(gamma)}
#'   \item{\code{s1.cg}}{changes of the scale(AAD)}
#' }

####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|
robGR <- function(x1, y1, g1=0, c1=8, rp.max=100, cg.rt=0.001){

  s1.cg <- rt.cg <- g1.cg  <- rep(NA, (rp.max*4))	# for preserving histories
  rp1     <- c(1, 0, 0, 0)  				# count of iteration for each step
  s0      <- 0	    					# initial value of scale parameter

####---------------------------------------------------------------------------------------(1)
   Trt1 <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))	# ratio estimation
   rs1.x  <- y1 - Trt1 * x1	                        # heteroscedastic residuals (2025.08.26)
   rs1   <- y1 / x1^g1 - Trt1 * x1^(1-g1)		# homoscedastic quasi-residuals

   s1 <- s1.cg[rp1[1]] <- mean(abs(rs1))	        # AAD scale of quasi-residuals
   g1.cg[rp1[1]]  <- g1
   rt.cg[rp1[1]]  <- Trt1
####---------------------------------------------------------------------------------------
#### Non-robust Estimation of rate (Trt1) and power (g1) (without weight)
####---------------------------------------------------------------------------------------
       rp1[1] <- rp1[1] + 1

       ix1 <- which((x1 !=0) & (rs1.x !=0))
       x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
       g1 <- (solve(t(x2) %*% x2) %*% t(x2) %*% log(abs(rs1.x[ix1])))[2]   # estimate power

       Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))		   # estimate ratio
       rs1.x  <- y1 - Trt1 * x1	                # heteroscedastic residuals 
       rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	# homoscedastic quasi-residuals 

       s0 <- s1					# preserve the previous value (AAD)
       s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))	# AAD of current quasi residuals

       if (s1 > s0) {        # abandon new g1 
          rp1[1] <- rp1[1] - 1
          s1 <- s0
          g1 <- g1.cg[sum(rp1)]  
          Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))
       } else {  
       g1.cg[sum(rp1)]  <- g1
       rt.cg[sum(rp1)]  <- Trt1
       }

####---------------------------------------------------------------------------------------(2)
####  Iterative robust estimation of the ratio "Trt1" 
####                               with a fixed power "g1" derived in step(1)
####---------------------------------------------------------------------------------------
  for (i in 1:rp.max) {
      rp1[2] <- rp1[2] + 1

    # renew weights 
      u1 <-rs1/(c1*s1)
      w1 <- (1-u1**2)**2
      w1[which(abs(u1)>=1)] <- 0

    # prevent all zero weights
      ix1 <- which((w1*x1 !=0) & (rs1.x !=0))    # remove observations that make Trt1 and rs1 NaN
       if (length(ix1) < (length(x1)/2)) {                     # reset w1 
          w1 <- rep(1, length(x1))
          ix1 <- which((w1*x1 !=0) & (rs1.x !=0))
      }

    # estimate Trt1 with fixed value of g1 obtained in the previous step
      Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
      rs1.x  <- y1 - Trt1 * x1		            # heteroscedastic residuals
      rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	    # homoscedastic quasi-residuals

      s0 <- s1					    # preserve previous AAD value
      s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))	    # preserve history of AAD scale
      g1.cg[sum(rp1)]  <- g1
      rt.cg[sum(rp1)]  <- Trt1
      if (s1 > s0) {        # abandon new g1 
          rp1[2] <- rp1[2] - 1
          s1 <- s0 
          g1 <- g1.cg[sum(rp1)]  
          Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))
       } else if (abs(1-s1/s0) < cg.rt) break 		    # convergence condition
  }

####---------------------------------------------------------------------------------------(3)
####  Iterative robust and simultaneous estimation of the ratio "Trt1" and the power "g1"
####---------------------------------------------------------------------------------------

  for (i in 1:rp.max) {
#  for (i in 1:4) {
      rp1[3] <- rp1[3] + 1

    # renew weights 
      u1 <-rs1/(c1*s1)			
      w1 <- (1-u1**2)**2		
      w1[which(abs(u1)>=1)] <- 0	

      ix1 <- which((w1*x1!=0) & (rs1.x !=0))  # remove observations that make Trt1 and rs1 NaN
      if (length(ix1) < (length(x1)/2)) {                     # reset w1 
          w1 <- rep(1, length(x1))
          ix1 <- which((w1*x1 !=0) & (rs1.x !=0))
      }
      x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))

      w2 <- matrix(0, ncol=length(ix1), nrow=length(ix1))
      diag(w2) <- w1[ix1]
      g1 <- (solve(t(x2) %*% w2 %*% x2) %*% t(x2) %*% w2 %*% log(abs(rs1.x[ix1])))[2]     # estimate power (g1)

      Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
      rs1.x  <- y1 - Trt1 * x1	                # heteroscedastic residuals (2025.08.26)
      rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	# homoscedastic quasi-residuals

      s0 <- s1
      s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))
      g1.cg[sum(rp1)]  <- g1
      rt.cg[sum(rp1)]  <- Trt1

      if (s1 > s0) {        # abandon new g1 
          rp1[3] <- rp1[3] - 1
          s1 <- s0 
          g1 <- g1.cg[sum(rp1)]  
          Trt1  <- sum(y1 * x1^(1-2*g1)) / sum(x1^(2*(1-g1)))
       } else if (abs(1-s1/s0) < cg.rt) break 		    # convergence condition
    }

####---------------------------------------------------------------------------------------(4)
####  If the estimation is not converge, reduce tuning constant slightly
####---------------------------------------------------------------------------------------

  if (abs(1-s1/s0) >= cg.rt) {  # when convergence condition is not met
     c1 <- c1 - 0.1		# reduce tuning constant slightly
     for (i in 1:rp.max) {
         rp1[4] <- rp1[4] + 1

       # renew weights 
         u1 <-rs1/(c1*s1)
         w1 <- (1-u1**2)**2
         w1[which(abs(u1)>=1)] <- 0

         ix1 <- which((x1*w1 !=0) & (rs1.x !=0))
         if (length(ix1) < (length(x1)/2)) {                     # reset w1 
            w1 <- rep(1, length(x1))		
            ix1 <- which((w1*x1 !=0) & (rs1.x !=0))
         }
         w2 <- matrix(0, ncol=length(ix1), nrow=length(ix1))
         diag(w2) <- w1[ix1]
         x2 <- cbind(rep(1, length(ix1)), log(x1[ix1]))
         g1 <- (solve(t(x2) %*% w2 %*% x2) %*% t(x2) %*% w2 %*% log(abs(rs1.x[ix1])))[2]     # estimate power (g1)
 
         Trt1  <- sum(w1[ix1] *y1[ix1] * (w1[ix1] * x1[ix1])^(1-2*g1)) / sum((w1[ix1] * x1[ix1])^(2*(1-g1))) 
         rs1.x  <- y1 - Trt1 * x1	                # heteroscedastic residuals 
         rs1    <- y1 / x1^g1 - Trt1 * x1^(1-g1)	# homoscedastic quasi-residuals

         s0 <- s1
         s1 <- s1.cg[sum(rp1)] <- mean(abs(rs1))
         g1.cg[sum(rp1)]  <- g1
         rt.cg[sum(rp1)]  <- Trt1

      if (abs(1-s1/s0) < cg.rt) break                   # convergence condition
    }
  }
    if (abs(1-s1/s0) < cg.rt) fg1 <- 0 else fg1 <- 1
    return(list(par=Trt1, g1=g1, res=rs1, wt=w1, rp=rp1, efg=fg1, rt.cg=rt.cg[1:length(sum(rp1))], g1.cg=g1.cg[1:length(sum(rp1))], s1.cg=s1.cg[1:length(sum(rp1))]))
}
