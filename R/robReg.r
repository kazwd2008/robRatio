####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|
#' @name robReg
#' @title \code{robReg} is for Robust regression by the IRLS algorithm.
#'    It integrates Tirls.r and Hirls.r. 
#-------------------------------------------------------------------------------
#' @param x1     explanatory variable in regression (a vector or a matrix)
#' @param y1     objective variable in regression (a vector)   
#' @param wf     weight function ("T" for Tukey's biweight, and "H" for Huber weight)
#' @param scale  scale for residuals. "AAD"(default) or "MAD".
#' @param rt     sample weight (default 1)
#' @param tp tuning parameter (tp=4, 6 or 8) for weight function. 
#'           Smaller figure is more robust.
#' @param rp.max The maximum number of iteration (default 150)  
#' @param cg.rt  convergence condition to stop iteration (default: cg.rt=0.001)
#-------------------------------------------------------------------------------
#' @return a list with the following elements
#' \describe{
#'   \item{\code{cond}}{Weight function, scale, and other arguments choosed}
#'   \item{\code{TK}}{robustly estimated regression coefficients using Tukey's biweight}
#'   \item{\code{HB}}{robustly estimated regression coefficients using Huber weight}
#'   \item{\code{wt}}{final robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{iterative changes in the sclae of residuals (AAD or MAD)}
#' }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#' @examples
#' \dontrun{
#' require(robRatio)
#'
#' set.seed(4)
#' cov1 <- matrix(c(3, 2.8, 2.8, 3), 2, 2)
#' cov2 <- matrix(c(2.5, 0, 0, 3), 2, 2)
#' dat1 <-  MASS::mvrnorm(n=400, mu=c(100, 100), Sigma=cov1, empirical=TRUE)
#' dat2 <-  cbind(runif(100,  min=96, max=104),  runif(50,  min=95, max=105))
#' dat3 <- matrix(c(103, 103.5, 104.5, 104.8, 96, 98, 94, 95), 4, 2) 
#' dat <- rbind(dat1, dat2, dat3)
#' plot(dat)
#' y1 <- dat[,2]
#' x1 <- dat[,1]
#' 
#' R0  <- lm(y1~x1)        # regression by OLS                                         
#' 
#' o1 <- robReg(x1, y1, tp=4)  # robust regression by IRLS (more robust)
#' o2 <- robReg(x1, y1, tp=8)  # robust regression by IRLS (less robust)
#' 
#' par(mfrow=c(2,2))
#' 
#' # non-robust regression
#'   plot(dat, pch=20, main="non-robust regression")
#'   abline(R0, col="red", lwd=2)
#' 
#' # robust regression with coloring robust weight
#'   f.o1 <- rep(1, length(x1))
#'   f.o1[which(o1$wt < 0.8)] <- 3
#'   f.o1[which(o1$wt < 0.5)] <- 7
#'   f.o1[which(o1$wt < 0.2)] <- 2
#'   f.o1[which(o1$wt == 0)] <- 8
#' 
#'   plot(x1, y1, pch=20, col=f.o1)
#'   abline(R0, col="red", lty=3)
#'   abline(o1$TK, col="blue", lwd=2)
#'   abline(o2$TK, col="cyan", lwd=2)
#' 
#' # robust weights (more robust)
#'   hist(o1$wt, main="tp=4(more robust)")
#' 
#' # robust weights (less robust)
#'   hist(o2$wt, main="tp=4(less robust)")
#' 
#' }
#'
#' @export
#'
####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|
robReg <- function(x1, y1, wf="T", scale="AAD", rt=rep(1, length(y1)), tp=8, 
                   rp.max=150, cg.rt=0.01) {

  vl.wf    <- c("T", "H")
  vl.scale <- c("AAD", "MAD")
  vl.tp    <- c(4, 6, 8)

  if (wf    %in% vl.wf    == FALSE) stop("Please choose 'T' or 'H' for wf")
  if (scale %in% vl.scale == FALSE) stop("Please choose 'AAD' or 'MAD' for scale")
  if (tp    %in% vl.tp    == FALSE) stop("Please choose 4, 6 or 8 for tp")
  if (length(rt) %in% c(1, length(y1)) ==FALSE) stop("Incorrect length of rt")

wfs <- paste(wf, scale, sep=".")
tp2 <- tp/2-1  # (4, 6, 8) => (1, 2, 3)

  if (wfs == "T.AAD") {          
        c1 <- switch(tp2, 4, 6, 8)
	rt <- Tirls.aad(x1, y1, rt, c1, rp.max, cg.rt)
        cd <- paste("wf=", wf, "scale=", scale, "tp=", tp, "c1=", c1, sep=",")
        return(list(cond=cd, TK=rt$TK, wt=rt$wt, rp=rt$rp, s1=rt$s1))
     } else if (wfs == "T.MAD") { 
        c1 <- switch(tp2, 5.01, 7.52, 10.03) 
        rt <- Tirls.mad(x1, y1, rt, c1, rp.max, cg.rt)
        cd <- paste("wf=", wf, "scale=", scale, "tp=", tp, "c1=", c1, sep=",")
        return(list(cond=cd, TK=rt$TK, wt=rt$wt, rp=rt$rp, s1=rt$s1))
     } else if (wfs == "H.AAD") { 
        c1 <- switch(tp2, 1.15, 1.72, 2.30)
        rt <- Hirls.aad(x1, y1, rt, c1, rp.max, cg.rt)
        cd <- paste("wf=", wf, "scale=", scale, "tp=", tp, "c1=", c1, sep=",")
        return(list(cond=cd, HB=rt$HB, wt=rt$wt, rp=rt$rp, s1=rt$s1))
     } else if (wfs == "H.MAD") { 
        c1 <- switch(tp2, 1.44, 2.16, 2.88)
        rt <- Hirls.mad(x1, y1, rt, c1, rp.max, cg.rt)
        cd <- paste("wf=", wf, "scale=", scale, "tp=", tp, "c1=", c1, sep=",")
        return(list(cond=cd, HB=rt$HB, wt=rt$wt, rp=rt$rp, s1=rt$s1))
     } else stop('Select an "T" or "H" for wf, and "AAD" or "MAD" for scale ')

}
