####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####
#' 
#' @title Robust estimator for ratio models
#'
#' @description This robRatio function integrates 12 functions below for easy 
#' use for users.  Please note that the values for the tuning parameter \code{tp} 
#' allowed in this function is standardized. 
#
####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|
#' @name robRatio
#' @param x1 single explanatory variable (a vector)
#' @param y1 objective variable to be imputed (a vector)
#' @param gm indication of gamma value as follows:  \cr
#'           gm="a": gamma=1 \cr
#'           gm="b": gamma=1/2	(conventional ratio model) \cr
#'           gm="c"; gamma=0    (regression model without intercept) \cr
#' @param wf weight function (wf=T : Tukey, wf=H : Huber) 
#' @param scale scale for residuals. "AAD"(default) or "MAD".
#' @param rt sample weight (default 1)
#' @param tp standardized tuning parameter. choose 4, 6 or 8. Smaller figure is 
#'           more robust (default tp=8).  See details.
##  @param T.tp tuning parameter for Tukey's biweight function. See details.
##  @param H.tp tuning parameter for Huber weight function. See details. 
#' @param rp.max maximum number of iteration (default: rp.max=50)
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#' @return a list with the following elements
#' \describe{
#'   \item{\code{cond}}{Weight function, scale, and other arguments choosed}
#'   \item{\code{par}}{robustly estimated ratio of y1 to x1 (beta)}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes of the scale (AAD or MAD)}
#'   \item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
#' }
####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|
#'
#' @examples
#' \dontrun{
#' require(robRatio)
#'
#' x1 <- seq(1, 10, by=0.1)
#' #e <- rnorm(length(x))
#' e <- rt(length(x1), df=3)   # error term following t distribution
#'
#' b <- 2		# true value of slope
#'
#' y1 <- b*x + x*e			# example 1: gamma=1
#' y2 <- b*x + sqrt(x)*e   # example 2: gamma=1/2
#'
#' o1 <- robRatio(x1, y1, gm="a")
#' o2 <- robRatio(x1, y2, gm="b")
#'
#' o1$par;  o2$par     # estimated slope
#'
#' cols = RColorBrewer::brewer.pal(11, "PiYG")
#' cl1 <- round((o1$wt)*10+1)
#' cl2 <- round((o2$wt)*10+1)
#'
#' par(mfrow=c(1,2))
#' plot(x, y1, col=cols[cl1], pch=20)
#' plot(x, y2, col=cols[cl2], pch=20)
#' }
#'
#' @export
####+####|####+####|####+####|####+####|####+####|####+####|####+####|####+####|

robRatio <- function(x1, y1, gm="b", wf="T", scale="AAD", rt=1, tp=8, 
                     rp.max=100, cg.rt=0.01){
#----------------------------------------------- value check arguments
  vl.gm    <- c("a", "b", "c")
  vl.wf    <- c("T", "H")
  vl.scale <- c("AAD", "MAD")
  vl.tp    <- c(4, 6, 8)

  if (gm    %in% vl.gm    == FALSE) stop("Please choose 'a', 'b' or 'c' for gm")
  if (wf    %in% vl.wf    == FALSE) stop("Please choose 'T' or 'H' for wf")
  if (scale %in% vl.scale == FALSE) stop("Please choose 'AAD' or 'MAD' for scale")
  if (tp    %in% vl.tp    == FALSE) stop("Please choose 4, 6 or 8 for tp")
  if (length(rt) %in% c(1, length(y1)) ==FALSE) stop("Incorrect length of rt")

#----------------------------------------------- select an appropriate function

 tp2 <- tp/2-1  # (4, 6, 8) => (1, 2, 3)

  if (wf=="T") {

     if (scale=="AAD"){

        c1 <- switch(tp2, 4, 6, 8)
        if (gm=="a") ot <- RrT.aad(x1, y1, g1=1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrT.aad(x1, y1, g1=1/2, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrT.aad(x1, y1, g1=0, c1, rp.max, cg.rt)

     } else {

        c1 <- switch(tp2, 5.01, 7.52, 10.03)  # tuning constant for SD (MAD in R)
        if (gm=="a") ot <- RrT.mad(x1, y1, g1=1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrT.mad(x1, y1, g1=1/2, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrT.mad(x1, y1, g1=0, c1, rp.max, cg.rt)

     }
     cd <- paste("gm=",gm, "scale=",scale,"wf=",wf, "tp=", tp, "c1=",c1, sep=",") 
     return(list(cond=cd, par=ot$par, res=ot$res, wt=ot$wt, rp=ot$rp, s1=ot$s1, efg=ot$efg))

  } else {
     if (scale=="AAD"){

        c1 <- switch(tp2, 1.15, 1.72, 2.30)
        if (gm=="a") ot <- RrH.aad(x1, y1, g1=1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrH.aad(x1, y1, g1=1/2, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrH.aad(x1, y1, g1=0, c1, rp.max, cg.rt)

     } else {

        c1 <- switch(tp2, 1.44, 2.16, 2.88)  # tuning constant for SD (MAD in R)
        if (gm=="a") ot <- RrH.mad(x1, y1, g1=1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrH.mad(x1, y1, g1=1/2, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrH.mad(x1, y1, g1=0, c1, rp.max, cg.rt)

     }
     cd <- paste("gm=",gm, "scale=",scale,"wf=",wf, "tp=", tp, "c1=",c1, sep=",") 
     return(list(cond=cd, par=ot$par, res=ot$res, wt=ot$wt, rp=ot$rp, s1=ot$s1, efg=ot$efg))
  }
}