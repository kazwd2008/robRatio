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
robReg <- function(x1, y1, wf="T", scale="AAD", rt=rep(1, length(y1)), tp=8, 
                   rp.max=150, cg.rt=0.01) {

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
