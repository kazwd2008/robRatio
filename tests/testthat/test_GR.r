#---------------------------------------
set.seed(105)

x1 <- rnorm(100, mean=100, sd=1) 	
ya100 <- yb100 <- rep(NA, 100)
 
e2 <- rnorm(100, sd=0.5)
ya100 <- 2 * x1 + e2 * x1^0.5         
yb100 <- 2 * x1 + e2 * x1         
o1a <- robGR(x1, ya100)
o1b <- robGR(x1, yb100)

o1a$par; o1a$g1
o1b$par; o1b$g1

#---------------------------------------
# set.seed(105)
# 
# x1 <- rnorm(100, mean=100, sd=1) 	
# ya100 <- yb100 <- matrix(NA, nr=1000, nc=100)
# ot1a <- ot1b  <- matrix(NA, nr=1000, nc=2)	
# f.ot1a <- f.ot1b <- rep(NA, 1000)
# 
# for (i in 1:1000){
#     e2 <- rnorm(100, sd=0.5)
#     ya100[i,] <- 2 * x1 + e2 * x1^0.5         
#     yb100[i,] <- 2 * x1 + e2 * x1         
#     o1a <- robGR(x1, ya100[i,])
#     o1b <- robGR(x1, yb100[i,])
#     ot1a[i,] <- c(o1a$par, o1a$g1)
#     ot1b[i,] <- c(o1b$par, o1b$g1)
#     f.ot1a[i] <- o1a$efg
#     f.ot1b[i] <- o1b$efg
# }
# 
# par(mfrow=c(2,2))
#   plot(ot1a[,1], pch=21, bg="gray", main=latex2exp::TeX(r"(data A: estimated $\beta$)"))
#    abline(a=2, b=0, col="red")
# 
#   plot(ot1a[,2], pch=21, bg="gray", main=latex2exp::TeX(r"(data A: estimated $\gamma$)"))
#    abline(a=0.5, b=0, col="red")
#    abline(a=apply(ot1a, 2, mean)[2], b=0, col="blue", lty=3, lwd=2)
#    abline(a=apply(ot1a, 2, median)[2], b=0, col="cyan", lty=3, lwd=2)
# 
#  plot(ot1b[,1], pch=21, bg="gray", main=latex2exp::TeX(r"(data B: estimated $\beta$)"))
#    abline(a=2, b=0, col="red")
# 
#  plot(ot1b[,2], pch=21, bg="gray", main=latex2exp::TeX(r"(data B: estimated $\gamma$)"))
#    abline(a=1, b=0, col="red")
#    abline(a=apply(ot1b, 2, mean)[2], b=0, col="blue", lty=3, lwd=2)
#    abline(a=apply(ot1b, 2, median)[2], b=0, col="cyan", lty=3, lwd=2)
# 
# apply(ot1a, 2, mean)
# apply(ot1a, 2, median)
# apply(ot1b, 2, mean)
# apply(ot1b, 2, median)
#----------------------------------------------