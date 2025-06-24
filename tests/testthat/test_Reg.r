
set.seed(4)
 cov1 <- matrix(c(3, 2.8, 2.8, 3), 2, 2)
 cov2 <- matrix(c(2.5, 0, 0, 3), 2, 2)

 dat1 <- MASS::mvrnorm(n=400, mu=c(100, 100), Sigma=cov1, empirical=TRUE)
 dat2 <- cbind(runif(100,  min=96, max=104),  runif(50,  min=95, max=105))
 dat3 <- matrix(c(103, 103.5, 104.5, 104.8, 96, 98, 94, 95), 4, 2) 

dat <- rbind(dat1, dat2, dat3)
plot(dat)

y1 <- dat[,2]
x1 <- dat[,1]

o1 <- robReg(x1, y1)
o2 <- lm(y1~x1)

    fg1 <- rep(1, length(y1))
    fg1[which(o1$wt < 0.8)] 	<- 3
    fg1[which(o1$wt < 0.5)] 	<- 7
    fg1[which(o1$wt < 0.2)] 	<- 2
    fg1[which(o1$wt == 0)]	<- 8

plot(x1, y1, pch=19, col=fg1)
 abline(a=o1$TK$coeff[1], b=o1$TK$coeff[2], col="green", lwd=3)
 abline(a=o2$coeff[1], b=o2$coeff[2], col="red", lwd=3, lty=3)
