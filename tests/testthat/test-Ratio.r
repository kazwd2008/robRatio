 set.seed(453)
 x1 <- rnorm(500,  mean=50, sd=5)
 e1 <- rt(500, df=3)
 r1 <- 2
 y1 <- r1 * x1 + e1 * sqrt(x1)

 out <- which(x1 > 55)
 y1[sample(out, 10)] <- y1[sample(out, 10)] * 2

 o1.Tk8 <- robRatio(x1, y1)
 o2 <- sum(y1)/sum(x1)

    fg1 <- rep(1, length(y1))
    fg1[which(o1.Tk8$wt < 0.8)] 	<- 3
    fg1[which(o1.Tk8$wt < 0.5)] 	<- 7
    fg1[which(o1.Tk8$wt < 0.2)] 	<- 2
    fg1[which(o1.Tk8$wt == 0)]	<- 8

 plot(x1, y1, pch=19, col=fg1, main="test robRatio")
    abline(a=0, b=o1.Tk8$par, col="green", lwd=3)
    abline(a=0, b=o2, col="red", lty=3, lwd=3)
