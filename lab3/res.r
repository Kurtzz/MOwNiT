results = read.csv("result.txt")
avg_results = aggregate( time ~ n:alg, data=results, FUN=mean)
avg_results$sd = aggregate( time ~ n:alg, data=results, FUN=sd)$time

y_range = range(0, results$time)


pdf(file="res.pdf", height=5, width=10, bg="white")

plot(avg_results[avg_results$alg=="alg1",]$n, avg_results[avg_results$alg=="alg1",]$time, lty = 1, ylim = y_range, xlab = "n", ylab = "Computing time in seconds", col="blue")
lines(avg_results[avg_results$alg=="alg2",]$n, avg_results[avg_results$alg=="alg2",]$time, lty = 2, col = "red")
lines(avg_results[avg_results$alg=="blas",]$n, avg_results[avg_results$alg=="blas",]$time, lty = 3, col = "green")
legend("topleft", c("alg1","alg2", "blas"),col=c("blue","red", "green"),lty=1:3)


errbar(avg_results[avg_results$alg=="alg1",]$n,
	avg_results[avg_results$alg=="alg1",]$time,
	avg_results[avg_results$alg=="alg1",]$time-avg_results[avg_results$alg=="alg1",]$sd,
	avg_results[avg_results$alg=="alg1",]$time+avg_results[avg_results$alg=="alg1",]$sd,
	xlab="n", ylab="Computing time in seconds", add=TRUE)


errbar(avg_results[avg_results$alg=="alg2",]$n,
	avg_results[avg_results$alg=="alg2",]$time,
	avg_results[avg_results$alg=="alg2",]$time-avg_results[avg_results$alg=="alg2",]$sd,
	avg_results[avg_results$alg=="alg2",]$time+avg_results[avg_results$alg=="alg2",]$sd,
	xlab="n", ylab="Computing time in seconds", add=TRUE)

errbar(avg_results[avg_results$alg=="blas",]$n,
	avg_results[avg_results$alg=="blas",]$time,
	avg_results[avg_results$alg=="blas",]$time-avg_results[avg_results$alg=="blas",]$sd,
	avg_results[avg_results$alg=="blas",]$time+avg_results[avg_results$alg=="blas",]$sd,
	xlab="n", ylab="Computing time in seconds", add=TRUE)

x1 = avg_results[avg_results$alg=="alg1",]$n
y1 = avg_results[avg_results$alg=="alg1",]$time
fit1 <- lm(y1 ~ poly(x1, 3, raw=TRUE))
xx1 <- seq(100,700, length.out=250)
lines(xx1, predict(fit1, data.frame(x1=xx1)), col='blue')

x2 = avg_results[avg_results$alg=="alg2",]$n
y2 = avg_results[avg_results$alg=="alg2",]$time
fit2 <- lm(y2 ~ poly(x2, 3, raw=TRUE))
xx2 <- seq(100,700, length.out=250)
lines(xx2, predict(fit2, data.frame(x2=xx2)), col='red')

xb = avg_results[avg_results$alg=="blas",]$n
yb = avg_results[avg_results$alg=="blas",]$time
fitb <- lm(yb ~ poly(xb, 3, raw=TRUE))
xxb <- seq(100,700, length.out=250)
lines(xxb, predict(fitb, data.frame(xb=xxb)), col='green')

dev.off()
