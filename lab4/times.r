library("Hmisc")
results = read.csv("times.txt")

avg_results = aggregate( time ~ n:alg, data=results, FUN=mean)
avg_results$sd = aggregate( time ~ n:alg, data=results, FUN=sd)$time

y_range = range(0,avg_results$time)
x_range = range(0,avg_results$n)

pdf(file="times.pdf", height=8, width=8, bg="white")


plot(x=avg_results[avg_results$alg=="Lagrange",]$n,y=avg_results[avg_results$alg=="Lagrange",]$time, type="l", lty = 1, ylim = y_range,xlim = x_range, xlab = "n", ylab = "Computing time in seconds", col="blue")
lines(x=avg_results[avg_results$alg=="Newton",]$n,y=avg_results[avg_results$alg=="Newton",]$time, lty = 1, col = "red")
lines(x=avg_results[avg_results$alg=="GSL_Polym_Interpol",]$n,y=avg_results[avg_results$alg=="GSL_Polym_Interpol",]$time, lty = 1, col = "green")
lines(x=avg_results[avg_results$alg=="Spline",]$n,y=avg_results[avg_results$alg=="Spline",]$time, type="l", lty = 1, col = "cadetblue1")
lines(x=avg_results[avg_results$alg=="Akima",]$n,y=avg_results[avg_results$alg=="Akima",]$time, type="l", lty = 1, col = "darkorange1")
legend("topleft", c("Lagrange","Newton","GSL_Polym_Interpol","Spline","Akima"),col=c("blue","red","green","cadetblue1","darkorange1"),lty=1)

plot(x=avg_results[avg_results$alg=="Lagrange",]$n,y=avg_results[avg_results$alg=="Lagrange",]$time, type="l", lty = 1, xlab = "n", ylab = "Computing time in seconds", col="blue")
legend("topleft", c("Lagrange"),col=c("blue"),lty=1)

plot(x=avg_results[avg_results$alg=="Newton",]$n,y=avg_results[avg_results$alg=="Newton",]$time, type="l", lty = 1, xlab = "n", ylab = "Computing time in seconds", col="red")
legend("topleft", c("Newton"),col=c("red"),lty=1)

plot(x=avg_results[avg_results$alg=="GSL_Polym_Interpol",]$n,y=avg_results[avg_results$alg=="GSL_Polym_Interpol",]$time, type="l", lty = 1, xlab = "n", ylab = "Computing time in seconds", col="green")
legend("topleft", c("GSL_Polym_Interpol"),col=c("green"),lty=1)

plot(x=avg_results[avg_results$alg=="Spline",]$n,y=avg_results[avg_results$alg=="Spline",]$time, type="l", lty = 1, xlab = "n", ylab = "Computing time in seconds", col="cadetblue")
legend("topleft", c("Spline"),col=c("cadetblue"),lty=1)

plot(x=avg_results[avg_results$alg=="Akima",]$n,y=avg_results[avg_results$alg=="Akima",]$time, type="l", lty = 1, xlab = "n", ylab = "Computing time in seconds", col="darkorange")
legend("topleft", c("Akima"),col=c("darkorange"),lty=1)

errbar( avg_results[avg_results$alg=="Lagrange",]$n, avg_results[avg_results$alg=="Lagrange",]$time, avg_results[avg_results$alg=="Lagrange",]$time + avg_results[avg_results$alg=="Lagrange",]$sd, avg_results[avg_results$alg=="Lagrange",]$time - avg_results[avg_results$alg=="Lagrange",]$sd, xlab = "n", ylab = "Computing time in seconds",type="o", lty = 1,col = "red")
legend("topleft", c("Lagrange"),col=c("red"),lty=1)

errbar( avg_results[avg_results$alg=="Newton",]$n, avg_results[avg_results$alg=="Newton",]$time, avg_results[avg_results$alg=="Newton",]$time + avg_results[avg_results$alg=="Newton",]$sd, avg_results[avg_results$alg=="Newton",]$time - avg_results[avg_results$alg=="Newton",]$sd, xlab = "n", ylab = "Computing time in seconds",type="o", lty = 1,col = "green")
legend("topleft", c("Newton"),col=c("green"),lty=1)

errbar( avg_results[avg_results$alg=="Spline",]$n, avg_results[avg_results$alg=="Spline",]$time, avg_results[avg_results$alg=="Spline",]$time + avg_results[avg_results$alg=="Spline",]$sd, avg_results[avg_results$alg=="Spline",]$time - avg_results[avg_results$alg=="Spline",]$sd, xlab = "n", ylab = "Computing time in seconds",type="o", lty = 1,col = "cadetblue")
legend("topleft", c("Spline"),col=c("cadetblue"),lty=1)

errbar( avg_results[avg_results$alg=="Akima",]$n, avg_results[avg_results$alg=="Akima",]$time, avg_results[avg_results$alg=="Akima",]$time + avg_results[avg_results$alg=="Akima",]$sd, avg_results[avg_results$alg=="Akima",]$time - avg_results[avg_results$alg=="Akima",]$sd, xlab = "n", ylab = "Computing time in seconds",type="o", lty = 1,col = "darkorange")
legend("topleft", c("Akima"),col=c("darkorange"),lty=1)



dev.off()
