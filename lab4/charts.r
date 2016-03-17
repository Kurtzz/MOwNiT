results = read.csv("result.txt")

avg_results = aggregate( cbind(x,y) ~ x:alg, data=results, FUN=mean)

x_range = range(0,avg_results$x+1)
y_range = range(min(avg_results$y),max(avg_results$y))

pdf(file="charts.pdf", height=8, width=8, bg="white")

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="GSL",]$y~avg_results[avg_results$alg=="GSL",]$x,col="red")
legend("topleft", c("GSL-Interpolation"),col=c("red"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="Lagrange",]$y~avg_results[avg_results$alg=="Lagrange",]$x,col="blue")
legend("topleft", c("Lagrange-Interpolation"),col=c("blue"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y" )
lines(avg_results[avg_results$alg=="Newton",]$y~avg_results[avg_results$alg=="Newton",]$x,col="green")
legend("topleft", c("Newton-Interpolation"),col=c("green"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="GSL",]$y~avg_results[avg_results$alg=="GSL",]$x,col="red")
lines(avg_results[avg_results$alg=="Lagrange",]$y~avg_results[avg_results$alg=="Lagrange",]$x,col="blue")
lines(avg_results[avg_results$alg=="Newton",]$y~avg_results[avg_results$alg=="Newton",]$x,col="green")
legend("topleft", c("GSL-Interpolation","Lagrange-Interpolation","Newton-Interpolation"),col=c("red","blue","green"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="Spline",]$y~avg_results[avg_results$alg=="Spline",]$x,col="dodgerblue")
legend("topleft", c("Spline"),col=c("dodgerblue"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="Akima",]$y~avg_results[avg_results$alg=="Akima",]$x,col="forestgreen")
legend("topleft", c("Akima"),col=c("forestgreen"),lty=1)

x_range = range(0,avg_results$x+1)
y_range = range(0,(avg_results[avg_results$alg=="Gen_points",]$y))

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="GSL",]$y~avg_results[avg_results$alg=="GSL",]$x,col="red")
legend("topleft", c("GSL-Interpolation"),col=c("red"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y")
lines(avg_results[avg_results$alg=="Lagrange",]$y~avg_results[avg_results$alg=="Lagrange",]$x,col="blue")
legend("topleft", c("Lagrange-Interpolation"),col=c("blue"),lty=1)

plot(avg_results[avg_results$alg=="Gen_points",]$y~avg_results[avg_results$alg=="Gen_points",]$x, ylim = y_range,xlim = x_range,xlab = "AXIS X", ylab = "AXIS Y" )
lines(avg_results[avg_results$alg=="Newton",]$y~avg_results[avg_results$alg=="Newton",]$x,col="green")
legend("topleft", c("Newton-Interpolation"),col=c("green"),lty=1)




dev.off()
