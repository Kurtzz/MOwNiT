pdf(file="result.pdf", height=8, width=8, bg="white")

plot(function(x) x^2, -5, 5)
plot(sin, -4, 4)

res = read.csv("f_osc.txt")

which = which(res$alg=="qag")
alg = res[which,]
plot(alg$a, alg$n, type="l")

which = which(res$alg=="qags")
alg = res[which,]
plot(alg$a, alg$n, type="l")

which = which(res$alg=="qagp")
alg = res[which,]
plot(alg$a, alg$n, type="l")

which = which(res$alg=="qawo")
alg = res[which,]
plot(alg$a, alg$n, type="l")

dev.off()
