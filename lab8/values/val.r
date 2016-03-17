pdf(file="vals.pdf", height=8, width=8, bg="white")

res = read.csv("values_10_0.5_50.txt")
barplot(res$fft_val, main="FREQ = 10, AMP = 0.5, N = 50")

res = read.csv("values_10_0.5_200.txt")
barplot(res$fft_val, main="FREQ = 10, AMP = 0.5, N = 200")

res = read.csv("values_20_0.5_200.txt")
barplot(res$fft_val, main="FREQ = 20, AMP = 0.5, N = 200")

res = read.csv("values_10_1_200.txt")
barplot(res$fft_val, main="FREQ = 10, AMP = 1, N = 200")


res = read.csv("values_fft_inv.txt")
barplot(res$fun, main="Wykres z zakloceniami")
barplot(res$fft_val, main="wykres transformaty Fouriera")
barplot(res$inv_fft, main="Wykres po transformacie odwrotnej")

res = read.csv("values_fun2.txt")
barplot(res$fun)
barplot(res$fft_val)

dev.off()
