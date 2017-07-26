# company asset data
pdf("D:/Master Thesis/MyThesis/Pictures/company.pdf",
    width = 8, height = 5) 
plot(sqrt(diag(assets1457$covar)), assets1457$mu, pch = 16, cex = 0.4,
     xlab = expression(sigma), ylab = expression(mu), 
     main = "Company Assets")
dev.off()
sum(assets1457$mu > 0)/1457
quantile(sqrt(diag(assets1457$covar)), c(0.1, 0.9))
quantile(assets1457$mu, c(0.1, 0.9))

# SP500
se <- sqrt(diag(trans.sp$covar))
i <- order(se, decreasing = TRUE)[1:2] #remove the outlier for better ploting..
pdf("D:/Master Thesis/MyThesis/Pictures/sp500.pdf",
    width = 8, height = 5)
plot(se[-i], trans.sp$mu[-i], pch = 16, cex = 0.4,
     xlab = expression(sigma), ylab = expression(mu), 
     main = "SP500")
dev.off()

quantile(trans.sp$mu, c(0.1, 0.9))
quantile(sqrt(diag(trans.sp$covar)), c(0.1, 0.9))

# NASDAQ
se <- sqrt(diag(trans.nasdaq$covar))
i <- order(se, decreasing = TRUE)[1:5] #remove the outlier for better ploting..
pdf("D:/Master Thesis/MyThesis/Pictures/nasdaq.pdf",
    width = 8, height = 5) 
plot(se[-i], trans.nasdaq$mu[-i], pch = 16, cex = 0.4,
     xlab = expression(sigma), ylab = expression(mu), 
     main = "NASDAQ")
dev.off()

quantile(trans.nasdaq$mu, c(0.1, 0.9))
quantile(sqrt(diag(trans.nasdaq$covar)), c(0.1, 0.9))
