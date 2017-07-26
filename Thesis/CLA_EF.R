assets50 <- GetIndex(1:50, assets1457, lB = rep(0, 50), uB = rep(1, 50))
r50.7.1 <- Env7.1$cla.solve(assets50, purge = FALSE, lam.tol = 1e-14)
ms50.7 <- Env7$cla.solve(assets50, purge= TRUE)$MS_weight
m <- nrow(ms50.7)
ms50.7.1 <- MS.insert(r50.7.1$weights_set, assets50$mu, assets50$covar)

# sample some points, first decide sample range
set.seed(2017)
w.r <- matrix(0, nrow = 50, ncol = 10000)
for(j in 1:10000){
  ran <- runif(50, 0, 1)
  w.r[, j] <- ran/sum(ran)
}
ms.r <- MS(w.r, assets50$mu, assets50$covar)
mu.min <- min(ms.r[, 2])
sig.max <- max(ms.r[ ,1])
x <- runif(100, 0.015, sig.max)
y <- runif(100, mu.min, 0.002)
######################################## not used... ##########
pdf("D:/Master Thesis/MyThesis/Pictures/EFeg.pdf", width = 10, height = 7) 
plot(ms50.7.1, type = "l", xaxt = "n", yaxt = "n",  
     main = "Efficient Frontier", ylab = expression(mu(w)),
     xlab = expression(sigma(w)))
points(ms50.7, pch = 16, col = c(rep("blue", m-1), "red"))
legend("bottomright", legend = c("turning points", "minmum variance portfolio"), 
       pch = 16, col = c("blue", "red"))
dev.off()

################################################################
# only changing upper bounds
pdf("D:/Master Thesis/MyThesis/Pictures/EF_u.pdf", width = 10, height = 5) 
plot_EF_bound(l = rep(0, 6), 
              u = c(1, 0.5, 0.25, 0.125, 0.0625, 0.03125), ind = 1:50)
dev.off()

# only changing lower bounds
pdf("D:/Master Thesis/MyThesis/Pictures/EF_l.pdf", width = 10, height = 5) 
plot_EF_bound(l = seq(0, 0.01, by = 0.002), 
              u = rep(1, 6), ind = 1:50)
dev.off()





