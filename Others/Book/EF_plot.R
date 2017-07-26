if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")


##########################################################################################
ind <- SP500[,1:20]  # define a index

Nassets <- ncol(ind)  
Nobs <- nrow(ind)
mu <- colMeans(ind)
S <- cov(ind)
SR <- sqrm(S) # square root of a matrix
delta <- sqrt(qchisq(0.9, Nassets)) # 90% quantile for the returns uncertainties, under ND
## Determining feasible risk aversion
SigMax <- max(colSds(ind))
SigMin <- min(colSds(ind))
#in case of overfitting
ra <- seq(SigMin * 1.1, SigMax * 0.9, length.out = 20) / SigMax 
SigSeq <- SigMin/ra


##list of risks, outlier removed and rescaled

## classical mean-variance(MV) and robust counterpart(RC) results
MV <- sapply(SigSeq, function(x) c(Porgr(SR, mu, x)$risk, Porgr(SR, mu, x)$mean))
rownames(MV) <- c("risk", "mean")


theta <- ra + (1 - ra) * delta / sqrt(Nobs) 
Robust_MV <- sapply(theta, function(x)c(Porgr(SR, mu, SigMin/x)$risk, Porgr(SR, mu, SigMin/x)$mean) )
rownames(Robust_MV) <- c("risk", "mean")

## replace theta by lambda
theta2lambda <- function(theta, Nassets, Nobs, level){
  delta <- sqrt(qchisq(level, Nassets))
  lambda <- theta / (1 + theta * (delta / sqrt(Nobs)))
  lambda
}
## robust allocation and equivalent mean-variance allocation, eg. theta = 0.7


theta0 <- 0.7 
MV0 <- c(Porgr(SR, mu, SigMin / theta0)$risk, Porgr(SR, mu, SigMin / theta0)$mean)

lam <- theta2lambda(theta0, Nassets = Nassets, Nobs = Nobs, level = 0.9)
Robust_MV0 <- c(Porgr(SR, mu = mu, SigMin / lam)$risk, Porgr(SR, mu = mu, SigMin / lam)$mean)

## plotting efficient frontier 
plot(MV["risk", ], MV["mean", ], type = "o",
     xlim = range(MV["risk", ], Robust_MV["risk", ]), 
     ylim = range(MV["mean", ], Robust_MV["mean", ]),
     xlab = expression(sigma), ylab = expression(mu), pch =16, col = "blue")
points(Robust_MV["risk", ], Robust_MV["mean", ], col = "red", pch = 19)

## Superimposing equivalence points
points(MV0[1], MV0[2], col = "green",  pch = 16)
points(Robust_MV0[1], Robust_MV0[2], col = "darkgreen", pch = 16)
## Legend
legend("bottomright", legend = c("EF", "MV", "Robusts MV",
                             expression(paste("Robust MV: ", theta == 0.7)),
                             "Equivalent MV"),
       lty = c(1, NA, NA, NA, NA), pch = c(NA, 16, 16, 16, 16),
       col = c("blue", "blue", "red", "darkgreen", "green"))


source("Functions/Ver1_previous.R")
#previous purged version, w>0
wp <- as.matrix(CLA$M1(assets50$mu, assets50$covar, 
                       as.matrix(rep(0,50)), as.matrix(rep(1,50)))[[3]])

#previous purged version, 1e-8 < w < 0.075
wplu <- as.matrix(CLA$M1(assets50$mu, assets50$covar, 
                                 assets50$lB, assets50$uB)[[3]])

#my purged version, 1e-8 < w < 0.075 
wmlu <- as.matrix(CLA$M3(assets50$mu, assets50$covar, 
                                  assets50$lB, assets50$uB)$weights_set)
#my purged version, 1e-8 < w  
wml <- as.matrix(CLA$M3(assets50$mu, assets50$covar, 
                            assets50$lB, as.matrix(rep(1,50)))$weights_set)

sig50 <- sqrm(assets50$covar)
isSymmetric(sig50)

mu50 <- assets50$mu


wq <- sapply(sqrt(diag(t(wmlu) %*% assets50$covar %*% wmlu)),
             function(x) Porgr(sig50, as.vector(assets50$mu), x)$weight)
colSums(wq) # good, always 1
# all weights are non-negative, not practical

MS.q <- MS(wq, assets50$mu, assets50$covar)

###

MS.plu <- rbind(sqrt(diag(t(wplu) %*% assets50$covar %*% wplu)),
                t(mu50) %*% wplu)

MS.ml <- rbind(sqrt(diag(t(wml) %*% assets50$covar %*% wml)),
                t(mu50) %*% wml)

MS.p <- rbind(sqrt(diag(t(wp) %*% assets50$covar %*% wp)),
               t(mu50) %*% wp)

plot(MS.q[,"Sig"], MS0[, "Mu"], type = "o", pch =16, col = "blue", 
     xlab = expression(sigma), ylab = expression(mu), 
     xlim = range(cbind(MS.p, MS.mlu, MS.plu, MS.ml, t(MS.q)[2:3,])[1,]), 
     ylim = range(cbind(MS.p, MS.mlu, MS.plu, MS.ml, t(MS.q)[2:3,])[2,]))

points(MS.p[1,], MS.p[2,], col = "orange", pch = 16)
points(MS.mlu[1,], MS.mlu[2,], col = "red", type = "o")
points(MS.plu[1,], MS.plu[2,], col = "green")
points(MS.ml[1,], MS.ml[2,], col = "yellow", pch=16)

legend("bottomright", 
       legend = c("MV", "cla-p", "cla-mlu", "cla-plu", "cla-mlu-chull", "cla-ml"),
       lty = 1, col = c("blue", "orange", "red", "green", "purple", "yellow"))

#the weight set of highest expected return, absolutely true?!
w.lu.max <- initAlgo(assets50$mu, assets50$lB, assets50$uB)$weights
points(sqrt(diag(t(w.lu.max) %*% assets50$covar %*% w.lu.max)),
       t(mu50) %*% w.lu.max, col = "black", pch = 16)


w.max <-initAlgo(assets50$mu, as.matrix(rep(0,50)), as.matrix(rep(1,50)))$weights
points(sqrt(diag(t(w.max) %*% assets50$covar %*% w.max)),
       t(mu50) %*% w.max, col = "black", pch = 16)

#analysis:
# 1) cla-plu and cla-mlu:
# previous version and my version come up with almost same weights set (almost same on the graph).
# the little differences might be caused by improvement of matrix inverse 
# (computeLambda and computeW)

#

# 3) cla-p and cla-ml:
# lower bounds make negligible difference (but upper bounds do: cla-ml vs cla-mlu)

# 4) MV
# part of MV efficient frontier






