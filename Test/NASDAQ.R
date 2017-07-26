if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
########################################################################
# NASDAQ
index.df["NASDAQ",] # 2196 assets, 265 observations
analysis.nasdaq <- IndexAnalysis(Index.trans = transIndex(NASDAQ), 
                                 name = "nasdaq", n = ncol(NASDAQ))

#######################################################################################
# plot
IndexPlot(analysis.nasdaq, name = "nasdaq", n = ncol(NASDAQ))  # open nasdaq_plot.pdf


#########################################################################################

# time of qp dropped sharply at 285, why?? SP500?
which(label.nasdaq == 285) #10
analysis.nasdaq$lambda.QP[[10]]  # numeric(0)
nasdaq <- GetIndex(ind.nasdaq[[10]][[1]], NASDAQ)
QP.solve(nasdaq$mu, nasdaq$covar, nasdaq$lB, nasdaq$uB, 0.3)
# error1: D(i.e covariance matrix not positive definite)
analysis.nasdaq$cov.check # is.positive.definite = FALSE for last 3 columns

if(!file.exists(d.file("near_check_NASDAQ.rds", exists = FALSE))){
  set.seed(2017)
  l <- list()
  k <- 10
  near.check.nasdaq <- matrix(0, ncol = 10, nrow = k)
  for(j in 1:10){
    l[[j]] <- lapply(1:k, function(x) sample(seq(n.nasdaq), j+260)) # 251-270 assets
    near.check.nasdaq[, j] <- sapply(1:k, function(x) 
      is.positive.definite(GetIndex(l[[j]][[x]], NASDAQ)$covar))
  }
  saveRDS(near.check.nasdaq, d.file("near_check_NASDAQ.rds", exists = FALSE))
}
near.check.nasdaq <- readRDS(d.file("near_check_NASDAQ.rds", exists = FALSE))
near.check.nasdaq # continuous?




# t_compare_cla_qp
plot(seq(l.nasdaq), t_mean[1,], type = "o", ylim = range(t_mean), col = "red", 
     pch = 16, xaxt = "n", xlab = "Number of Assets", 
     ylab = "milliseconds", main = "Time Comparison of CLA and QP")
lines(seq(l.nasdaq), t_mean[2,], col = "blue", pch = 16, type = "o")
axis(1, at = seq(l.nasdaq),labels = label.nasdaq)
legend("topleft", legend = c("cla", "qp"), col = c("red", "blue"), lwd = 1, pch = 16)

# qp method is faster than cla when the number of assets is small, (55-95)
# but much slower when number of assets is larger (114-410)

nlambda.QP <- matrix(0, ncol = l.nasdaq, nrow = 10)
for(i in seq(l.nasdaq)){
  nlambda.QP[,i] <- sapply(1:10, function(x) length(analysis.nasdaq$lambda.QP[[i]][[x]]))
}
colnames(nlambda.QP) <- label.nasdaq 
nlambda.QP
#eg. ind.nasdaq[[1]][[1]]
nw <- analysis.nasdaq$nweights[1,1]
lam.err <- setdiff(((1:nw)/nw)^2, lambda.nasdaq[[1]][[1]])
nasdaq <- GetIndex(ind.nasdaq[[1]][[1]], NASDAQ)
QP.solve(nasdaq$mu, nasdaq$covar, nasdaq$lB, nasdaq$uB, lam.err[1]) 
# error 2: constraints are inconsistent, no solution


