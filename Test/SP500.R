if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
########################################################################
# SP500
index.df["SP500",] # 476 assets, 265 observations
analysis.sp <- IndexAnalysis(Index.trans = transIndex(SP500), 
                             name = "sp", n = ncol(SP500))

#######################################################################################
# plot
IndexPlot(analysis.sp, name = "sp", ncol(SP500))  # open sp_plot.pdf










#######################################################################################
  # time of qp dropped sharply at 285, why?? 
which(label.sp == 285) #10
analysis.sp$lambda.QP[[10]]  # numeric(0)
sp <- GetIndex(ind.sp[[10]][[1]], SP500)
QP.solve(sp$mu, sp$covar, sp$lB, sp$uB, 0.3)
# error1: D(i.e covariance matrix not positive definite)
analysis.sp$cov.check # is.positive.definite = FALSE for last 3 columns

if(!file.exists(d.file("near_check_SP.rds"))){
  set.seed(2017)
  l <- list()
  near.check.sp <- matrix(0, ncol = 10, nrow = 20)
  for(j in 1:10){
    l[[j]] <- lapply(1:20, function(x) sample(seq(n.sp), j+260)) # 261-270 assets
    near.check.sp[, j] <- sapply(1:20, function(x) 
      is.positive.definite(GetIndex(l[[j]][[x]], SP500)$covar))
  }
  colnames(near.check.sp) <- 261:270
  saveRDS(near.check.sp, d.file("near_check_SP.rds", exists = FALSE))
}
near.check.sp <- readRDS(d.file("near_check_SP.rds", exists = FALSE))
near.check.sp # continuous?

# t_compare_cla_qp


# qp method is faster than cla when the number of assets is small, (55-95)
# but much slower when number of assets is larger (114-410)

nlambda.QP <- matrix(0, ncol = l.sp, nrow = 10)
for(i in seq(l.sp)){
  nlambda.QP[,i] <- sapply(1:10, function(x) length(analysis.sp$lambda.QP[[i]][[x]]))
}
colnames(nlambda.QP) <- label.sp 
nlambda.QP
#eg. ind.sp[[1]][[1]]
nw <- analysis.sp$nweights[1,1]
lam.err <- setdiff(((1:nw)/nw)^2, lambda.sp[[1]][[1]])
sp <- GetIndex(ind.sp[[1]][[1]], SP500)
QP.solve(sp$mu, sp$covar, sp$lB, sp$uB, lam.err[1]) 
# error 2: constraints are inconsistent, no solution

