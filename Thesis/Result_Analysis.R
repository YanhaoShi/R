if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
#############################################################################
## 1) Convexity of Efficient Frontier
l <- length(label.nasdaq)
rep.nasdaq <- matrix(0, nrow = 10, ncol = l)
for(i in seq(l)){
  rep.nasdaq[, i] <- sapply(1:10, function(j) 
    testRep(readRDS(d.file(paste0("NASDAQ", label.nasdaq[i], "_set", j, ".rds")))))
}
all(rep.nasdaq == 1) # TRUE

l <- length(label.sp)
rep.sp <- matrix(0, nrow = 10, ncol = l)
for(i in seq(l)){
  rep.sp[, i] <- sapply(1:10, function(j) 
    testRep(readRDS(d.file(paste0("SP", label.nasdaq[i], "_set", j, ".rds")))))
}
all(rep.sp == 1) # TRUE

l <- length(label.company)
rep.company <- matrix(0, nrow = 10, ncol = l)
for(i in seq(l)){
  rep.company[, i] <- sapply(1:10, function(j) 
    testRep(readRDS(d.file(paste0("Company", label.nasdaq[i], "_set", j, ".rds")))))
}
all(rep.company == 1) # TRUE

#########################################################################################
# 2 ) portfolio on efficient frontier
result <- result50.list$Env7
mu <- assets50$mu; covar <- assets50$covar

findSig <- function(Mu0, result, covar){ 
  ms.w <- result$MS_weight
  n <- nrow(ms.w)
  ## revert order for all (w, mu, sig):
  w <- result$weights_set[, n:1] 
  mu.w <- ms.w[n:1, "Mu"]
  sig.w <- ms.w[n:1, "Sig"]
  if(Mu0 < min(mu.w) | Mu0 > max(mu.w)) stop("not on EF")
  i <- findInterval(Mu0, mu.w)
  ## We have to catch the case mu[i] ~= mu[i+1]  (interpolation would divide by zero)
  if(isTRUE(all.equal(Mu0, mu.w[i]))) {
    list(Sig = sig.w[i], weight = w[, i])
  }else if(isTRUE(all.equal(Mu0, mu.w[i+1]))) {
    list(Sig = sig.w[i+1], weight = w[, i+1])
  }else{
    # solve for a : Mu0 = a* Mu1 + (1-a)* Mu2 :
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    w0 <- a* w[, i] + (1-a)*w[, i+1]
    list(Sig = colSums(w0 *(covar %*% w0) ), #return Sig0
         weight = w0)}
}


findSig0 <- function(Mu0, result, covar){ 
  ms.w <- result$MS_weight
  n <- nrow(ms.w)
  ## revert order for all (w, mu, sig):
  w <- result$weights_set[, n:1] 
  mu.w <- ms.w[n:1, "Mu"]
  sig.w <- ms.w[n:1, "Sig"]
  if(Mu0 < min(mu.w) | Mu0 > max(mu.w)) stop("not on EF")
  i <- findInterval(Mu0, mu.w)
  ## We have to catch the case mu[i] ~= mu[i+1]  (interpolation would divide by zero)
  if(isTRUE(all.equal(Mu0, mu.w[i]))) {
    sig.w[i]
  }else if(isTRUE(all.equal(Mu0, mu.w[i+1]))) {
    sig.w[i+1]
  }else{
    # solve for a : Mu0 = a* Mu1 + (1-a)* Mu2 :
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    w0 <- a* w[, i] + (1-a)*w[, i+1]
    colSums(w0 *(covar %*% w0))
    #return Sig0
  }    
}


findMu <- function(Sig0, result, covar, tol.unir = 1e-6){
  n <- nrow(ms.w)
  mu.w <- ms.w[n:1, "Mu"]
  sig.w <- ms.w[n:1, "Sig"]
  if(Sig0 < min(sig.w) | Sig0 > max(sig.w)) stop(sprintf("Sig0 must be in [%d,%d]",
                                                         min(sig.w), max(sig.w)))
  i <- findInterval(Sig0, sig.w)
  ## We have to catch the case mu[i] ~= mu[i+1]  (interpolation would divide by zero)
  if(isTRUE(all.equal(Mu0, mu.w[i]))) {
    list(mu = mu.w[i], weight = w[, i])
  } else if(isTRUE(all.equal(Mu0, mu.w[i+1]))) {
    list(mu = mu.w[i+1], weight = w[, i+1])
  } else {
    r <- uniroot(function(mu) findSig0(mu, result, covar) - Sig0, 
                 interval = sig.w[c(i,i+1)], tol=tol.unir)
    mu <- r$root
    # solve for a : mu = a* Mu1 + (1-a)* Mu2 :
    a <- (mu - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    w0 <- a* w[, i] + (1-a)*w[, i+1]
    list(mu = mu, weight = w0)
  }
}
findMu(0.02, result, covar)

findSig(Mu0 = 0.002, result, covar)


findMu <- function(Sig0, result, covar){ 
  n <- nrow(ms.w)
  w <- result$weights_set[, n:1]
  mu.w <- sort(ms.w[, "Mu"])
  sig.w <- sort(ms.w[, "Sig"])
  if(Sig0 < min(sig.w) | Sig0 > max(sig.w)) stop("not on EF")
  i <- findInterval(Sig0, sig.w)
  
  Cholesky(sig.w[i+1])
  
  
}

##############################################################################

r20 <- Env6$cla.solve(assets20$mu, assets20$covar, assets20$lB, assets20$uB)

Result_table(r20, assets20$mu, assets20$covar)
Result_table(r20, assets20$mu, assets20$covar)



Mul_mucov(assets = assets20, t.mu = 5, t.cov = 7)

#############################################################################
# Appendix B
company.tcompare <- cbind(company.cla.mt[,1], company.cla.mt.pre[,1])
colnames(company.tcompare) <- c("new version", "old version")
xtable(company.tcompare, digits = 3, caption = "Time Comparison of Company Assets")

sp.tcompare <- cbind(sp.cla.mt[,1], sp.cla.mt.pre[,1])
colnames(sp.tcompare) <- c("new version", "old version")
xtable(sp.tcompare, digits = 3, caption = "Time Comparison of SP500")

nasdaq.tcompare <- cbind(nasdaq.cla.mt[,1], nasdaq.cla.mt.pre[,1])
colnames(nasdaq.tcompare) <- c("new version", "old version")
xtable(nasdaq.tcompare, digits = 3, caption = "Time Comparison of NASDAQ")
