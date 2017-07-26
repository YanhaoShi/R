## Weights
MS <- function(weights_set, mu, covar){ 
    Sig2 <- colSums(weights_set *(covar %*% weights_set) )
    cbind(Sig = sqrt(Sig2), Mu = as.vector(t(weights_set) %*% mu))
}

# add between each two weights columns.
MS.insert <- function(weights_set, mu, covar, n.insert = 20){ 
  ms <- MS(weights_set, mu, covar)
  #ind <- testPurgeChull(result, mu, covar)[[paste0("ind.", tolower(EF))]]
 # l <- length(ind)
 # weights <- weights_set[, ind]
 # ms <- ms[ind, ]
  
  d <- sapply(seq(nrow(ms)-1), function(j) dist(ms[j:(j+1),])) # distance of each to points
  nd <- ceiling(d/max(d) * n.insert) # insert 20 points in longest distance
  w.new <- weights_set[, 1]
  for(k in seq_along(nd)){
    d.seq <- seq(0, 1, length = nd[k] + 2)[-1]
    w.seq <- weights_set[, k]%*% t(1 - d.seq) + weights_set[, k + 1] %*% t(d.seq)
    w.new <- cbind(w.new, w.seq, deparse.level = 0L)
  }
  MS(w.new, mu, covar) # return inserted weights_set
}

# piecewise hyperbolic on Efficient Frontier, return cov of mu_in
hyperbolic <- function(mu_in, w1, w2, Mu, Cov){  
  mu1 <- sum(Mu * w1)   # mu1 > mu2
  mu2 <- sum(Mu * w2)
  lambda <- (mu_in - mu1)/(mu2 - mu1)
  w_in <- (1-lambda)*w1 + lambda*w2
  sum(w_in * (Cov %*% w_in))  # sig2 of mu_in
}
## free indices
f_change <- function(free_indices){ # free index change of result$free_indices
  n <- length(free_indices)
  f.label <- c()
  for(i in 2: n){
    if(!anyNA(free_indices[[i]])){
      s <- sign(length(free_indices[[i]]) - length(free_indices[[i-1]]))
      if(s == 0){
        f.label[i] <- NA
      } else if(s > 0){
        f.label[i] <- setdiff(free_indices[[i]], free_indices[[i-1]])*s
        #f.label[i] <- paste0("(", paste(diff, collapse=","), ")") 
      }else if(s < 0){
        f.label[i] <- setdiff(free_indices[[i-1]], free_indices[[i]])*s
        #f.label[i] <- paste0("(", paste(diff, collapse=","), ")")
      }
    }else f.label[i] <- NA
  }
  f.label
} 

Mul_mucov <- function(assets, t.mu, t.cov){
  r1 <- Env6$cla.solve(assets$mu, assets$covar, assets$lB, assets$uB)
  r2 <- Env6$cla.solve(assets$mu * t.mu, assets$covar, assets$lB, assets$uB)
  r3 <- Env6$cla.solve(assets$mu, assets$covar * t.cov, assets$lB, assets$uB)
  r4 <- Env6$cla.solve(assets$mu * t.mu, assets$covar * t.cov, assets$lB, assets$uB)
  r <- list(r1, r2, r3, r4)
  m.name <- c("mu, cov", paste0("mu*", t.mu, ", cov"),
              paste0("mu, cov*", t.cov), paste0("mu*", t.mu, ", cov*", t.cov))
  n <- ncol(r1$weights_set)
  lam <- sapply(1:4, function(i) r[[i]]$lambdas[-c(1,n)]) #remove NA and 0
  gam <- sapply(1:4, function(i) r[[i]]$gammas[-1]) #remove NA
  ratio.lam <- apply(lam/lam[, 1], 2, mean)
  ratio.gam <- apply(gam/gam[, 1], 2, mean)
  names(ratio.lam) <- m.name
  names(ratio.gam) <- m.name
  is.free_indices <- all(sapply(1:3, function(i) all.equal(r[[i]]$free_indices, 
                                                           r[[i+1]]$free_indices))) 
  is.weights_set <- all(sapply(1:3, function(i) all.equal(r[[i]]$weights_set, 
                                                          r[[i+1]]$weights_set))) 
  list(ratio.lam = ratio.lam, ratio.gam = ratio.gam, 
       is.free_indices = is.free_indices, is.weights_set = is.weights_set)
}

testLam <- function(lam, equal.tol){
  nlam <- length(lam)
  lam.equal <- sapply(seq(nlam -1), 
                      function(m) isTRUE(all.equal(lam[m], lam[m+1], tol = equal.tol)))
  list(lam.equal = lam.equal, n.equal = sum(lam.equal))
}

plotNLam <- function(lam, is.label = FALSE){
  
  n.equal <- sapply(1:23, function(x) testLam(lam, equal.tol = 10^-(x))$n.equal )
  plot(0:23, c(length(lam), n.equal), type = "o", col = "blue", pch = 16,
       ylim = c(min(n.equal), max(n.equal + 50)))
  if(is.label) text(0:23, c(length(lam), n.equal), c(length(lam), n.equal), pos = 3)
}

mLam.tol <- function(lam.list){
  l <- length(lam.list)
  lam.name <- names(lam.list)
  n.equal <- matrix(0, nrow = 23, ncol = l)
  nlam <- sapply(lam.list, length)
  for(j in seq(l)){
    n.equal[, j] <- sapply(1:23, function(x) testLam(lam.list[[j]], 
                                                     equal.tol = 10^-(x))$n.equal)
  }
  n.equal <- rbind(nlam, n.equal)
  colnames(n.equal) <- lam.name
  rownames(n.equal) <- c("nlambda", paste0("equal.tol-", 1:23))
  n.equal
}

plotLam.tol <- function(lam.mat, start.tol = -3){
  l <- ncol(lam.mat)
  lam.tol <- (-18:0)[-((l+start.tol+1):(l-1))]
  lam.mat <- lam.mat[ ,-((l+start.tol+1):(l-1))]
  l <- ncol(lam.mat)
  plot(0:23, lam.mat[, l], type = "n", col = l, pch = 16, xaxt = "n",
       ylim = c(min(lam.mat) - 50, max(lam.mat + 50)), 
       ylab = "# duplicated lambda pairs",
       xlab = "tolerance of all.equal()",
       main = "Number of duplicated lambda pairs on tolerance scale")
  axis(1, at = 0:23, labels = c("nlambda", paste0("1e-", 1:23)))
  col <- sapply(c("black", "yellow", "blue", "red", "green", "orange", "purple"), adjustcolor, 0.8)
  ind.unique <- which(!duplicated(t(lam.mat))) # unique columns of lambda numbers
  lam.unique <- lam.mat[, ind.unique]
  n.unique <- length(ind.unique)
  label <- rep(0, n.unique)
  for(i in seq(n.unique)){
    lam.col <- lam.mat[, ind.unique[i]]
    points(0:23, lam.col, type = "o", col = col[i], pch = 16)
    label[i] <- paste0("tol 1e", paste(lam.tol[which(colSums(lam.mat == lam.col) == 24)], 
                                       collapse = ","))
  }
  min.ind <- which.min(lam.unique[nrow(lam.unique), ])
  max.ind <- which.max(lam.unique[nrow(lam.unique), ])
  text(0:23, lam.unique[, min.ind], lam.unique[, min.ind], col = "black",  pos = 1)
  text(0:23, lam.unique[, max.ind], lam.unique[, max.ind], col = "black", pos = 3)
  
  legend("topright", legend = label, col = col[seq(n.unique)], 
         bg = "transparent", bty = "n", pch = 16, lty = 1)
  abline(v = 16, col = "darkgrey", lty = 2)
}

testRep <- function(r){
  d <- setdiff(seq_along(r$lambdas), purgeChull(r$MS_weight))
  nd <- length(d)
  rep <- rep(0, nd)
  for(i in seq(nd)){
    rep[i] <- any(c(isTRUE(all.equal(r$MS_weight[d[i], ], r$MS_weight[d[i]-1, ])),
                    isTRUE(all.equal(r$MS_weight[d[i], ], r$MS_weight[d[i]+1, ]))))
  }
  all(rep == 1)
}








