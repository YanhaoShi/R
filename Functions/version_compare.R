#version.compare <- function(assets.number, times, versions){
#  assets <- GetAssets(1:assets.number, assets1457)
#  micro <- microbenchmark(r1 <- CLA[[paste("M", versions[1], sep = "")]](assets$mu, assets$covar, assets$lB, assets$uB),
#                          r2 <- CLA[[paste("M", versions[2], sep = "")]](assets$mu, assets$covar, assets$lB, assets$uB),
#                          times = times)
#  list(micro = micro, is.equal = all.equal(r1, r2),
#       weight.equal = all.equal(r1$weights_set_purge, r2$weights_set_purge),
#       free_indices.equal = all.equal(r1$free_indices, r2$free_indices))
# }

# (1) purgeNumErr
testPurgeNumErr <- function(result, lB, uB){
  free_indices <- result$free_indices
  gammas <- result$gammas
  lambdas <- result$lambdas
  weights_set <- result$weights_set
  ind.old <- Env1.2$purgeNumErr(free_indices, gammas, lambdas, 
                                lB, uB, weights_set, 1e-09)$index
  ind.old <- seq_along(ind.old) - 1 + ind.old
  
  ind <- which(Env2$purgeNumErr(lB = lB, uB = uB, weights_set = weights_set, 
                                tol.s = 1e-09, tol.b = 1e-10))
  ind.new <- seq_len(ncol(weights_set))[-ind]
  list(ind.old = ind.old, ind.new = ind.new) # return indices should remove
}

plot_weight <- function(result, lB, uB, plot.show){ # plot.show = "f" or "sum"
  free_indices <- result$free_indices
  weights_set <- result$weights_set
  f.label <- f_change(free_indices)
  n <- length(free_indices)
  y <- log(1e-10 + abs(colSums(weights_set) - 1)) # deviation from 1
  
  if(plot.show == "f"){
    plot(seq(n), y, ylim = c(min(y), max(y) +10), type = "o", pch =16,
         col = c(1, sign(f.label[-1])+3 ),
         xlab = "steps of CLA", main = "PurgeNumErr: log deviation from 1",
         ylab = paste0("deviation: ", expression(10^y)))
    text(2:(n-1), log(1e-10 + abs(colSums(weights_set) - 1))[-1], 
         f.label[-1], pos = 3, col = sign(f.label[-1])+3  )
    legend("topleft", legend = c("add f", "remove f"), col = c(4,2), pch = 16,
           bg = "transparent", bty = "n")
  }
  if(plot.show == "sum"){
    ind <- testPurgeNumErr(result, lB, uB)
    ind.old <- ind$ind.old
    ind.new <- ind$ind.new
    col2 <- rep("grey", n)
    col2[ind.old] <- 
      col2[ind.new] <- adjustcolor("purple", 0.8)
    plot(y, ylim = c(min(y), max(y) +10), type = "o", pch =16,
         col = "grey",
         xlab = "steps of CLA", main = "PurgeNumErr: Steps to be purged",
         ylab = paste0("deviation: ", expression(10^y)))
    points(seq_len(n)[ind.old], y[ind.old], col = "blue3", pch = 15, cex = 1.5)
    points(seq_len(n)[ind.new], y[ind.new], col = "red", pch = 18)
    text(ind.new, y[ind.new], 
         ind.new, pos = 3, col = "red"  )
    legend("topleft", legend = c("old PurgeNumErr", "new PurgeNumErr"), 
           col = c("blue3", "red"), 
           pch = c(15, 18), bg = "transparent", bty = "n")
  }
}

######################################################################################
# (2) purgeChull
testPurgeChull <- function(result, mu, covar){
  free_indices <- result$free_indices
  gammas <- result$gammas
  lambdas <- result$lambdas
  weights_set <- result$weights_set
  ind <- Env1.2$purgeExcess(mu = mu, lambdas = lambdas, gammas = gammas, 
                                free_indices = free_indices, 
                                weights_set = weights_set)$k
  ind.excess <- seq_len(ncol(weights_set))[-ind]
  ind.chull <- Env6$purgeChull(weights_set , mu, 
                  covar = covar)
  ind.slope <- Env4$purgeSlope(weights_set = weights_set ,mu = mu, 
                               covar = covar)
  list(ind.excess = ind.excess, ind.chull = ind.chull,
       ind.slope = ind.slope, ind.unpurge = seq_along(gammas)) #return keep
}

plot_EF_purge <- function(result, mu, covar, EF = c("unpurge", "Excess", "Chull"), 
                    Smoo = FALSE, n.insert = 20, point.label = FALSE){ 
  # smooth = T or F, line = T or F, EF = "unpurge" or "Excess" or "Chull"
  weights_set <- result$weights_set
  ms <- MS(weights_set, mu, covar)
  ind <- testPurgeChull(result, mu, covar)
  MS_plot(ms, type = "n")
  EF.label <- c("unpurge", "Excess", "Chull")
  color <- c(adjustcolor("black", 0.8), adjustcolor("blue", 0.5), adjustcolor("red", 0.8))
  pch <- c(1, 16, 17)
  type <- c("p", "p", "p")
  lty <- c(0, 2, 1)
  j <- c()
  for(pur in EF){
    i <- grep(pur, EF.label)
    j <- c(j, i)
    ind.pur <- ind[[paste0("ind.", tolower(pur))]]
    if(!Smoo){type[i] <- "o"}
    MS_points(ms[ind.pur, ], type = type[i], pch = pch[i], col = color[i], lty = lty[i])
    if(Smoo){
      ms.smooth <- MS.insert(result$weights_set[, ind.pur], mu = mu, covar = covar)
      MS_points(ms.smooth, type = "l", pch = pch[i], col = color[i], lty = lty[i])
    }
  }
  if(point.label){
    n <- nrow(ms)
    pos <- rep(1, n)
    pos[ind$ind.chull] <- 3
    text(ms[,1], jitter(ms[,2]), seq(n), pos = pos ) 
  }
  
  legend("bottomright", legend = EF.label[j],
         col = color[j], pch = pch[j], 
         bg = "transparent", bty = "n", lty = lty[j])
}


purgeChull <- function(ms){
  ch <- chull(ms)
  sort(ch)
}

purgeUnique <- function(ms){
  which(!duplicated(ms))
}



testUnique <- function(cla.input){
  result <- Env6.1$cla.solve(cla.input, purge = FALSE)
  w <- result$weights_set
  ms <- result$MS_weight
  j1 <- purgeChull(round(ms, 8)) ## round at digits == 8
  j2 <- purgeUnique(round(ms, 8))
  t.chull <- system.time(replicate(100, purgeChull(ms)))[1] ## replicate for 100 times
  t.unique <- system.time(replicate(100, purgeUnique(ms)))[1]
  is.equal <- ifelse(isTRUE(all.equal(w[, j1], w[, j2])), "T", "F")
  list(time = c(chull = t.chull, unique = t.unique),
       diff.num = c(is.equal, length(j2)-length(j1))) # chull purge more than unique
}

plot_unipurge <- function(uni_purge, index.name = "Index"){
  mt <- uni_purge$time.mean
  lab <- row.names(mt)
  plot(seq_along(lab), mt[, 1], ylim = range(mt), type = "l", 
       col = "red", xaxt = "n", xlab = "number of assets",
       ylab = "100 times", 
       main = paste0(index.name,": time of purgeChull and purgeUnique"))
  axis(1, at = seq_along(lab), labels = lab)
  points(mt[, 2], col = "blue", type = "l")
  legend("topleft", legend = c("purgeChull", "purgeUnique"), 
         col = c("red", "blue"), lwd = 1)
}
#############################################################################
# (3) solve covarF
# time costs of inverse methods on a given covarF, ...
Inv_covarF <- function(get){
  covarF <- get$covarF
  covarFB <- get$covarFB
  muF <- get$muF
  wB <- get$wB
  
  inv1 <- replicate(5, system.time(chol2inv(chol(covarF)) %*% cbind(1, muF, covarFB %*% wB)))
  inv2 <- replicate(5, system.time(solve(covarF) %*% cbind(1, muF, covarFB %*% wB)))
  inv3 <- replicate(5, system.time(solve(covarF, cbind(1, muF, covarFB %*% wB))))
  inv <- rbind(apply(t(inv1), 2, mean)[1:3], apply(t(inv2), 2, mean)[1:3],
               apply(t(inv3), 2, mean)[1:3])
  inv
}

Inv_t <- function(ind, total_assets){
  assets <- GetIndex(ind, total_assets)
  #result <- readRDS(d.file())
  result <- Env6$cla.solve(assets$mu, assets$covar, assets$lB, assets$uB)
  free_indices <- result$free_indices
  weights_set <- result$weights_set
  n <- ncol(weights_set)
  
  get <- list()
  for(j in seq(n-1)){
    get[[j]] <- Env6$getMatrices(assets$mu, assets$covar, weights_set[, j], free_indices[[j]])
  }
  
  inv.list <- lapply(seq(n-1), function(i) Inv_covarF(get[[i]]))
  inv.t <- c()
  for(i in seq(n-1)){
    inv.t <- rbind(inv.t, inv.list[[i]])
  }
  inv1.t <- apply(inv.t[seq(n-1)*3-2, ], 2, sum)
  inv2.t <- apply(inv.t[seq(n-1)*3-1, ], 2, sum)
  inv3.t <- apply(inv.t[seq(n-1)*3, ], 2, sum)
  rbind(inv1.t, inv2.t, inv3.t)
}

testInv <- function(ind.list, total_assets){
  t.m <- c()
  for(i in seq_along(ind.list)){
    t.list <- lapply(1:10, function(j) Inv_t(ind.list[[i]][[j]], total_assets))
    
    for(j in seq(10)){
      t.m <- rbind(t.m, t.list[[j]])
    }
  }
  n <- length(ind.list) * 10
  ind.label <- exp_ind(1.2, 50, length(total_assets$mu))$ind
  inv.name <-  paste0(rep(ind.label, each = 10), "set", each = 1:10)
  inv1.matrix <- t.m[seq(n)*3 -2, ]
  inv2.matrix <- t.m[seq(n)*3 -1, ]
  inv3.matrix <- t.m[seq(n)*3   , ]
  row.names(inv1.matrix) <- inv.name
  row.names(inv2.matrix) <- inv.name
  row.names(inv3.matrix) <- inv.name
  
  inv.m1 <- sapply(seq_along(ind.label), function(i) mean(inv1.matrix[1:10 + (i - 1)*10, 1]))
  inv.m2 <- sapply(seq_along(ind.label), function(i) mean(inv2.matrix[1:10 + (i - 1)*10, 1]))
  inv.m3 <- sapply(seq_along(ind.label), function(i) mean(inv3.matrix[1:10 + (i - 1)*10, 1]))
  inv.mean <- cbind(inv.m1, inv.m2, inv.m3)
  row.names(inv.mean) <- ind.label
  list(inv1 = inv1.matrix, inv2 = inv2.matrix, inv3 = inv3.matrix, inv.mean = inv.mean)
}

plot_Inv <- function(covF_inv, index.name = "Index"){
  inv.mt <- covF_inv$inv.mean
  lab <- row.names(inv.mt)
  l <- seq_along(lab)
  plot(l, inv.mt[, 1], ylim = range(inv.mt), type = "l", 
       col = "black", xaxt = "n", 
       main = paste0(index.name, ": Time cost of covF inverse"),
       xlab = "number of assets", ylab = "system.time")
  lines(l, inv.mt[, 2], col = "blue")
  lines(l, inv.mt[, 3], col = "red")
  axis(1, at = l, labels = lab)
  legend("topleft", legend = c("method 1", "method 2", "method 3"),
         col = c("black", "blue", "red"), lty = 1)
}

#############################################################################
# (4) lambdas
plot_lambda <- function(result7.1, result, mu, covar, smoo = FALSE){
  ms1 <- result7.1$MS_weight
  ms2 <- result$MS_weight
  MS_plot(ms1, col = adjustcolor("blue", 0.5), type = "p")
  MS_points(ms2, col = adjustcolor("red", 0.5), type = "p")
  
  if(smoo){
    ms1 <- MS.insert(result7.1$weights_set, mu, covar)
    ms2 <- MS.insert(result$weights_set, mu, covar)
  }
  MS_points(ms1, col = adjustcolor("blue", 0.5), type = "l")
  MS_points(ms2, col = adjustcolor("red", 0.5), type = "l")
  
  legend("bottomright", legend = c("EF", "EF+"), lwd = 1, pch = 16,
         col = c(adjustcolor("red", 0.5), adjustcolor("blue", 0.5)))
}