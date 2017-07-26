MS_plot <- function(ms, col = adjustcolor("blue", alpha.f = 0.5), 
                    type = "o", pch =16, ...){  #list of weights_set, legend...
  plot(ms[,"Sig"], ms[,"Mu"], type = type, pch = pch, col = col,
       main = "Efficient Frontier", ylab = expression(mu(w)),
       xlab = expression(sigma(w)), ...)
}

MS_points <- function(ms, col =  adjustcolor("red", alpha.f = 0.5), pch = 16, type = "o", ...){  #list of weights_set, legend...
  points(ms[,"Sig"], ms[,"Mu"], type = type, pch = pch, 
         col =col, ...)
}


plot_EF <- function(result, mu, covar, EF = c("unpurge", "Excess", "Chull"), 
                    Smoo = FALSE, Ln = TRUE, n.insert = 20, point.label = FALSE){ 
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
    if(Ln){type[i] <- "o"}
    MS_points(ms[ind.pur, ], type = type[i], pch = pch[i], col = color[i], lty = lty[i])
    if(Smoo){
      ms.smooth <- insertW(result = result, mu = mu, covar = covar, EF = pur)
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

Total.plot <- function(ms.list){
  MS_plot(ms.list[[1]])
  if(length(ms.list) >1){
    for (i in 2:length(ms.list)){
      points(ms.list[[i]][,"Sig2"], ms.list[[i]][,"Mu"] , pch = 16, col = adjustcolor(col = i, alpha.f = 0.5))
    }
  }
  legend("bottomright", 
         legend = names(ms.list),
         lty = 1, col = c("blue", 2:length(ms.list)), pch = 16)
  
} 

IndexPlot <- function(analysis, name, ind.label, n, 
                      show_plot = c("nwbox", "nwplot", "kappa","rcond")){
  show_all <- c("nwbox", "nwplot", "kappa", "rcond")
  l <- length(exp_ind(1.2, 50, n)$ind)
  ind.label <- exp_ind(1.2, 50, n)$ind
  is.plot <- rep(FALSE, 6)
  is.i <- sapply(seq_along(show_plot), function(i) grep(show_plot[i], show_all, fixed = TRUE))
  is.plot[is.i] <- TRUE
  # 1. nweights_ind_boxplot
  if(is.plot[1]){ 
    boxplot(analysis$nweights, xlab = "Number of Assets", ylab = "Number", 
            main = paste(name, ": Number of Weights Sets", sep = ""), xaxt = "n") 
    axis(1, at = seq(l), labels = ind.label)
  }
  
  # 2. nweights_ind_plot
  if(is.plot[2]){ 
    plot(apply(analysis$nweights, 2, mean), 
         ylim = range(c(analysis$nweights, analysis$nweights_unpurge)), 
         pch = 16, col = "red",
         xlab = "Number of Assets", ylab = "Number", type = "o",
         main = paste(name, ": Number of Weights Sets", sep = ""), xaxt = "n")
    axis(1, at = seq(l), labels = ind.label)
    lapply(1:10, function(x) points(analysis$nweights[x,], pch = 16, 
                                    col = adjustcolor("blue", 0.5)))
    lapply(1:10, function(x) points(analysis$nweights_unpurge[x,], pch = 16, 
                                    col = adjustcolor("yellow", 0.5)))
    legend("topleft", legend = c("mean-purged", "purged", "unpurged"), 
           col = c("red", "blue", "yellow"), lwd = 1, pch =16)
  } 
  
  # 3. kappa
  if(is.plot[3]){ 
    boxplot(analysis$covF.kappa, xlab = "Number of Assets", ylab = "kappa", 
            main = paste(name, ": Max-kappa of CovF", sep= ""), xaxt = "n")
    axis(1, at = seq(l), labels = ind.label)
  }
  
  # 4. rcond
  if(is.plot[4]){ 
    boxplot(analysis$covF.rcond, xlab = "Number of Assets", ylab = "rcond", 
            main = paste(name, ": Min-rcond of CovF", sep =""), xaxt = "n")
    axis(1, at = seq(l), labels = ind.label)
  }
}

plot_method_tcompare <- function(analysis, name, cla.mt, scale.log = FALSE){
    # compare t qp, cccp from analysis and t of ver6
    l1 <- length(analysis$t_compare) #length of qp, cccp
    l <- nrow(cla.mt)
    t_mean <- t(sapply(seq(l1), function(x) apply(analysis$t_compare[[x]], 1, mean)))/1000
  
    k <- ncol(t_mean)
    color <- c("red", "blue", "green")[seq(k)]
    color.adjust <- sapply(1:k, function(k) adjustcolor(color[k], 0.3))
    
    is.log <- ifelse(scale.log, "y", "")
    plot(cla.mt[,1], ylim = range(c(cla.mt[,1], t_mean)), type = "o",
         xaxt = "n", pch = 16, col = color[1],
         xlab = "Number of Assets", ylab = "seconds", 
         main = paste(name, ": Time Comparison", sep = ""), log = is.log)
    lapply(2:k, function(k) points(t_mean[, k], col = color[k], pch = 16, type = "o"))
    axis(1, at = seq(l), labels = row.names(cla.mt))
    legend("topleft", legend = c("cla", "qp", "cccp")[1:k], 
           col = color[1:k], pch = 16, lwd = 1)
}

plot_cla_tcompare <- function(cla.mt, cla.mt.pre, scale.log = FALSE){
  l <- nrow(cla.mt)
  is.log <- ifelse(scale.log, "y", "")
  label.log <- ifelse(scale.log, "Log", "")
  plot(seq(l), cla.mt.pre[, 1], type = "o", col = adjustcolor("blue", 0.8), 
       pch = 16, xaxt = "n", main = paste("CLA time", label.log, "Comparison"), 
       xlab = "number of assets", ylab = "time(s)", 
       ylim = range(cbind(cla.mt[, 1], cla.mt.pre[, 1])), log = is.log)
  lines(seq(l), cla.mt[, 1], type = "o", col = adjustcolor("red", 0.8), pch = 16)
  axis(1, at = seq(l), labels = row.names(cla.mt))
  legend("topleft", legend = c("Norring's version", "improved version"), pch = 16,
         col = c(adjustcolor("blue", 0.8), adjustcolor("red", 0.8)), lwd = 1)
}

time_lm <- function(analysis, n, vector.predict, ind.label){ # vector to fit
  x <- seq_along(ind.label)
  t.mean.log <- log(sapply(x, 
                           function(x) apply(analysis$t_compare[[x]], 1, mean)))
  k <- nrow(t.mean.log)
  fit <- lapply(1:k, function(i) lm(t.mean.log[i, ] ~ x ))
  names(fit) <- c("cla", "qp", "cccp")[1:k]
  new <- data.frame(x = vector.predict)
  p <- sapply(1:k, function(x) predict.lm(fit[[x]], new, se.fit = TRUE)$fit)
  if(is.matrix(p)){ 
    colnames(p) <- c("cla", "qp", "cccp")[1:k]
    rownames(p) <- vector.predict
  }else{names(p) <- c("cla", "qp", "cccp")[1:k]}
  list(fit = fit, predict = p)
}

t_lm <- function(mt.list){ # vector to fit
  l <- length(mt.list)
  x <- exp_ind(1.2, 50, 2196)$power
  fit <- lapply(seq(l), function(i) lm(mt.list[[i]] ~ x[seq_along(mt.list[[i]])] ))
  names(fit) <- names(mt.list)
  #new <- data.frame(x = n.predict)
  #p <- sapply(seq(l), function(x) predict.lm(fit[[x]], new, se.fit = TRUE)$fit)
  #list(fit = fit, predict = p)
  fit
}

plot_EF_bound <- function(l, u, ind = 1:50) {
  n <- length(l)
  ms <- list()
  for(i in seq(n)){
    assets.lb <- GetIndex(ind, assets1457, lB = rep(l[i], length(ind)), 
                          uB = rep(u[i], length(ind)))
    ms[[i]] <- MS.insert(Env7$cla.solve(assets.lb, purge = F)$weights_set,
                         assets.lb$mu, assets.lb$covar)
  }
  y.lim <- range(sapply(seq(n), function(i) range(ms[[i]][, "Mu"])))
  label <- rep(0, n)
  MS_plot(ms[[1]], type = "l", col = 1, ylim = y.lim)
  for(i in 2:n){
    points(ms[[i]], type = "l", col = i)
    label[i] <- paste0("(", l[i], ",", u[i], ")")
  }
  
  legend("bottomright", legend = label, col = seq(n), lty = 1,
         title = "weight bounds")
}






