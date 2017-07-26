plot_fit <- function(fit.list, analysis, name, n, n.predict, ind.label){
  pdf(d.file(paste(name, "_lmplot.pdf", sep = ""), exists = FALSE)) 
  l <- length(ind.label)
  label <- exp_ind(1.2, 50, n)$ind
  t.mean.log <- log(sapply(seq(l), 
                           function(x) apply(analysis$t_compare[[x]], 1, mean)))
  k <- nrow(t.mean.log)
  ## 1.
  plot(c(1,l), range(t.mean.log), xaxt = "n", col = "green", pch = 16, 
       type = "n", xlab = "number of assets", ylab = "log(milliseconds)",
       main = paste(name, ": linear models of computation time"))
  axis(1, at = seq(l), labels = label)
  
  color <- c("red", "blue", "green")[1:k]
  lapply(1:k, function(x) {points(seq(l), t.mean.log[x,], col = color[x], pch = 16)
    abline(fit.list[[x]], col = color[x])})
  legend("topleft", legend = c("cla", "qp", "cccp")[1:k], 
         col = color[1:k],  pch = 16, lwd = 1)
  
  l.predict <- length(exp_ind(1.2, 50, n.predict)$ind)
  label.predict <- exp_ind(1.2, 50, n.predict)$ind
  
  pf <- function(k){
    fit <- fit.list[[k]]
    new <- data.frame(x = seq.int(l.predict))
    predict.lm(fit, new, se.fit = TRUE) 
    pred.pi <- predict.lm(fit, new, interval = "prediction")
    pred.ci <- predict.lm(fit, new, interval = "confidence")
    matplot(new$x, exp(cbind(pred.ci, pred.pi[,-1])),
            lty = 1, type = "l", ylab = "predicted y",
            xaxt = "n", xlab ="number of assets", 
            col = c("black", "blue", "blue", "red", "red"),
            main = paste(name, ":regression analysis of", names(fit.list)[k]))
    
    axis(1, at = seq(label.predict), labels = label.predict)
    legend("topleft", legend = c("confidence interval", "prediction interval"),
           col = c("blue", "red"), lty = 1)
  }
  lapply(1:k, pf)
  ##
  new <- data.frame(x = seq.int(l.predict))
  pm <- sapply(1:k, function(x) predict.lm(fit.list[[x]], new, se.fit = TRUE)$fit)
  color <- c("blue", "red", "green")
  plot(pm[,1], type = "l", ylim = range(pm), col = color[1],
       xlab = "number of assets", ylab = "log(milliseconds)" ,
       main = "prediction", xaxt = "n")
  axis(1, at = seq(label.predict), labels = label.predict)
  lapply(2:k, function(x) lines(pm[,x], 
                                type = "l", col = color[x]))
  legend("topleft", legend = c("cla", "qp", "cccp"), col = color, lwd = 1)
  dev.off()
  path <- gsub("/", "\\\\", d.file(paste(name, "_lmplot.pdf", sep = ""), exists = FALSE), 
               perl=TRUE)
  system(paste0('open "', path, '"'))
}

