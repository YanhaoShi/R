Index.ms2 <- function(Index){
  nr <- nrow(Index) 
  nc <- ncol(Index)
  lr <- log(Index[-1, ]/Index[-nr,]) #log return
  loss <- -lr
  gfit <- lapply(loss, garchFit,  formula = ~ garch(1,1),
         cond.dist = "std", trace = FALSE)
  gpredict <- t(sapply(gfit, function(x) predict(x, n.ahead = 1)))
  mu <- -unlist(gpredict[, "meanForecast"])
  std <- unlist(gpredict[, "standardDeviation"])
  index.stand <- sapply(seq(nc), function(x) lr[,x]/gfit[[x]]@sigma.t)
  diag.std <- diag(std)
  covar <- diag.std %*% cor(index.stand) %*% diag.std
  list(mu = mu, covar = covar)
}

