# 1) fiiting.plot
# a function returns graphs of the kth stock in index, 
# comparing methods GHD, HYP, NIG and Normal Distribution fitting the density distribution 
fitting.plot <- function(index, k){ 
  y <- timeSeries(index[,k], charvec = row.names(index))
  name <- colnames(index)[k]
  y1 <- na.omit(diff(log(y))*100) #take logarithm of raw time series and differencing
  y1.den <- density(y1)
  #GHD
  GHD.fit <- fit.ghypuv(y1, symmetric = FALSE, silent = TRUE) # non-zero skewness 
  GHD.den <- dghyp(y1.den$x, GHD.fit)
  #HYP
  HYP.fit <- fit.hypuv(y1, symmetric = FALSE, silent = TRUE)
  HYP.den <- dghyp(y1.den$x, HYP.fit)    
  #NIG
  NIG.fit <- fit.NIGuv(y1, symmetric = FALSE, silent = TRUE)
  NIG.den <- dghyp(y1.den$x, NIG.fit)
  
  #Normal distribution:
  ND.den <- dnorm(y1.den$x, mean = mean(y1), sd = sd(y1))
  
  # graph1
  plot(y1.den, xlab="", ylab="", ylim = c(0, max(GHD.den, HYP.den, NIG.den, ND.den)),
       main = "1.Fitted density functions")
  lines(y1.den$x, GHD.den, col = "red")
  lines(y1.den$x, HYP.den, col = "green")
  lines(y1.den$x, NIG.den, col = "orange")
  lines(y1.den$x, ND.den, col = "blue")
  legend("topleft",
         legend = c("empirical", "GHD", "HYP", "NIG", "NORM"),
         col = c("black", "red", "green", "orange", "blue" ), 
         lty = 1)
  
  # graph2
  qqghyp(GHD.fit, line = TRUE, plot.legend = FALSE,
         gaussian = FALSE, main = "", cex = 0.8, ghyp.col = "red")
  title("2.qqplot")
  qqghyp(HYP.fit, add = TRUE, plot.legend = FALSE,
         gaussian = FALSE, main = "", cex = 0.8, ghyp.col = "green")
  qqghyp(NIG.fit, add = TRUE, plot.legend = FALSE,
         gaussian = FALSE, main = "", cex = 0.8, ghyp.col = "orange")
  legend("topleft", legend = c("GHD", "HYP", "NIG"), col = c("red", "green", "orange"), lty = 1)
  
  # graph3
  i <- seq(0.001, 0.05, 0.001)
  GHD.VaR <- qghyp(i, GHD.fit)
  HYP.VaR <- qghyp(i, HYP.fit)
  NIG.VaR <- qghyp(i, NIG.fit)
  ND.VaR <- qnorm(i, mean = mean(y1), sd = sd(y1)) 
  y.VaR <- quantile(x = y1, probs = i)
  plot(y.VaR, ylim = range(y.VaR, GHD.VaR, HYP.VaR, NIG.VaR, ND.VaR), 
       type = "l", main = "VaR", xaxt = "n")
  axis(1, at = 10*(0:10),labels = paste(0:10,"%"))
  lines(GHD.VaR, col = "red")
  lines(HYP.VaR, col = "green")
  lines(NIG.VaR, col = "orange")
  lines(ND.VaR, col = "blue")
  legend("bottomright",
         legend = c("empirical", "GHD", "HYP", "NIG", "NORM"),
         col = c("black", "red", "green", "orange", "blue" ), lty = 1)
  
  # graph4
  GHD.ES <- ESghyp(i, GHD.fit)
  HYP.ES <- ESghyp(i, HYP.fit)
  NIG.ES <- ESghyp(i, NIG.fit)
  ND.ES <- mean(y1) - sd(y1) *dnorm(qnorm(1 - i))/i
  y.ES <- sapply(ceiling(i * length(y1)), function(x) mean(sort(c(y1))[1:x]))#empirical quantile
  
  plot(y.ES, type = "l", xlab = "", ylab = "ES",
       ylim = range(c(GHD.ES, HYP.ES, NIG.ES, ND.ES, y.ES)),
       main = "4.Expected Shortfall", xaxt = "n")
  axis(1, at = 10*(0:10),labels = paste(0:10,"%"))
  lines(GHD.ES, col = "red")
  lines(HYP.ES, col = "green")
  lines(NIG.ES, col = "orange")
  lines(ND.ES, col = "blue")
  legend("bottomright",
         legend = c("empirical", "GHD", "HYP", "NIG", "Normal"),
         col = c("black", "red", "green", "orange", "blue"), lty = 1)
  
  # determine which model (among GHD, HYP and NIG) is the best in AIC selection
  AIC <- stepAIC.ghyp(y1, dist = c("ghyp", "hyp", "NIG"),
                      symmetric = FALSE,
                      control = list(maxit = 1000), silent = TRUE)
  list(Stock = paste("index name:", name), AIC = AIC$fit.table)
}
########################################################################
index.df # brief infomation of FRAPO index

#for example
fitting.plot(SP500, 21) #weekly observations
fitting.plot(ESCBFX, 3) #daily observations
fitting.plot(NASDAQ, 10) #weekly observations


fitting.plot(StockIndexAdj,1)# weekly: HYP best; (Company dataset: weekly)
fitting.plot(StockIndexAdjD,1)# daily: GHD best
