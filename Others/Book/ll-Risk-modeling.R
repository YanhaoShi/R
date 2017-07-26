## 3.2
library(zoo)
data(EuStockMarkets)
## Time Series plot of Levels
EuStockLevel <- as.zoo(EuStockMarkets)[, c("DAX", "CAC", "FTSE")]
plot(EuStockLevel, xlab = "", main = "")
## Perecntage returns
EuStockRet <- diff(log(EuStockLevel)) * 100
plot(EuStockRet, xlab = "", main = "")
## Cross correlations
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE))
ccf(EuStockRet[, 1], EuStockRet[, 2], ylab = "", xlab = "",
    lag.max = 20, main = "Returns DAX vs CAC")
ccf(abs(EuStockRet)[, 1], abs(EuStockRet)[, 2], ylab = "",
    xlab = "", lag.max = 20, main = "Absolute returns DAX vs CAC")
ccf(EuStockRet[, 1], EuStockRet[, 3], ylab = "", xlab = "",
    lag.max = 20, main = "Returns DAX vs FTSE")
ccf(abs(EuStockRet)[, 1], abs(EuStockRet)[, 3], ylab = "",
    xlab = "", lag.max = 20, main = "Absolute returns DAX vs FTSE")
ccf(EuStockRet[, 2], EuStockRet[, 3], ylab = "", xlab = "",
    lag.max = 20, main = "Returns CAC vs FTSE")
ccf(abs(EuStockRet)[, 2], abs(EuStockRet)[, 3], ylab = "",
    xlab = "", lag.max = 20, main = "Absolute returns CAC vs FTSE")
## Rolling correlations
rollc <- function(x){
  dim <- ncol(x)
  rcor <- cor(x)[lower.tri(diag(dim), diag = FALSE)]
  return(rcor)
}
rcor <- rollapply(EuStockRet, width = 250, rollc,
                  align = "right", by.column = FALSE)
colnames(rcor) <- c("DAX & CAC", "DAX & FTSE", "CAC & FTSE")
plot(rcor, main = "", xlab = "")

tail(rcor)
cor(EuStockRet)
cor(EuStockRet[,c("DAX", "CAC")])
colnames(EuStockRet)
rollapply(SP500[,1:4], width = 250, rollc,
          align = "right", by.column = FALSE)


## Ch 6.1
library(ghyp)
library(timeSeries)
library(fBasics)
## Return calculation
data(DowJones30)
y <- timeSeries(DowJones30[, "HWP"], charvec =
                  as.character(DowJones30[, 1]))  ##take the first column
yret <- na.omit(diff(log(y)) * 100)
## Fitting
ef <- density(yret)
ghdfit <- fit.ghypuv(yret, symmetric = FALSE,
                     control = list(maxit = 1000))
hypfit <- fit.hypuv(yret, symmetric = FALSE,
                    control = list(maxit = 1000))
nigfit <- fit.NIGuv(yret, symmetric = FALSE,
                    control = list(maxit = 1000))
## Densities
ghddens <- dghyp(ef$x, ghdfit) #density of GHD
hypdens <- dghyp(ef$x, hypfit)
nigdens <- dghyp(ef$x, nigfit)
nordens <- dnorm(ef$x, mean = mean(yret), sd = sd(c(yret[, 1])))
col.def <- c("black", "red", "blue", "green", "orange")
plot(ef, xlab = "", ylab = expression(f(x)), ylim = c(0, 0.25))
lines(ef$x, ghddens, col = "red")
lines(ef$x, hypdens, col = "blue")
lines(ef$x, nigdens, col = "green")
lines(ef$x, nordens, col = "orange")
legend("topleft",
       legend = c("empirical", "GHD", "HYP", "NIG", "NORM"),
       col = col.def, lty = 1)
## QQ-Plots
qqghyp(ghdfit, line = TRUE, ghyp.col = "red", plot.legend = FALSE,
       gaussian = FALSE, main = "", cex = 0.8)
qqghyp(hypfit, add = TRUE, ghyp.pch = 2, ghyp.col = "blue",
       gaussian = FALSE, line = FALSE, cex = 0.8)
qqghyp(nigfit, add = TRUE, ghyp.pch = 3, ghyp.col = "green",
       gaussian = FALSE, line = FALSE, cex = 0.8)
legend("topleft", legend = c("GHD", "HYP", "NIG"),
       col = col.def[-c(1,5)], pch = 1:3)
## Diagnostics
AIC <- stepAIC.ghyp(yret, dist = c("ghyp", "hyp", "NIG"),
                    symmetric = FALSE,
                    control = list(maxit = 1000))
LRghdnig <- lik.ratio.test(ghdfit, nigfit)
LRghdhyp <- lik.ratio.test(ghdfit, hypfit)



## Ch 6.2
## Probabilities
p <- seq(0.001, 0.05, 0.001)
## VaR
ghd.VaR <- abs(qghyp(p, ghdfit))
hyp.VaR <- abs(qghyp(p, hypfit))
nig.VaR <- abs(qghyp(p, nigfit))
nor.VaR <- abs(qnorm(p, mean = mean(yret), sd = sd(c(yret[, 1])))) 
emp.VaR <- abs(quantile(x = yret, probs = p)) 
# Plot of VaR
plot(emp.VaR, type = "l", xlab = "", ylab = "VaR", axes = FALSE,
     ylim = range(c(hyp.VaR, nig.VaR, ghd.VaR, nor.VaR, emp.VaR)))
box()
axis(1, at = seq(along = p), labels = names(emp.VaR), tick = FALSE) # x:level 99.9% - 95%
axis(2, at = pretty(range(emp.VaR, ghd.VaR, hyp.VaR,
                          nig.VaR, nor.VaR)))
lines(seq(along = p), ghd.VaR, col = "red")
lines(seq(along = p), hyp.VaR, col = "blue")
lines(seq(along = p), nig.VaR, col = "green")
lines(seq(along = p), nor.VaR, col = "orange")
legend("topright",
       legend = c("Empirical", "GHD", "HYP", "NIG", "Normal"),
       col = col.def, lty = 1)
## ES
ghd.ES <- abs(ESghyp(p, ghdfit))
hyp.ES <- abs(ESghyp(p, hypfit))
nig.ES <- abs(ESghyp(p, nigfit))
nor.ES <- abs(mean(yret) - sd(c(yret[, 1])) *
                dnorm(qnorm(1 - p)) / p)
obs.p <- ceiling(p * length(yret))
emp.ES <- sapply(obs.p, function(x) abs(mean(sort(c(yret))[1:x])))
## Plot of ES
plot(emp.ES, type = "l", xlab = "", ylab = "ES", axes = FALSE,
     ylim = range(c(hyp.ES, nig.ES, ghd.ES, nor.ES, emp.ES)))
box()
axis(1, at = 1:length(p), labels = names(emp.VaR), tick = FALSE)
axis(2, at = pretty(range(emp.ES, ghd.ES, hyp.ES, nig.ES, nor.ES)))
lines(1:length(p), ghd.ES, col = "red")
lines(1:length(p), hyp.ES, col = "blue")
lines(1:length(p), nig.ES, col = "green")
lines(1:length(p), nor.ES, col = "orange")
legend("topright",
       legend = c("Empirical", "GHD", "HYP", "NIG", "Normal"),
       col = col.def, lty = 1)



## Ch 6.3 triangle of points with different frequencies
rd <- c(1, 5, 10, 20, 40)
yrets <- na.omit(matrix(unlist(lapply(rd,
                                      function(x) diff(log(y), lag = x))), ncol = 5))
## Function for xi/chi coefficients
xichi <- function(x){
  param <- coef(x, type = "alpha.delta")
  rho <- param[["beta"]] / param[["alpha"]]
  zeta <- param[["delta"]] * sqrt(param[["alpha"]]^2 -
                                    param[["beta"]]^2)
  xi <- 1 / sqrt(1 + zeta)
  chi <- xi * rho
  result <- c(chi, xi)
  names(result) <- c("chi", "xi")
  return(result)
}
## HYP Fitting
hypfits <- apply(yrets, 2, fit.hypuv, symmetric = FALSE)
points <- matrix(unlist(lapply(hypfits, xichi)),
                 ncol = 2, byrow = TRUE)
## Shape triangle
col.def <- c("black", "blue", "red", "green", "orange")
leg.def <- paste(rd, rep("day return", 5))
plot(points, ylim = c(-0.2, 1.2), xlim = c(-1.2, 1.2),
     col = col.def, pch = 16, ylab = expression(xi),
     xlab = expression(chi))
lines(x = c(0, -1), y = c(0, 1))
lines(x = c(0, 1), y = c(0, 1))
lines(x = c(-1, 1), y = c(1, 1))
legend("bottomright", legend = leg.def, col = col.def, pch = 16)
text(x = 0.0, y = 1.05, label = "Laplace", srt = 0)
text(x = -1.0, y = 1.05, label = "Exponential", srt = 0)
text(x = 1.0, y = 1.05, label = "Exponential", srt = 0)
text(x = 0.0, y = -0.1, label = "Normal", srt = 0)
text(x = -0.6, y = 0.5, label = "Hyperbolic, left skewed",
     srt = 302)
text(x = 0.6, y = 0.5, label = "Hyperbolic, right skewed",
     srt = 57)
abline(v = 0, col = "grey")



## Ch 6.4 back-test
## Loading of packages
library(lmomco)
library(FRAPO)
## Data loading
data(SP500)
Idx <- SP500[, "QCOM"] #length265
L <- -1 * returnseries(Idx, method = "discrete", trim = TRUE) # percentaged return
## Computing VaR (Normal & GLD) 99%, moving window 
ep <- 104:length(L) #end points
sp <- 1:length(ep) #start points
level <- 0.99
VaR <- matrix(NA, ncol = 2, nrow = length(ep))
for(i in 1:length(sp)){
  x <- L[sp[i]:ep[i]]
  lmom <- lmom.ub(x)
  fit <- pargld(lmom)
  VaRGld <- quagld(level, fit)
  VaRNor <- qnorm(level, mean(x), sd(x))
  VaR[i, ] <- c(VaRGld, VaRNor)
  print(paste("Result for", ep[i], ":", VaRGld, "and", VaRNor)) 
}
## Summarising results
Res <- cbind(L[105:length(L)], VaR[-nrow(VaR), ])
colnames(Res) <- c("Loss", "VaRGld", "VaRNor")
## Plot of backtest results
plot(Res[, "Loss"], type = "p", xlab = "Time Index",
     ylab = "Losses in percent", pch = 19, cex = 0.5,
     ylim = c(-15, max(Res)))
abline(h = 0, col = "grey")
lines(Res[, "VaRGld"], col = "blue", lwd = 2)
lines(Res[, "VaRNor"], col = "red", lwd = 2)
legend("bottomleft", legend = c("Losses", "VaR GLD", "VaR Normal"),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1), pch = c(19, NA, NA), bty = "n")


## Ch 6.5
library(FRAPO)
library(fBasics)
## Loading of data
data(INDTRACK3)
P <- INDTRACK3[, -1]
R <- returnseries(P, method = "discret", trim = TRUE)
## Fitting and calculating beta and lambda
Fit <- apply(R, 2, gldFit, method = "rob", doplot = FALSE,
             trace = FALSE)
DeltaBetaParam <- matrix(unlist(lapply(Fit, function(x){
  l <- x@fit$estimate[c(3, 4)]
  res <- c(l[2] - l[1], l[1] + l[2])
  res})), ncol = 2, byrow = TRUE)
## Shape triangle
plot(DeltaBetaParam, xlim = c(-2, 2), ylim = c(-2, 0),
     xlab = expression(delta == lambda[4] - lambda[3]),
     ylab = expression(beta == lambda[3] + lambda[4]),
     pch = 19, cex = 0.5)
segments(x0 = -2, y0 = -2, x1 = 0, y1 = 0,
         col = "grey", lwd = 0.8, lty = 2)
segments(x0 = 2, y0 = -2, x1 = 0, y1 = 0,
         col = "grey", lwd = 0.8, lty = 2)
segments(x0 = 0, y0 = -2, x1 = 0, y1 = 0, col = "blue",
         lwd = 0.8, lty = 2)
segments(x0 = -0.5, y0 = -0.5, x1 = 0.5, y1 = -0.5,
         col = "red", lwd = 0.8, lty = 2)
segments(x0 = -1.0, y0 = -1.0, x1 = 1.0, y1 = -1.0,
         col = "red", lwd = 0.8, lty = 2)


## Ch 8.1
library(AER)
library(fGarch)
data(NYSESW)
NYSELOSS <- timeSeries(-1.0 * diff(log(NYSESW)) * 100,  #loss , -1
                       char.vec = time(NYSESW))
## Function for ES of t-GARCH
ESgarch <- function(y, p = 0.99){
  gfit <- garchFit(formula = ~garch(1, 1), data = y,
                   cond.dist = "std", trace = FALSE)
  sigma <-  predict(gfit, n.ahead = 1)[3] # $std
  df <- coef(gfit)["shape"]
  ES <- sigma * (dt(qt(p, df), df)/(1 - p)) *
    ((df + (qt(p, df))^2)/(df - 1))
  return(ES)
}
## Date vectors for backtest
from <- time(NYSELOSS)[-c((nrow(NYSELOSS) - 999) : nrow(NYSELOSS))]
to <- time(NYSELOSS)[-c(1:1000)]
NYSEES <- fapply(NYSELOSS, from = from, to = to, FUN = ESgarch)
NYSEESL1 <- lag(NYSEES, k = 1)
res <- na.omit(cbind(NYSELOSS, NYSEESL1))
colnames(res) <- c("NYSELOSS", "ES99")
plot(res[, 2], col = "red", ylim = range(res),
     main = "NYSE: t-GARCH(1,1) ES 99%",
     ylab = "percentages", xlab = "")
points(res[, 1], type = "p", cex = 0.2, pch = 19, col = "blue")
legend("topleft", legend = c("Loss", "ES"),
       col = c("blue", "red"), lty = c(NA, 1), pch = c(19, NA))


## Ch 9.1
library(QRM)
library(fGarch)
## Losses
data(EuStockMarkets)
loss <-  as.data.frame(na.omit(-1.0 * diff(log(EuStockMarkets)) *
                                 100.0))
gain <-  as.data.frame(na.omit(1.0 * diff(log(EuStockMarkets)) *
                                 100.0))
## GARCH
gfit <- lapply(loss, garchFit,  formula = ~ garch(1,1),
               cond.dist = "std", trace = FALSE)
gprog <- unlist(lapply(gfit, function(x)
  predict(x, n.ahead = 1)[3]))
gshape <- unlist(lapply(gfit, function(x) x@fit$coef[5]))
gresid <- as.matrix(data.frame(lapply(gfit,
                                      function(x) x@residuals / sqrt(x@h.t))))
## Copula
U <- sapply(1:4, function(y) pt(gresid[, y], df = gshape[y]))# t-distribution
cop <- fit.tcopula(Udata = U, method = "Kendall")
rcop <- rcopula.t(100000, df = cop$nu, Sigma = cop$P)
qcop <- sapply(1:4, function(x) qstd(rcop[, x], nu = gshape[x]))
ht.mat <- matrix(gprog, nrow = 100000, ncol = ncol(loss),
                 byrow = TRUE)
pf <- qcop * ht.mat
## ES 95 percent
weights <- c(0.4, 0.2, 0.2, 0.2)
pfall <- (qcop * ht.mat) %*% weights
pfall.es95 <- median(tail(sort(pfall), 5000))
pfall.es95
