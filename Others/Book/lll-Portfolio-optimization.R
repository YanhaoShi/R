#Ch 10.1
## Loading of packages
library(copula)
library(quadprog)
library(rrcov)
## Creating copula objects
ncop <- normalCopula(param = 0.5, dim = 5)
tcop <- tCopula(param = 0.5, dim = 5, df = 5, df.fixed = TRUE)
## Creating DGPs(data-generating processes)
NcopMargN <- mvdc(ncop, margins = "norm",
                  paramMargins = list(list(mean = 0, sd = 1)),
                  marginsIdentical = TRUE)
NcopMargT <- mvdc(ncop, margins = "t",
                  paramMargins = list(df = 5),
                  marginsIdentical = TRUE)
TcopMargT <- mvdc(tcop, margins = "t",
                  paramMargins = list(df = 5),
                  marginsIdentical = TRUE)
## Initialising list objects for DGP
Lobj <- list()
length(Lobj) <- 1000
## Setting a seed
set.seed(12345)
## Generating random samples
rNcopMargN <- lapply(Lobj, function(x) rMvdc(240, NcopMargN))
rNcopMargT <- lapply(Lobj, function(x) rMvdc(240, NcopMargT))
rTcopMargT <- lapply(Lobj, function(x) rMvdc(240, TcopMargT))
# 1000 samples of each DGPs with 240 rows and five columns


##Ch 10.2
## Function for Moment Estimation
Moments <- function(x, method = c("CovClassic", "CovMcd",
                                  "CovMest", "CovMMest", "CovMve", "CovOgk",
                                  "CovSde", "CovSest"), ...){
  method <- match.arg(method)
  ans <- do.call(method, list(x = x, ...))
  return(getCov(ans))
}



##Ch 10.3
## Dimensions of Simulation
DGP <- c("rNcopMargN", "rNcopMargT", "rTcopMargT")
EST <- c("CovClassic", "CovMcd", "CovMest", "CovMMest",
         "CovMve", "CovOgk", "CovSde", "CovSest")
SAMPLE <- c(60, 120, 240)
## Creating list objects for combinations of
## DGP and sample sizes
## initialising vector for data objects
datnames <- NULL
for(i in DGP){
  for(j in SAMPLE){
    objname <- paste(i, j, sep = "")
    datnames <- c(datnames, objname)
    cat(paste("Creating list object", objname, "\n"))
    assign(objname, lapply(eval(as.name(i)), function(x) x[1:j, ]))
  }
}
## Creating list objects with estimates of 
## location and dispersion for combinations of
## DGP, sample sizes and estimators
## initialising vector for list objects
objnames <- NULL
for(i in datnames){
  for(j in EST){
    objname <- paste(j, i, sep = "")
    objnames <- c(objnames, objname)
    cat(paste("Creating list object", objname, "\n"))
    assign(objname, lapply(eval(as.name(i)), Moments, method = j))
  }
}



##Ch 10.4
## Function for minimum-variance portfolio
## Constraints: Fully invested, long-only
PortMinVar <- function(x){
  Dmat <- x   # Matrix to minimize: covariance matrix
  k <- ncol(Dmat)
  dvec <- rep.int(0, k)
  a1 <- rep.int(1, k) # weights sum up to one
  b1 <- 1
  a2 <- diag(k)
  b2 <- rep.int(0, k) # weights non-negative
  Amat <- t(rbind(a1, a2))
  bvec <- c(b1, b2)
  opt <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat,
                  bvec = bvec, meq = 1)
  return(opt$solution)
}
## Conduct optimisation
portnames <- NULL
idx <- 1:1000
for(i in objnames){
  objname <- paste("Port", i, sep = "")
  portnames <- c(portnames, objname)
  obj <- eval(as.name(i))
  weights <- lapply(obj, PortMinVar)
  assign(objname, sapply(idx, function(x)
    sqrt(t(weights[[x]]) %*% obj[[x]] %*%
           weights[[x]])))
}
## Caluculate median and IQR of portfolio risks
mednames <- NULL
iqrnames <- NULL
for(i in portnames){
  objname1 <- paste("Med", i, sep = "")
  objname2 <- paste("IQR", i, sep = "")
  mednames <- c(mednames, objname1)
  iqrnames <- c(iqrnames, objname2)
  assign(objname1, median(eval(as.name(i))))
  assign(objname2, IQR(eval(as.name(i))))
}



## Ch 10.5
## Loading of packages
library(FRAPO)
library(PerformanceAnalytics)
library(quadprog)
library(rrcov)
library(zoo)
## Loading data, calculate returns
data(StockIndex)
pzoo <- zoo(StockIndex, order.by = rownames(StockIndex)) #organized ts
rzoo <- (pzoo / lag(pzoo, k = -1) - 1) * 100 # return (percentage)
## boxplot and descriptive statistics
boxplot(coredata(rzoo), ylab = "return") #多个assets用不上boxplot
abline(h=0,col = "grey")

rstats <- rbind(apply(rzoo, 2, summary),
                skewness(rzoo),
                kurtosis(rzoo)
)
rstats


## Ch 10.6
## Conduct back test
EST <- c("CovClassic", "CovMcd", "CovMest", "CovMMest", "CovMve",
         "CovOgk", "CovSde", "CovSest")
## Function for back test
PfBack <- function(x, method = c("CovClassic", "CovMcd", "CovMest",
                                 "CovMMest", "CovMve", "CovOgk", "CovSde",
                                 "CovSest"), ...){
  cov <- Moments(x, method = method)   ## Moments() Ch 10.2
  return(PortMinVar(cov)) ## PortMinVar() Ch 10.4
}
## Conducting back test
PfWeights <- lapply(EST, function(x)
  rollapply(rzoo, width = 120, FUN = PfBack,
            method = x, by.column = FALSE,
            align = "right"))

periods <- as.Date(index(PfWeights[[1]]))
## Calculate portfolio returns / relative performance
PfReturns <- lapply(PfWeights, function(x)
  rowSums(lag(x, k = -1) * rzoo))
PfReturns <- zoo(matrix(unlist(PfReturns),
                        ncol = length(PfReturns)), periods)
colnames(PfReturns) <- EST
PortOut <- (PfReturns[, -1] - PfReturns[, 1])
## plot relative peformance
plot(PortOut, type = "h",
     xlab = "",
     ylab = EST[-1],
     main = "Relative Performance",
     ylim = range(PortOut))
## statistics on relative performance
PortRelStats <- rbind(apply(PortOut, 2, summary),
                      skewness(PortOut)
)
PortRelStats 



## Ch 10.7
## Defining function for points on efficient frontier
PMV <- function(SRoot, mu, SigTerm,   #SRoot: square root of cov matrix
                optctrl = ctrl(trace = FALSE)){
  N <- nrow(SRoot)
  ## Portfolio risk constraint
  soc1 <- socc(F = SRoot, g = rep(0, N), d = rep(0, N), f = SigTerm)
  ## non-negativity constraint
  nno1 <- nnoc(G = -diag(N), h = rep(0, N))
  ## Budget constraint
  A1 <- matrix(rep(1, N), nrow = 1)
  b1 <- 1.0
  ## optimization
  ans <- cccp(q = -mu, A = A1, b = b1, cList = list(nno1, soc1),
              optctrl = optctrl)
  getx(ans)
}



## Ch 10.8
library(cccp)
## Setting of parameters
Nassets <- ncol(rzoo)  
Nobs <- nrow(rzoo)
mu <- colMeans(rzoo)
S <- cov(rzoo)
SR <- sqrm(S)
delta <- sqrt(qchisq(0.9, Nassets)) # 90% quantile for the returns uncertainties
## Determining feasible risk aversion
SigMax <- max(colSds(rzoo))
SigMin <- min(colSds(rzoo))

#risk aversion parameter
ra <- seq(SigMin * 1.1, SigMax * 0.9, length.out = 10) / SigMax  ##remove outlier, rescale
## Initializing objects for classical mean-variance(MV)
## and robust counterpart(RC) results
RCans <- MVans <- matrix(NA,
                         nrow = 10,
                         ncol = Nassets + 2)
## Computing points on efficient frontier and allocations
for(i in 1:10){
  ## minimum-variance
  wmv <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / ra[i]) # MV weights
  MVans[i, ] <- c(sqrt(t(wmv) %*% S %*% wmv), # each row contains risk, return and weights
                  crossprod(mu, wmv),
                  wmv)
  ## robust counterpart
  theta <- ra[i] + (1 - ra[i]) * delta / sqrt(Nobs)
  wrc <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / theta) # RC weights
  RCans[i, ] <- c(sqrt(t(wrc) %*% S %*% wrc),
                  crossprod(mu, wrc),
                  wrc)
}



## Ch 10.9
## Equivalent weighting parameter for theta
theta2lambda <- function(theta, Nassets, Nobs, level){
  delta <- sqrt(qchisq(level, Nassets))
  lambda <- theta / (1 + theta * (delta / sqrt(Nobs)))
  lambda
}
## robust allocation and equivalent
## mean-variance allocation
theta <- 0.7
wrc <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / theta)
## RC point on efficient frontier
rceq <- c(sqrt(t(wrc) %*% S %*% wrc), crossprod(mu, wrc))
## Equivalent risk weighting
rweq <- theta2lambda(theta, Nassets = Nassets, Nobs = Nobs, level = 0.9)
## Equivalent MV point on efficient frontier
wmv <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / rweq)
mveq <- c(sqrt(t(wmv) %*% S %*% wmv), crossprod(mu, wmv))




## Ch 10.10
## Efficient Frontier
## determining ranges
xlims <- c(SigMin * 0.9,
           max(cbind(MVans[, 1], RCans[, 1])))
ylims <- c(min(cbind(MVans[, 2], RCans[, 2])),
           max(cbind(MVans[, 2], RCans[, 2])))
## plotting efficient frontier for MV
plot(MVans[, 1], MVans[, 2], type = "l",
     xlim = xlims,
     ylim = ylims,
     xlab = expression(sigma),
     ylab = expression(mu))
## Superimposing points
for(i in 1:nrow(MVans)){
  points(x = MVans[i, 1], y = MVans[i, 2], col = "blue", pch = 19)
  points(x = RCans[i, 1], y = RCans[i, 2], col = "red", pch = 19)
}
## Superimposing equivalence points
points(x = rceq[1], y = rceq[2], col = "darkgreen", bg = "darkgreen", pch = 23)
points(x = mveq[1], y = mveq[2], col = "green", bg = "green", pch = 23)
## Legend
legend("topleft", legend = c("Efficient Frontier", "MV points", "RC points",
                             expression(paste("RC allocation with ", theta == 0.7)),
                             "Equivalent MV-Portfolio"),
       lty = c(1, NA, NA, NA, NA), pch = c(NA, 19, 19, 23, 23),
       pt.bg = c(NA, NA, NA, "darkgreen", "orange"),
       col = c("black", "blue", "red", "darkgreen", "orange"))

