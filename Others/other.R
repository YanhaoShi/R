library("PerformanceAnalytics")
library("fPortfolio")
library("fAssets")
library("xts")
library("MarkowitzR")
library("quadprog")
frontierPlot()


is.timeSeries(SP500)
# fAssets package
SPzoo <- timeSeries(zoo(SP500, order.by=as.Date(row.names(SP500), format= "%Y-%m-%d")))
assetsHistPlot(SPzoo[,1])
assetsLogDensityPlot(SPzoo[,1])
assetsHistPairsPlot(SPzoo[,1:2])
assetsBoxPlot(SPzoo[,1])
assetsBoxPercentilePlot(SPzoo[,1])
assetsQQNormPlot(SPzoo[,1])

assetsTreePlot(SPzoo[,1:50])
assetsStarsPlot(SPzoo[,1:5])

frontierPlot(SPzoo[,1], return = "mean", risk = "Sigma")


data(assay)
class(assay)
data(dow.jan.2005)
map.market(id = dow.jan.2005$symbol,
           area = dow.jan.2005$price,
           group = dow.jan.2005$sector,
           color = 100 * dow.jan.2005$month.ret)




# portfolio.optimize function--ghyp
data(indices)
ind.mu <- unname(apply(indices, 2, mean))
ind.cov <- unname(cov(indices))
sapply(diag(ind.cov), function(x) Porgr(ind.cov, ind.mu, x, as.matrix(rep(0,5)),
                                        as.matrix(rep(1,5))))


t.object <- fit.tmv(-indices, silent = TRUE)
# given a targeted return
t.ptf1 <- portfolio.optimize(t.object,
                             risk.measure = "sd",
                             type = "target.return",
                             target.return = 0.00859,
                             distr = "return",
                             silent = TRUE)