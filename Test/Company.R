if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
########################################################################
# company
# 1457 assets
analysis.company <- IndexAnalysis(Index.trans = assets1457, 
                             name = "company", n = 591)

# latex Appendix B
tab <- t(analysis.company$t_mean)
xtable(tab) # add "assets"!!
#######################################################################################
# plot
IndexPlot(analysis.company, name = "company", n = 591)  # open sp_plot.pdf




ans <- time_lm(analysis.company, n = 591, vector.predict = c(999, 100))

plot_fit(ans$fit, analysis.company, "company", n = 591, n.predict = 1457) 















##draft
for(i in 12:14){  ##change to 15-18
  for(j in 1:10){
    nas <- readRDS(d.file(paste("Company", length(ind.company[[i]][[1]]), "_set", j, sep = ""),exists = F))
    saveRDS(nas, d.file(paste("Company", length(ind.company[[i]][[1]]), "_set", j, ".rds", sep = ""),exists = F))
  }
}
l.company = 11
label.company <- label.company[1:11]
###

# qp method is faster than cla when the number of assets is small, (55-95)
# but much slower when number of assets is larger (114-410)


## draft here

m <- covF[[2]][[1]][[22]]
kappa(m)

rcond(m)
plot(cov.kappa, type = "o", main = "kappa")

cov.m <- lapply(cf, as.matrix)
sapply(m, rcond)
lapply(cov.m, nearPD)

