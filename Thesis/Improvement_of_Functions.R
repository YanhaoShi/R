if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
#############################################################################
Ver.info #view version information
verind

# Example: the first 20 assets
assets20 <- GetIndex(1:20, assets1457)
assets50 <- GetIndex(1:50, assets1457)
result20.list <- lapply(CLA, function(fun) 
          fun(assets20$mu, assets20$covar, assets20$lB, assets20$uB))
result50.list <- lapply(CLA, function(fun) 
  fun(assets50$mu, assets50$covar, assets50$lB, assets50$uB))

############################################################################
# (1) PurgeNumErr 
# compare function from Ver 1.2 & Ver 2
# test using results from Ver 2

# return indices purged by PurgeNumErr function from Ver1.2 and Ver 2
testPurgeNumErr(result50.list$Env2, assets50$lB, assets50$uB) #??where is the function?

# index change from each step (NA for first and last step)
ncol(result50.list$Env2$weights_set)

pdf("D:/Master Thesis/MyThesis/Pictures/purgeNumErr_wsum20.pdf",
    width = 10, height = 5) 
plot_weight(result50.list$Env2, assets50$lB, assets50$uB, plot.show = "sum")
abline(h = log(1e-10), col = "grey")
dev.off()

pdf("D:/Master Thesis/MyThesis/Pictures/purgeNumErr_wf20.pdf",
    width = 10, height = 5) 
plot_weight(result50.list$Env2, assets50$lB, assets50$uB, plot.show = "f")
dev.off()

############################################################################
# (2) purgeChull/ purgeSlope(best) / purgeChull
# compare purgeExcess(Ver 1.2), purgeChull(Ver 3) & purgeSlope(Ver 4)
# test using results from Ver 2
ind20 <- testPurgeChull(result20.list$M2, assets20$mu, assets20$covar) 
#return ind kept

# visualization
plot_EF_purge(result50.list$Env2, assets50$mu, assets50$covar, 
        EF = c("unpurge", "Excess", "Chull"), Smoo = TRUE)

# when weights sum to 1, EF is convex naturally!
plot_EF_purge(result50.list$Env3, assets50$mu, assets50$covar) 

pdf("D:/Master Thesis/MyThesis/Pictures/purgeEF_20.pdf",
    width = 10, height = 7) 
plot_EF_purge(result20.list$Env2, assets20$mu, assets20$covar, 
        EF = c("Chull", "Excess", "unpurge"), point.label = TRUE, Smoo = TRUE)
dev.off()
 ####    #####
uni_purge_sp <- readRDS(d.file("purge_unique_sp.rds", exists = FALSE))

pdf("D:/Master Thesis/MyThesis/Pictures/purgeunique.pdf",
    width = 10, height = 7) 
plot_unipurge(uni_purge_nasdaq, "NASDAQ")
dev.off()


#############################################################################
# (3) solve() from Ver 2 and chol() from Ver 1.2
# compare time cost 
pdf("D:/Master Thesis/MyThesis/Pictures/covF_nasdaq.pdf",
    width = 10, height = 7) 
plot_Inv(covF.nasdaq, "NASDAQ")
dev.off()


pdf("D:/Master Thesis/MyThesis/Pictures/covF_company.pdf",
    width = 10, height = 7) 
IndexPlot(analysis.company, "company", ind.company, 1457,  
                      show_plot = c("nwbox", "nwplot", "kappa",
                                    "rcond"))
dev.off()

#############################################################################
# lambda
pdf("D:/Master Thesis/MyThesis/Pictures/EF_lambda.pdf",
    width = 10, height = 7) 
re7 <- Env7$cla.solve(assets50, purge=TRUE)
re7.1 <- Env7.1$cla.solve(assets50, purge=TRUE, lam.tol = 1e-14)
plot_lambda(re7.1, re7,  assets50$mu, assets50$covar, smoo = TRUE)
dev.off()



#############################################################################
# (5) Time Improvements
# improving.....
company.cla.mt1 <- readRDS(d.file("company_cla_mt1.rds"))
company.cla.mt7 <- readRDS(d.file("company_cla_mt7.rds"))
nasdaq.cla.mt1 <- readRDS(d.file("nasdaq_cla_mt1.rds"))
nasdaq.cla.mt7 <- readRDS(d.file("nasdaq_cla_mt7.rds"))
sp.cla.mt7 <- readRDS(d.file("sp_cla_mt7.rds"))
sp.cla.mt1 <- readRDS(d.file("sp_cla_mt1.rds"))


pdf("D:/Master Thesis/MyThesis/Pictures/nasdaq_tcompare.pdf",
    width = 10, height = 7) 
plot_cla_tcompare(nasdaq.cla.mt7[1:18, ], nasdaq.cla.mt1[1:18, ], scale.log = T)
dev.off()

pdf("D:/Master Thesis/MyThesis/Pictures/company_tcompare.pdf",
    width = 10, height = 7) 
plot_cla_tcompare(company.cla.mt7, company.cla.mt1, scale.log = T)
dev.off()

t_lm(list(t1 = log(nasdaq.cla.mt1[,1]), t2 = log(nasdaq.cla.mt7[,1])))




pdf("D:/Master Thesis/MyThesis/Pictures/index_t7.pdf",
    width = 10, height = 7) 
mt1 <- log(company.cla.mt7[, 1])
mt2 <- log(sp.cla.mt7[, 1])
mt3 <- log(nasdaq.cla.mt7[, 1])
l1 <- length(mt1)
l2 <- length(mt2)
l3 <- length(mt3)
lmax <- max(l1, l2, l3)
plot(seq(l1), mt1, type = "l", col = adjustcolor("red", 0.5), 
     pch = 16, xaxt = "n", main = "Log Time Comarison (improved version)", 
     xlab = "number of assets", ylab = "time(s)", 
     ylim = range(c(mt1, mt2, mt3)), xlim = c(1, lmax))
axis(1, at = seq(lmax), labels = names(mt3))
lines(seq(l2), mt2,  col = adjustcolor("green", 0.8))
lines(seq(l3), mt3, col = "blue4")
legend("topleft", legend = c("company", "sp", "nasdaq"), lwd = 1, 
       col = c(adjustcolor("red", 0.5), adjustcolor("green", 0.8), "blue4"))
dev.off()
t_lm(list(t1 = log(company.cla.mt7[,1]), t2 = log(sp.cla.mt7[,1]),
          t3 = log(nasdaq.cla.mt7[, 1])))


result10 <- Env7$cla.solve(GetIndex(1:10, assets1457, lB = rep(0,10), uB = rep(0.5,10)))
matplot(t(result10$weights_set), type ="o")






