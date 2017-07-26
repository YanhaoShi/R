if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
#############################################################################
lam.tol <- -18:0
ntol <- length(lam.tol)

# compute list of lam.lmatrix
if(!file.exists(d.file("NASDAQ2116_lam.rds", exists = FALSE))){
  lam.lmatrix <- list() 
  lam.tol <- -18:0
  for(j in 1:10){ 
    lam.list <- lapply(lam.tol, function(i) 
      readRDS(d.file(paste0("NASDAQ2116_set", j, "_tol", i, ".rds")))$lambdas)
    names(lam.list) <- paste0("tol", lam.tol)
    lam.lmatrix[[j]] <- mLam.tol(lam.list)
  }
  saveRDS(lam.lmatrix, d.file("NASDAQ2116_lam.rds", exists = FALSE)) 
}
lam.lmatrix <- readRDS(d.file("NASDAQ2116_lam.rds", exists = FALSE)) 

#############################################################
lapply(1:10, function(j) plotLam.tol(lam.lmatrix[[j]]))

pdf("D:/Master Thesis/MyThesis/Pictures/lam_ana op++lysis.pdf",
    width = 10, height = 7) 
plotLam.tol(lam.lmatrix[[10]])
dev.off()

r.nas <- readRDS(d.file("NASDAQ2116_set10_tol0.rds"))
lam0 <- r.nas$lambdas
fchange <- f_change(r.nas$free_indices)
ind.equal <- which(testLam(lam0, equal.tol = -18)$lam.equal)

# latex: plot original lambdas from NASDAQset10
pdf("D:/Master Thesis/MyThesis/Pictures/lam_all.pdf",
    width = 10, height = 7) 
plot(lam0[-1], type = "l", ylab = expression(lambda), xlab = "CLA steps",
     main = "lambdas of result", ylim = c(-2, max(lam0[-1])))
points(ind.equal, lam0[ind.equal], pch =16, col = "red")
text(ind.equal, lam0[ind.equal], fchange[ind.equal], pos = c(1,3))
legend("topright", legend = "lambda in D", col = "red", pch = 16)
dev.off()

# latex: plot original lambdas from NASDAQset10, -/+, first 50 steps
ind50 <- ind.equal[ind.equal < 50]
ind50.b <- ind50 - 1
ind50.f <- ind50 + 1
ind.t <- c(ind50, ind50.b, ind50.f)

pdf("D:/Master Thesis/MyThesis/Pictures/lam_50.pdf",
    width = 10, height = 7) 
plot(lam0[2:51], type = "o", ylab = expression(lambda), xlab = "first 50 CLA steps",
     main = "lambdas of result", ylim = c(-2, max(lam0[2:51])),
     pch =16, col = "grey")
points(ind50, lam0[ind50], pch =16, col = "red")
points(ind50.b, lam0[ind50.b], pch =16, col = "blue")
points(ind50.f, lam0[ind50.f], pch =16, col = "green")
text(ind.t, lam0[ind.t], fchange[ind.t], pos = c(1,1,1,3,3,3,3,3,3))
legend("topright", 
       legend = c("lambda in D", "lambda in D-","lambda in D+", "others"), 
       col = c("red", "blue", "green", "grey"), pch = 16)
dev.off()



length(ind.equal); length(lam0) #latex 44/375
# same weight columns
all.equal(r.nas$weights_set[,ind.equal], r.nas$weights_set[,ind.equal +1])
all.equal(r.nas$weights_set[,ind.equal], r.nas$weights_set[,ind.equal -1])
# same gammas
all.equal(r.nas$gammas[ind.equal], r.nas$gammas[ind.equal +1])
all.equal(r.nas$gammas[ind.equal], r.nas$gammas[ind.equal -1])
lam0[1:10]

lam7 <- readRDS(d.file("NASDAQ2116_set10_tol-7.rds"))$lambdas
sum(duplicated(lam7))
setdiff(lam7, lam0)
setdiff(lam0, lam7)
all.equal(setdiff(lam0, lam7), lam0[ind.equal])
# proof that, purge at tol-7 doesnt harm different weigths columns
# proof EF is convex


lam7[-1]-lam7[length(lam7)]
