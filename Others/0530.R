CLA <- list(
  M1 = function(mean, covar, lB, uB){  # Alexander's version
    source("Functions/Ver1_previous.R")  #change the working directory
    Method1$cla.solver(mean, covar, lB, uB)
  },
  M2 = function(mean, covar, lB, uB){  # my version as benchmark
    source("Functions/Ver2_benchmark.R")  #change the working directory
    Method2$cla.solver(mean, covar, lB, uB)
  },
  M3 = function(mean, covar, lB, uB){  # the fastest version so far
    source("Functions/Ver3.R")  #change the working directory
    Method3$cla.solver(mean, covar, lB, uB)
  })

## To call a cla version
## eg. for first 50 assets:
assets50 <- readRDS(d.file("assets50.rds", exists = FALSE))

## call cla version 3
results <- CLA$M3(assets50$mu, assets50$covar, assets50$lB, assets50$uB)

names(results)
# after purge: "weights_set_purge", "free_indices_purge", "gammas_purge", "lambdas_purge"
# before purge: "weights_set", "free_indices"      "gammas"     "lambdas"
# covarF: to see the condition number of covarF


CLA_time <- readRDS(d.file("CLA_time.rds", exists = FALSE))  #change the working directory

t <- CLA_time[c(1, 51, 72:90),]
nrow(t)
plot((1:21)*50, t[,1], type = "l", pch = 16, xlab = "number of assets", ylab = "milliseconds")
title("Efficient Frontier for CLA Versions")
lines((1:21)*50, t[,2], col = "red")
lines((1:21)*50, t[,3], col = "blue")
legend("topleft", legend = c("Ver1", "Ver2", "Ver3"), lwd = 1, col = c("black", "red", "blue"))

#sessionInfo()
#R version 3.3.2 (2016-10-31)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 7 x64 (build 7601) Service Pack 1

#locale:
#[1] LC_COLLATE=German_Switzerland.1252  LC_CTYPE=German_Switzerland.1252   
#[3] LC_MONETARY=German_Switzerland.1252 LC_NUMERIC=C                       
#[5] LC_TIME=German_Switzerland.1252   