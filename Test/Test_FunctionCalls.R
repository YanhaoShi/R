if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
########################################################################

#Created data for validating
s <- validate(10,20)

########################################################################
#Compare list of results of different classes
compare<-function(l1,l2){    
  n<-length(l2)
  for(i in 1:n){
    unl1 <- unlist(l1[[i]])
    unl2 <- unlist(l2[[i]])
    stopifnot(all.equal(unl1, unl2, check.attributes = FALSE) ) 
  }
  TRUE
}

#initAlgo
all.equal(initAlgo(s$mu,data.frame(s$lB),data.frame(s$uB)),
        initAlgo(s$mu,s$lB,s$uB),
        check.attributes = FALSE)

compare(initAlgo(s$mu,data.frame(s$lB),data.frame(s$uB)),
          initAlgo1(s$mu,s$lB,s$uB))

#getB
all.equal(getB(s$mu,s$f),getB1(s$mu,s$f))

#reduceMatrix
all.equal(reduceMatrix(s$covar, s$meanF, s$f),
          reduceMatrix1(s$covar, s$meanF, s$f), 
          check.attributes = FALSE  )

all.equal(reduceMatrix(data.frame(s$mu), s$f, 0),
          reduceMatrix1(s$mu, s$f),
          check.attributes = FALSE)

#getMatrices
all.equal(getMatrices(as.matrix(s$mu), s$covar, s$solution_set, s$f), 
          getMatrices1(s$mu, s$covar, s$solution_set, s$f),
          check.attributes = FALSE)

#getMinVar (later)
#compare(getMinVar(covar,data.frame(s$solution_set)),getMinVar1(covar,s$solution_set))

#efFrontier (later)

#computeW 
cw <- computeW(s$lambdas,s$covarF_inv,s$covarFB,s$meanF,s$wB)
cw1 <- computeW1(s$lambdas, s$covarFB, s$meanF, s$wB, s$covarF)
cw2 <- computeW(s$lambdas, s$covarFB, s$meanF, s$wB, s$covarF)
all.equal(cw, cw1, check.attributes = FALSE)
all.equal(cw, cw2, check.attributes = FALSE)

#computeLambda
cl <- computeLambda(solve(s$covarF), s$covarFB, s$meanF, s$wB, i = 2, bi = 1:2)
cl1 <- computeLambda(solve(s$covarF), s$covarFB, s$meanF, s$wB, i = 2, bi = 1:2)
cl2 <- computeLambda2(s$covarFB, s$meanF, s$wB,i = 2, bi = 1:2, s$covarF)

all.equal(cl, cl1, check.names = FALSE)
all.equal(cl, cl2, check.names = FALSE)

 #compare a list of cl2 and cl3
i <- 1:5
b <- list(lB = 2:6, uB = 10:14)
cl2_list <- list()
for(i in 1:5){
  cl2_list$lambda[i] <- computeLambda2(s$covarFB, s$meanF, s$wB,
                                      i, list(b$lB[[i]],b$uB[[i]]), s$covarF)$lambda
  cl2_list$bi[i] <- computeLambda2(s$covarFB, s$meanF, s$wB,
                                 i, list(b$lB[[i]],b$uB[[i]]), s$covarF)$bi
}
cl3 <- computeLambda3(s$covarFB, s$meanF, s$wB, i = 1:5, b, s$covarF)

all.equal(cl2_list, cl3)

# purgeNummErr1 and purgeNummErr2 get same results
index1 <- purgeNummErr1(assets[[50]]$lB, assets[[50]]$uB, 
                        as.data.frame(r2_50_unpurge$weights_set), 1e-9, 1e-10)

index2 <- purgeNummErr2(assets[[50]]$lB, assets[[50]]$uB,
                        as.data.frame(r2_50_unpurge$weights_set), 1e-9, 1e-10)

pur1 <- list( r2_50_unpurge$weights_set[, -index1], r2_50_unpurge$free_indices[-index1],
              r2_50_unpurge$lambdas[-index1], r2_50_unpurge$gammas[-index1])

pur2 <- list(r2_50_unpurge$weights_set[, index2], r2_50_unpurge$free_indices[index2],
             r2_50_unpurge$lambdas[index2], r2_50_unpurge$gammas[index2])

all.equal(pur1, pur2)







