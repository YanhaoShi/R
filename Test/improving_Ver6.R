if(!file.exists("main.R")) stop("R directory does *not* contain main.R")
source("main.R")
#############################################################################
assets <- GetIndex(1:50, assets1457)

lapply((1:5)*50, function(n){ assets <- GetIndex(1:n, assets1457)
all.equal(Env7$cla.solve(assets, purge = F), Env8$cla.solve(assets, purge = F))
}
)
all.equal(Env7$cla.solve(assets, purge = F), 
          Env8$cla.solve(assets, purge = F))




lapply((2:10)*10, function(n){ assets <- GetIndex(1:n, assets1457)
all.equal(Env7$cla.solve(assets, purge = F), Env7.2$cla.solve(assets, purge = F))
}
)

sapply(1:3, function(x) c(x,x))


microbenchmark(Env7$cla.solve(assets, purge = F),
Env7.2$cla.solve(assets, purge = F), times = 5)

system.time(Env7$cla.solve(assets, purge = F))
system.time(Env8$cla.solve(assets, purge = F))

computeInv(getMatrices(assets$mu, assets$covar, 1:50, 1:24))

m <- matrix(1:9, 3, 3)
m1 <- matrix(sample(1:1000000,2500),50,50)
microbenchmark(assets$mu[1:10],
assets$mu[-(1:10)])

assets <- GetIndex(1:10, assets1457)
get <- getMatrices(assets$mu, assets$covar, rep(1,10), 1:4)
i <- 5
s <- solve(getMatrices(assets$mu, assets$covar, rep(1,10), 1:5)$covarF)
covar <- assets$covar; f <-1:4; mu <- assets$mu
solve(get$covarF, cbind(1, get$muF, get$covarFB %*% get$wB, deparse.level=0L))
####################
computeInv = function(covF.inv.pre , f, i, add = FALSE){
  fi <- c(f, i) # add = TRUE
  a <- covar[f, i]
  cc <- covF.inv.pre %*% a
  b <- covar[i, i] - sum(a * cc)
  e1 <-  cc%*%t(cc)/b
  m <- cc/b
  cov.inv <- rbind(cbind(covF.inv.pre + e1, -m),
                   cbind(-t(m), 1/b))
  
  fi <- f
  cov.inv %*% cbind(1, mu[fi], covar[fi, -fi] %*% w[-fi], deparse.level=0L)
}
#####################3
w<- rep(1,10)
all.equal(s, s2)
(inv1 <- computeInv(getMatrices(assets$mu, assets$covar, rep(1,10), c(1:5))))
inv.pre <- computeInv(getMatrices(assets$mu, assets$covar, rep(1,10), 1:4))
m%*%c(1,1,1)

rowSums(m)
###################
#solve(covF+, 1):
rbind(inv.pre[,1,drop =F]+rowSums(beta* cc%*% t(cc)) -beta*cc,
      -beta*sum(cc) + beta)

#solve(  , muF)
muF <- mu[1:4]
rbind(inv.pre[,2,drop=F] + e1%*%mu[1:4] + e2 *mu[i],
      t(e2)%*%muF + beta * mu[i])
#solve(, covFB wB)




