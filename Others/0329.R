## Part 1, my functions (instead of cla.solver1, but cla.solver.remove1)
# Initialize the weight -- Find first free weight
initAlgo1 <- function(mu, lB, uB ){
  #New-ordered return, lB, uB with decreasing return  
  w <- c()
  index.new <- order(mu,decreasing = T) # new order with decreasing return
  lB.new <- lB[index.new]
  uB.new <- uB[index.new]
  # free weight - starting solution
  i.new <- 0
  w.new <- lB.new # initialy
  while(sum(w.new) < 1) {
    i.new <- i.new + 1
    w.new[i.new] <- uB.new[i.new]
  }
  w.new[i.new] <- 1 - sum(w.new[-i.new])
  w[index.new] <- w.new                #back to original order
  i<-index.new[i.new] 
  list(index=i,weights=w)         # return the index of first free asset and vector w
} 

# getB1 ---  --------------------------------------------------------------------
getB1<-function(mean,f) (1:length(mean))[-f]

# reduceMatrix ---  ------------------------------------------------------------
reduceMatrix1 <- function(input, X, Y){
  if  (is.vector(input)){
    return(input[X])
  }
  else { return(input[X,Y,drop=F])}     
} 

# purgeNumErr ---  -------------------------------------------------------------
# Purge violations: remove w if sum(w)!=1, w<lB, w>uB
# tolerance, e.g 10e-10
purgeNummErr1 <- function(free_weights, gammas, lambdas, lB, uB, solution_set, tol=1e-9){ 
  i <- 1
  # k <- which(abs(colSums(solution_set)-1) > tol)
  # which(abs(colSums(re[[3]])-1) > tol)
  k <- c()
  l <- c()
  while(TRUE){# 
    flag <- FALSE
    if(i == ncol(solution_set)){
      return(list(solution_set,lambdas,gammas,free_weights))
      #return(list(k,l))
    }else{##
      if((abs(sum(solution_set[,i]))-1) > tol){ ###22
        flag <- TRUE
        
        k<-c(k,i)
      } else {
        for(j in seq(length(solution_set[,i]))){
          if((solution_set[j,i] - lB[j]) < (-tol) || (solution_set[j,i] - uB[j]) > (tol)){
            flag <- TRUE
            l <- c(l,i)
          }
        }
      } ###
      if(flag==TRUE){
        solution_set <- solution_set[,-i] 
        lambdas <- lambdas[-i]
        gammas <- gammas[-i]
        free_weights <- free_weights[-i] 
        
      }else{
        i <- i+1
      }
    }##
  }#
  return(list(solution_set,lambdas,gammas,free_weights,j))
  #return(list(k,l))
  
}

purgeNummErr1.draft <- function(free_weights, gammas, lambdas, lB, uB, solution_set1, tol=1e-9){
  # ncol=1时返回原值!!!!!
  k <- which(abs(colSums(solution_set1)-1) > tol)  ## k=1?
  
  for(i in 1:ncol(solution_set1)){
    #if(any(c(solution_set1[,i]-lB < -tol,solution_set1[,i]-uB > tol))){
    k<-c(which(solution_set1[,i]-lB < -tol),which(solution_set1[,i]-uB > tol),k)
    
    
  }
  k<- k[!duplicated(k)]
  solution_set1 <- solution_set1[,-k] 
  lambdas <- lambdas[-k]
  gammas <- gammas[-k]
  free_weights <- free_weights[-k] 
  
  return(list(solution_set1,lambda=lambdas,lgamma=gammas,weights_free=free_weights))
}

# purgeExcess ---   -------------------------------------------------------------
# Remove violations of the convex hull, ensure mu2>mu3, ...
purgeExcess1 <- function(mean, lambdas, gammas, free_weights, solution_set1){ 
  n <- ncol(solution_set1)
  mu <- t(mean) %*% solution_set1
  
  k<-which( mu[2:(n-1)] < mu[3:n] )  # 0 or tol?
  
  solution_set1 <- solution_set1[,-k] 
  lambdas <- lambdas[-k] 
  gammas <- gammas[-k]
  free_weights <- free_weights[-k] 
  
  return(list(solution_set1,lambdas, gammas, free_weights))
}

# getMinVar ---  ---------------------------------------------------------------
# Get the minimum variance solution
getMinVar1 <- function(covar,solution_set1){          ## getMinStd !
  n<-ncol(solution_set1)
  var <- rep(1,n)
  for(i in 1:n){
    w<-solution_set1[,i]
    var[i] <- t(w) %*% covar %*% w
  }
  
  list(index_minvar= sqrt(min(var)), weights_minvar=solution_set1[ ,which.min(var)] )
} 

# efFrontier -----------------------------------------------------------------
# Get the efficient frontier
efFrontier1 <- function(mean, covar, solution_set1, dataPoints){ 
  n<-ncol(solution_set1)
  mu <- c()
  sigma <- c()
  weights <- c()
  
  a <- seq(0,1,length=dataPoints+1)[-1]  #dataPoints-2 in between
  
  for(i in 1:(n-1)){
    w0 <- solution_set[i]
    w1 <- solution_set[i+1]
    
    weights <- w0%*%t(a)+w1 %*% t(1-a)   #weights between
    w <- cbind(weights,w)
    mu <- c(mu,colSums(weights))
    sigma <- c(sigma, sqrt(t(w) %*% covar %*% w))
  }
  return(list(mu, sigma, weights))
} 

# computeW -----------------------------------------------------------------
# compute gamma #(Bailey and Lopez de Prado 2013, 11)
computeW1 <- function(lambdas,covarF_inv,covarFB,meanF,wB){
  #1) compute gamma
  g1 <- sum(covarF_inv %*% meanF)  #rrr
  g2 <- sum(covarF_inv)
  if(is.null(wB)){
    g <- as.numeric(-lambdas[length(lambdas)]*g1/g2+1/g2) 
    w1 <- 0
  } else {
    g3 <- sum(wB)
    g4 <- covarF_inv %*% covarFB
    w1 <- g4 %*% wB
    g4 <- sum(w1)
    g <- as.numeric(-lambdas[length(lambdas)]*g1/g2+(1-g3+g4)/g2)
  }
  
  #2) compute weights
  w2 <- apply(covarF_inv,1,sum)
  w3 <- covarF_inv %*% meanF
  
  return(list(wF = -w1+g*w2+lambdas[length(lambdas)]*w3, gamma = g))
}


# computeLambda --------------------------------------------------------------
# (Niedermayer and Niedermayer 2007, 10)
computeLambda1 <- function(covarF_inv,covarFB,meanF,wB,i,bi){
  #1) C
  onesF <- matrix(1,length(meanF)) # (kx1)-vector
  c1 <- (t(onesF) %*% covarF_inv) %*% onesF
  c2i <- covarF_inv[i,] %*% meanF
  c3 <- t(onesF) %*% covarF_inv %*% meanF
  c4 <- covarF_inv %*% onesF
  c <- -c1*c2i+ c3*c4[i]  
  if(c==0){
    return(list(lambda = 0,bi = 0))
  }
  #2) bi
  if(is.list(bi)){
    
    bi <- ifelse(c>0, bi[[2]],bi[[1]])
  }
  #3) lambda
  if(length(wB)==0){
    # All free assets
    return(list(lambda = as.numeric((c4[i]-c1*bi)/c),bi = bi))
  }else{
    onesB <- matrix(1,length(wB))
    l1 <- t(onesB) %*% wB
    l2 <- covarF_inv %*% covarFB
    l3 <- l2 %*% wB
    l2 <- t(onesF) %*% l3
    return(list(lambda = as.numeric(((1-l1+l2)*c4[i]-c1*(bi+l3[i]))/c),bi = bi))
  }
}


# getMatrices --- -------------------------------------------------------------
getMatrices1 <- function(mean, covar, solution_set, f){
  # Slice covarF,covarFB,covarB,meanF,meanB,wF,wB
  covarF <- covar[f,f]
  meanF <- mean[f]
  b <- (1:length(mean))[-f]
  covarFB <- covar[f,b]
  wB <- solution_set[,ncol(solution_set)][b] 
  return(list(covarF=covarF,covarFB=covarFB,meanF=meanF,wB=wB))
}

cla.solver.remove1 <- function(mean, covar,lB, uB){                                             
  # Compute the turning points,free sets and weights
  f <- initAlgo1(mean, lB, uB)$index
  w <- as.matrix(initAlgo1(mean, lB, uB)$weights)
  solution_set <- w # store solution
  lambdas <-  NA_integer_
  gammas <- NA_integer_
  free_weights <- list(f) ### 
  
  while ( TRUE ) {
    # 1) case a): Bound one free weight 
    l_in <- 0 
    
    if (length(f) > 1) {
      
      temp_getMat1 <- getMatrices1(mean, covar,solution_set, f)
      covarF <- temp_getMat1$covarF
      covarFB <- temp_getMat1$covarFB
      meanF <- temp_getMat1$meanF
      wB <- temp_getMat1$wB
      covarF_inv <- chol2inv(chol(covarF))
      
      j <- 1
      for (i in f) {
        temp_compLam <- computeLambda1(covarF_inv,covarFB,meanF,wB,j,as.list(c(lB[i],uB[i])))
        l <- temp_compLam$lambda
        bi <- temp_compLam$bi
        if(l>l_in){  
          l_in <- l
          i_in <- i
          bi_in <- bi
        }
        j <- j + 1
      }
    }
    
    # 2) case b): Free one bounded weight
    l_out <- 0
    if(length(f) < length(mean)) {
      b <- getB(mean, f)
      for(i in b){
        temp_getMat1 <- getMatrices1(mean, covar,solution_set, c(f,i))
        covarF <- temp_getMat1$covarF
        covarFB <- temp_getMat1$covarFB
        meanF <- temp_getMat1$meanF
        #meanF <- c(meanF,mean[i])
        #stopifnot(all.equal(meanF,meanF1))
        wB <- temp_getMat1$wB
        covarF_inv <- chol2inv(chol(covarF))
        
        temp_compLam <- computeLambda1(covarF_inv,covarFB,meanF,wB,length(meanF),solution_set[,ncol(solution_set)][[i]])
        l <- temp_compLam$lambda 
        bi <- temp_compLam$bi
        
        if ((is.na(lambdas[length(lambdas)]) || l<lambdas[length(lambdas)]) && (l>l_out)){
          l_out <- l
          i_out <- i
        }
      }
    }
    
    if( (l_in == 0 || l_in < 0) & (l_out == 0 || l_out < 0) ) { # "stop" when at the min var solution! 
      
      # 3) compute minimum variance solution
      lambdas <- c(lambdas,0)
      temp_getMat1 <- getMatrices1(mean, covar,solution_set, f)
      covarF <- temp_getMat1$covarF
      covarFB <- temp_getMat1$covarFB
      meanF <- temp_getMat1$meanF
      wB <- temp_getMat1$wB
      
      covarF_inv <- chol2inv(chol(covarF))
      meanF <- matrix(0,length(meanF))
    } else {
      
      # 4) decide lambda
      if(l_in>l_out){ 
        lambdas <- c(lambdas,l_in)
        f <- f[f != i_in] # remove i_in within f
        w[i_in]=bi_in # set value at the correct boundary
      }else{
        lambdas <- c(lambdas,l_out)
        f <- c(f,i_out) # append i_out into f
      }
      temp_getMat1 <- getMatrices1(mean, covar,solution_set, f)
      covarF <- temp_getMat1$covarF
      covarFB <- temp_getMat1$covarFB
      meanF <- temp_getMat1$meanF
      wB <- temp_getMat1$wB
      covarF_inv <- chol2inv(chol(covarF))
    }
    
    # 5) compute solution vector
    temp_compW <- computeW1(lambdas,covarF_inv,covarFB,meanF,wB)
    wF <- temp_compW$wF
    g <- temp_compW$gamma
    
    for(i in seq(length(f))){
      w[f[i]]=wF[i]
    }
    
    solution_set <- cbind(solution_set, w) # store solution
    gammas <- c(gammas, g)
    free_weights <- c(free_weights, list(f))
    if(lambdas[length(lambdas)] == 0) {
      break 
    }
  } #end While
  
  
  # temp_pNE <- purgeNummErr(free_weights, gammas, lambdas, lB, uB, solution_set, 1e-09) 
  #solution_set <- temp_pNE[[1]]
  #lambdas <- temp_pNE[[2]]
  #gammas <- temp_pNE[[3]]
  #free_weights <- temp_pNE[[4]]
  
  #temp_purgeEx <- purgeExcess(mean, lambdas, gammas, free_weights, solution_set)
  #solution_set <- temp_purgeEx[[1]]
  #lambdas <- temp_purgeEx[[2]]
  #gammas <- temp_purgeEx[[3]]
  #free_weights <- temp_purgeEx[[4]]
  
  list(mean = mean, covar = covar, solution_set = solution_set, 
       free_weights=free_weights, gammas = gammas, lambdas = lambdas)
  
}

## Part 2: Previous cla.solver.remove (only remove the last ten lines from original cla.solver)
cla.solver.remove <- function(mean, covar,lB, uB)
{
  #initialize data.frames and vectors
  lambdas <- numeric(0) # lambdas 
  gammas <- numeric(0) # gammas
  free_weights <- data.frame() # free weights
  solution_set <- data.frame() # solution
  
  # Compute the turning points,free sets and weights
  temp_initAlgo <- initAlgo(mean, lB, uB)
  f <- temp_initAlgo[[1]]
  w <- temp_initAlgo[[2]]
  solution_set <- rbind(solution_set, w) # store solution
  lambdas <- append(lambdas, NA_integer_) 
  gammas <- append(gammas, NA_integer_) 
  free_weights <- as.data.frame(f)
  
  
  while ( TRUE ) {
    
    # 1) case a): Bound one free weight 
    l_in <- 0 
    if (length(f) > 1) {
      temp_getMat <- getMatrices(mean, covar,solution_set, f)
      covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; meanF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
      covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
      j <- 1
      for (i in f) {
        temp_compLam <- computeLambda(covarF_inv,covarFB,meanF,wB,j,as.list(c(lB[i],uB[i])))
        l <- temp_compLam[[1]]
        bi <- temp_compLam[[2]]
        if(l>l_in){  
          l_in <- l
          i_in <- i
          bi_in <- bi
        }
        j <- j + 1
      }
    }
    
    # 2) case b): Free one bounded weight
    l_out <- 0
    if(length(f) < length(mean)) {
      b <- getB(mean, f)
      for(i in b){
        temp_getMat <- getMatrices(mean, covar,solution_set,c(f,i))
        covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; meanF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
        covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
        temp_compLam <- computeLambda(covarF_inv,covarFB,meanF,wB,length(meanF),solution_set[,ncol(solution_set)][[i]])
        l <- temp_compLam[[1]]; bi <- temp_compLam[[2]]
        
        if ((is.na(lambdas[length(lambdas)]) || l<lambdas[length(lambdas)]) && (l>l_out)){
          l_out <- l
          i_out <- i
        }
      }
    }
    
    if( (l_in == 0 || l_in < 0) & (l_out == 0 || l_out < 0) ) { # "stop" when at the min var solution! 
      
      # 3) compute minimum variance solution
      lambdas <- append(lambdas,0)
      temp_getMat <- getMatrices(mean, covar,solution_set,f)
      covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; meanF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
      covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
      meanF <- matrix(0,length(meanF))
    } else {
      
      # 4) decide lambda
      if(l_in>l_out){ 
        lambdas <- append(lambdas,l_in)
        f <- f[f != i_in] # remove i_in within f
        w[i_in]=bi_in # set value at the correct boundary
      }else{
        lambdas <- append(lambdas,l_out)
        f <- append(f,i_out) # append i_out into f
      }
      temp_getMat <- getMatrices(mean, covar,solution_set, f)
      covarF <- temp_getMat[[1]]
      covarFB <- temp_getMat[[2]]
      meanF <- temp_getMat[[3]]
      wB <- temp_getMat[[4]]
      covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
    }
    
    # 5) compute solution vector
    temp_compW <- computeW(lambdas,covarF_inv,covarFB,meanF,wB)
    wF <- temp_compW[[1]]
    g <- temp_compW[[2]]
    for(i in seq(length(f))){
      w[f[i]]=wF[i]
    }
    solution_set <- cbind(solution_set, w) # store solution
    gammas <- append(gammas, g)
    free_weights <- append(free_weights, as.data.frame(f))
    if(lambdas[length(lambdas)] == 0) {
      break 
    }
  } #end While
  
  
  # temp_pNE <- purgeNummErr(free_weights, gammas, lambdas, lB, uB, solution_set, 1e-09) 
  #solution_set <- temp_pNE[[1]]
  #lambdas <- temp_pNE[[2]]
  #gammas <- temp_pNE[[3]]
  #free_weights <- temp_pNE[[4]]
  
  #temp_purgeEx <- purgeExcess(mean, lambdas, gammas, free_weights, solution_set)
  #solution_set <- temp_purgeEx[[1]]
  #lambdas <- temp_purgeEx[[2]]
  #gammas <- temp_purgeEx[[3]]
  #free_weights <- temp_purgeEx[[4]]
  
  ans <- list(mean, covar, solution_set, free_weights, gammas, lambdas)
  return(ans)
}

## Part 3, compare cla.solver.remove and cla.solver.remove1
#prepare the data, load original data from your own working directory!
company_assets <- read.table(d.file('03-02_CLA_Data_Tot.csv'), header = TRUE, sep = ',')

transData <- function(d, n = 1457){
  assets <- as.matrix(t(d[1:3, 1:n]))
  colnames(assets) <- c("mu", "lB", "uB")
  covar <- as.matrix(d[4:(n+3), 1:n ])
  rownames(covar) <- colnames(covar)
  list(assets = assets, covar = covar)
}

trans_company_assets <- transData(company_assets)

# for function input
getData <- function(d){
  row.names(d$assets) <- NULL
  mu <- as.matrix(d$assets[, "mu"])
  lB <- as.matrix(d$assets[, "lB"])
  uB <- as.matrix(d$assets[, "uB"])
  colnames(d$covar) <- NULL
  rownames(d$covar) <- NULL
  covar <- d$covar
  list(mu = mu, lB = lB, uB = uB, covar = covar)
}
assets_50 <- getData(transData(company_assets, 50))

cla.solver.remove1(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB)
cla.solver.remove(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB)

all.equal(cla.solver.remove(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB), 
          cla.solver.remove1(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB),
          check.attributes = F, check.names = F)

#Compare the time
#install.packages("microbenchmark")
require("microbenchmark")
microbenchmark( cla.solver.remove(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB),
                cla.solver.remove1(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB), 
                times = 3)


## Part 4, validating purgeNummErr 
re <- cla.solver.remove1(assets_50$mu, assets_50$covar, assets_50$lB, assets_50$uB)
dim(re[[3]]) # input, solution set with 45 columns 
pu <- purgeNummErr(re[[4]], re[[5]], re[[6]], assets_50$lB, assets_50$uB, 
             as.data.frame(re[[3]]), 1e-9)

dim(pu[[1]]) # output, solution set with 40 columns 
# could find out column numbers purged: 8 9 16 20 39

# re[[3]] is the input solution set
# purgeNummErr should remove columns which: sum of weights unequal to 1, or 
# any weights out of the bounds, with given tolerance.

colSums(re[[3]]) # take the sum of weights in each column

which(abs(colSums(re[[3]])-1) > 1e-9) # columns that sum unequal to 1
# which is different from 8 9 16 20 39 





