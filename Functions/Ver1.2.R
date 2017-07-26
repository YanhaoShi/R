Env1.2 <- funEnv( 
  info = "chol=T, change purgeNum",
  initAlgo = function(mu, lB, uB ){
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
  }, 
  
  # getB ---  --------------------------------------------------------------------
  getB = function(mu,f) (1:length(mu))[-f],
  
  # reduceMatrix ---  ------------------------------------------------------------
  reduceMatrix = function(input, X, Y){
    if  (is.vector(input)){
      return(input[X])
    }
    else { return(input[X,Y,drop=F])}     
  }, 
  
  # purgeNumErr ---  -------------------------------------------------------------
  # Purge violations: remove w if sum(w)!=1, w<lB, w>uB
  # tolerance, e.g 10e-10
  purgeNumErr = function(free_indices, gammas, lambdas, lB, uB, weights_set, tol=1e-9){ 
    i <- 1
    # k <- which(abs(colSums(weights_set)-1) > tol)
    # which(abs(colSums(re[[3]])-1) > tol)
    k <- c()
    l <- c()
    while(TRUE){# 
      flag <- FALSE
      if(i == ncol(weights_set)){
        return(list(weights_set,lambdas,gammas,free_indices, index = c(k,l)))
        #return(list(k,l))
      }else{##
        if((abs(sum(weights_set[,i]))-1) > tol){ ###22
          flag <- TRUE
          
          k<-c(k,i)
        } else {
          for(j in seq(length(weights_set[,i]))){
            if((weights_set[j,i] - lB[j]) < (-tol) || (weights_set[j,i] - uB[j]) > (tol)){
              flag <- TRUE
              l <- c(l,i)
            }
          }
        } 
        if(flag==TRUE){
          weights_set <- weights_set[,-i] 
          lambdas <- lambdas[-i]
          gammas <- gammas[-i]
          free_indices <- free_indices[-i] 
          
        }else{
          i <- i+1
        }
      }##
    }#
    list(weights_set,lambdas,gammas,free_indices, index = c(k,l))
    
    #return(list(k,l))
    
  },
  
  purgeNumErr.draft = function(free_indices, gammas, lambdas, lB, uB, weights_set1, tol=1e-9){
    # ncol=1时返回原值!!!!!
    k <- which(abs(colSums(weights_set1)-1) > tol)  ## k=1?
    
    for(i in 1:ncol(weights_set1)){
      #if(any(c(weights_set1[,i]-lB < -tol,weights_set1[,i]-uB > tol))){
      k<-c(which(weights_set1[,i]-lB < -tol),which(weights_set1[,i]-uB > tol),k)
      
      
    }
    k<- k[!duplicated(k)]
    weights_set1 <- weights_set1[,-k] 
    lambdas <- lambdas[-k]
    gammas <- gammas[-k]
    free_indices <- free_indices[-k] 
    
    return(list(weights_set1,lambda=lambdas,lgamma=gammas,weights_free=free_indices))
  },
  
  # purgeExcess ---   -------------------------------------------------------------
  # Remove violations of the convex hull, ensure mu2>mu3, ...
  purgeExcess = function(mu, lambdas, gammas, free_indices, weights_set){ 
    n <- ncol(weights_set)
    Mu <- t(mu) %*% weights_set
    
    k<-which( Mu[2:(n-1)] < Mu[3:n] )  # 0 or tol?
    
    weights_set <- weights_set[,-k] 
    lambdas <- lambdas[-k] 
    gammas <- gammas[-k]
    free_indices <- free_indices[-k] 
    
    return(list(weights_set,lambdas, gammas, free_indices, k = k))
  },
  
  # getMinVar ---  ---------------------------------------------------------------
  # Get the minimum variance solution
  getMinVar = function(covar,weights_set1){          ## getMinStd !
    n<-ncol(weights_set1)
    var <- rep(1,n)
    for(i in 1:n){
      w<-weights_set1[,i]
      var[i] <- t(w) %*% covar %*% w
    }
    
    list(index_minvar= sqrt(min(var)), weights_minvar=weights_set1[ ,which.min(var)] )
  }, 
  
  # efFrontier -----------------------------------------------------------------
  # Get the efficient frontier
  efFrontier = function(mu, covar, weights_set1, dataPoints){ 
    n<-ncol(weights_set1)
    mu <- c()
    sigma <- c()
    weights <- c()
    
    a <- seq(0,1,length=dataPoints+1)[-1]  #dataPoints-2 in between
    
    for(i in 1:(n-1)){
      w0 <- weights_set[i]
      w1 <- weights_set[i+1]
      
      weights <- w0%*%t(a)+w1 %*% t(1-a)   #weights between
      w <- cbind(weights,w)
      mu <- c(mu,colSums(weights))
      sigma <- c(sigma, sqrt(t(w) %*% covar %*% w))
    }
    return(list(mu, sigma, weights))
  }, 
  
  # computeW -----------------------------------------------------------------
  # compute gamma #(Bailey and Lopez de Prado 2013, 11)
  computeW = function(lambdas,covarF_inv,covarFB,muF,wB){
    #1) compute gamma
    g1 <- sum(covarF_inv %*% muF)  #rrr
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
    w3 <- covarF_inv %*% muF
    
    return(list(wF = -w1+g*w2+lambdas[length(lambdas)]*w3, gamma = g))
  },
  
  
  # computeLambda --------------------------------------------------------------
  # (Niedermayer and Niedermayer 2007, 10)
  computeLambda = function(covarF_inv,covarFB,muF,wB,i,bi){
    #1) C
    onesF <- matrix(1,length(muF)) # (kx1)-vector
    c1 <- (t(onesF) %*% covarF_inv) %*% onesF
    c2i <- covarF_inv[i,] %*% muF
    c3 <- t(onesF) %*% covarF_inv %*% muF
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
  },
  
  
  # getMatrices --- -------------------------------------------------------------
  getMatrices = function(mu, covar, weights_set, f){
    # Slice covarF,covarFB,covarB,muF,muB,wF,wB
    covarF <- covar[f,f]
    muF <- mu[f]
    b <- (1:length(mu))[-f]
    covarFB <- covar[f,b]
    wB <- weights_set[,ncol(weights_set)][b] 
    return(list(covarF=covarF,covarFB=covarFB,muF=muF,wB=wB))
  },
  
  ##
  cla.solve = function(cla.input){        
    mu <- cla.input$mu
    covar <- cla.input$covar
    lB <- cla.input$lB
    uB <- cla.input$uB
    # Compute the turning points,free sets and weights
    f <- initAlgo(mu, lB, uB)$index
    w <- as.matrix(initAlgo(mu, lB, uB)$weights)
    weights_set <- w # store solution
    lambdas <-  NA_integer_
    gammas <- NA_integer_
    free_indices <- list(f) ### 
    
    while ( TRUE ) {
      # 1) case a): Bound one free weight 
      l_in <- 0 
      
      if (length(f) > 1) {
        
        temp_getMat1 <- getMatrices(mu, covar,weights_set, f)
        covarF <- temp_getMat1$covarF
        covarFB <- temp_getMat1$covarFB
        muF <- temp_getMat1$muF
        wB <- temp_getMat1$wB
        covarF_inv <- chol2inv(chol(covarF))
        
        j <- 1
        for (i in f) {
          temp_compLam <- computeLambda(covarF_inv,covarFB,muF,wB,j,as.list(c(lB[i],uB[i])))
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
      if(length(f) < length(mu)) {
        b <- getB(mu, f)
        for(i in b){
          temp_getMat1 <- getMatrices(mu, covar,weights_set, c(f,i))
          covarF <- temp_getMat1$covarF
          covarFB <- temp_getMat1$covarFB
          muF <- temp_getMat1$muF
          #muF <- c(muF,mu[i])
          #stopifnot(all.equal(muF,muF1))
          wB <- temp_getMat1$wB
          covarF_inv <- chol2inv(chol(covarF))
          
          temp_compLam <- computeLambda(covarF_inv,covarFB,muF,wB,length(muF),weights_set[,ncol(weights_set)][[i]])
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
        temp_getMat1 <- getMatrices(mu, covar,weights_set, f)
        covarF <- temp_getMat1$covarF
        covarFB <- temp_getMat1$covarFB
        muF <- temp_getMat1$muF
        wB <- temp_getMat1$wB
        
        covarF_inv <- chol2inv(chol(covarF))
        muF <- matrix(0,length(muF))
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
        temp_getMat1 <- getMatrices(mu, covar,weights_set, f)
        covarF <- temp_getMat1$covarF
        covarFB <- temp_getMat1$covarFB
        muF <- temp_getMat1$muF
        wB <- temp_getMat1$wB
        covarF_inv <- chol2inv(chol(covarF))
      }
      
      # 5) compute solution vector
      temp_compW <- computeW(lambdas,covarF_inv,covarFB,muF,wB)
      wF <- temp_compW$wF
      g <- temp_compW$gamma
      
      for(i in seq(length(f))){
        w[f[i]]=wF[i]
      }
      
      weights_set <- cbind(weights_set, w) # store solution
      gammas <- c(gammas, g)
      free_indices <- c(free_indices, list(f))
      if(lambdas[length(lambdas)] == 0) {
        break 
      }
    } #end While
    
    temp_pNE <- purgeNumErr(free_indices, gammas, lambdas, lB, uB, weights_set, 1e-09) 
    weights_set_purge <- temp_pNE[[1]]
    lambdas_purge <- temp_pNE[[2]]
    gammas_purge <- temp_pNE[[3]]
    free_indices_purge <- temp_pNE[[4]]
    
    temp_purgeEx <- purgeExcess(mu, lambdas_purge, gammas_purge, 
                                free_indices_purge, weights_set_purge)
    weights_set_purge <- temp_purgeEx[[1]]
    lambdas_purge <- temp_purgeEx[[2]]
    gammas_purge <- temp_purgeEx[[3]]
    free_indices_purge <- temp_purgeEx[[4]]
    
    list(weights_set_purge = weights_set_purge, free_indices_purge = free_indices_purge, 
         gammas_purge = gammas_purge, lambdas_purge = lambdas_purge,
         weights_set = weights_set, free_indices=free_indices, 
         gammas = gammas, lambdas = lambdas)
    
  }
)