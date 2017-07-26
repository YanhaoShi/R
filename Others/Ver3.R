Method3 <-list( 
  # Initialize the weight -- Find first free weight
  initAlgo = (initAlgo = function(mu, lB, uB ){
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
  }), 
  
  # getB1 ---  --------------------------------------------------------------------
  getB = (getB = function(mean,f) (1:length(mean))[-f]),
  
  # reduceMatrix ---  ------------------------------------------------------------
  reduceMatrix = (reduceMatrix = function(input, X, Y){
    if  (is.vector(input)){
      return(input[X])
    }
    else { return(input[X,Y,drop=F])}     
  }),
  
  
  # getMinVar ---  ---------------------------------------------------------------
  # Get the minimum variance solution
  getMinVar = (getMinVar = function(covar,solution_set1){          ## getMinStd !
    n<-ncol(solution_set1)
    var <- rep(1,n)
    for(i in 1:n){
      w<-solution_set1[,i]
      var[i] <- t(w) %*% covar %*% w
    }
    
    list(index_minvar= sqrt(min(var)), weights_minvar=solution_set1[ ,which.min(var)] )
  }), 
  
  # efFrontier -----------------------------------------------------------------
  # Get the efficient frontier
  efFrontier = (efFrontier = function(mean, covar, solution_set, dataPoints){ 
    n<-ncol(solution_set)
    mu <- c()
    sigma <- c()
    weights <- c()
    a <- seq(0,1,length=dataPoints+1)[-1]  #dataPoints-2 in between
    
    for(i in 1:(n-1)){
      
      w0 <- solution_set[,i]
      w1 <- solution_set[,i+1]
      
      w <- w0 %*% t(a) + w1 %*% t(1-a)   #weights between
      weights <- cbind(weights,w)
    }
    weights <- cbind(weights, solution_set[,n])
    mu <- t(mean) %*% weights
    
    sigma <- apply(weights, 2, function(x) sqrt(t(x) %*% covar %*% x))
    
    return(list(mu = mu, weights = weights, sigma = sigma))
  }),  
  
  # computeW -----------------------------------------------------------------
  # compute gamma #(Bailey and Lopez de Prado 2013, 11)
  
  computeW = (computeW = function(lambdas,covarFB,meanF,wB, covarF){
    #1) compute gamma
    onesF <- rep(1, length(meanF))
    inv <- solve(covarF, cbind(onesF, meanF, covarFB %*% wB))
    g1 <- sum(inv[,2])  
    g2 <- sum(inv[,1])
    l.tail <- lambdas[length(lambdas)]
    if(is.null(wB)){
      g <- as.numeric(-l.tail*g1/g2+1/g2) 
      w1 <- 0
    } else {
      g3 <- sum(wB)
      w1 <- inv[,3]
      g4 <- sum(w1)
      g <- as.numeric(-l.tail*g1/g2+(1-g3+g4)/g2)
    }
    #2) compute weights
    w2 <- inv[,1]
    w3 <- inv[,2]
    list(wF = -w1+g*w2+l.tail*w3, gamma = g)
  }),
  
  # computeLambda --------------------------------------------------------------
  # (Niedermayer and Niedermayer 2007, 10)
  computeLambda = (computeLambda = function(covarFB,meanF,wB,i,bi,covarF){
    #1) C
    onesF <- meanF
    onesF[] <- 1
    inv <- solve(covarF, cbind(onesF, meanF, covarFB %*% wB))
    c1 <- sum(inv[,1])
    c2i <-  inv[i, 2]
    c3 <- sum(inv[,2])
    c4i <- inv[i, 1]
    Ci <- -c1*c2i + c3*c4i  
    if(Ci == 0){
      return(list(lambda = 0,bi = 0))
    }
    #2) bi
    if(is.list(bi)){
      
      bi <- ifelse(Ci > 0, bi[[2]], bi[[1]])
    }
    #3) lambda
    if(length(wB)==0){
      # All free assets
      return(list(lambda = (c4i - c1 * bi)/Ci,bi = bi))
    }else{
      l1 <- sum(wB)
      l3 <- inv[,3]
      l2 <- sum(l3)
      return(list(lambda = ((1-l1+l2)*c4i-c1*(bi+inv[i,3]))/Ci, bi = bi))
    }
  }),
  
  #return a list of lambda and b with length(b)
  lcomputeLambda = (lcomputeLambda = function(covarFB,meanF,wB,i,bi,covarF){
    #1) C
    lambda <- c()
    bi.output <- c()
    onesF <- meanF
    onesF[] <- 1
    inv <- solve(covarF, cbind(onesF, meanF, covarFB %*% wB))
    c1 <- sum(inv[,1])
    c2i <-  inv[i, 2]
    c3 <- sum(inv[,2])
    c4i <- inv[i, 1]
    Ci <- -c1*c2i + c3*c4i  
    
    k <- which(Ci == 0)
    lambda[k] <- 0
    bi.output[k] <- 0
    
    #2) bi
    if(is.list(bi)){
      for (tt in seq(length(Ci)))
        bi.output[tt] <- ifelse(Ci[tt] > 0, bi$uB[tt], bi$lB[tt])
    }
    #3) lambda
    if(length(wB)==0){
      # All free assets
      return(list(lambda = (c4i - c1 * bi.output)/Ci,bi = bi.output))
    }else{
      l1 <- sum(wB)
      l3 <- inv[,3]
      l2 <- sum(l3)
      return(list(lambda = ((1-l1+l2)*c4i-c1*(bi.output+inv[i,3]))/Ci, bi = bi.output))
    }
  }),
  # getMatrices --- -------------------------------------------------------------
  getMatrices = (getMatrices = function(mean, covar, solution_set, f){
    # Slice covarF,covarFB,covarB,meanF,meanB,wF,wB
    covarF <- covar[f,f]
    meanF <- mean[f]
    b <- (1:length(mean))[-f]
    covarFB <- covar[f,b]
    wB <- solution_set[,ncol(solution_set)][b] 
    return(list(covarF=covarF,covarFB=covarFB,meanF=meanF,wB=wB))
  }),
  
  # purgeNumErr ---  -------------------------------------------------------------
  # Purge violations: remove w if sum(w)!=1, w<lB, w>uB
  purgeNummErr = (purgeNummErr = function(lB, uB, solution_set, 
                           tol.s, tol.b){ 
    n <- ncol(solution_set) 
    index.s <- abs(colSums(solution_set)-1) < tol.s
    m_lB <- matrix(rep(lB, n), ncol = n)
    m_uB <- matrix(rep(uB, n), ncol = n) 
    m_logical <- (solution_set - m_lB > -tol.b) | (solution_set - m_uB < tol.b)
    index.b <- apply(m_logical, 2, any)
    index.s & index.b
  }), 
  
  ## convex hull
  purgeChull = (purgeChull <- function(solution_set, mean, covar){
    Sig2 <- colSums(solution_set *(covar %*% solution_set) )
    Mu <- t(solution_set) %*% mean
    ch <- sort(chull(cbind(Sig2, Mu)))
    index.chull <- ch[Sig2[ch] <= Sig2[which.max(Mu)]]
    index.chull
  }),
  
  cla.solver = function(mean, covar,lB, uB){                                             
    # Compute the turning points,free sets and weights
    f <- initAlgo(mean, lB, uB)$index
    w <- as.matrix(initAlgo(mean, lB, uB)$weights)
    weights_set <- w # store solution
    lambdas <-  c()
    gammas <- c()
    free_indices <- list(f) ### 
    list_temp <- list() 
    while ( TRUE ) {
      # 1) case a): Bound one free weight 
      l_in <- 0
      temp_getMat1 <- getMatrices(mean, covar,weights_set, f)
      covarF1 <- temp_getMat1$covarF
      covarFB1 <- temp_getMat1$covarFB
      meanF1 <- temp_getMat1$meanF
      wB1 <- temp_getMat1$wB
      
      if(length(f) > 1){  
        compl <- lcomputeLambda(covarFB1,meanF1,wB1,seq(length(f)),
                                list(lB = lB[f],uB = uB[f]),covarF1)
        lam <- compl$lambda
        bi <- compl$bi
        i_in <- ifelse(max(lam)>0, f[which.max(lam)], 0)
        l_in <- max(lam, 0) 
        bi_in <- ifelse(max(lam)>0, bi[which.max(lam)], 0)
      }  
      
      # 2) case b): Free one bounded weight
      l_out <- 0
      if(length(f) < length(mean)) {
        
        b <- (1:length(mean))[-f]
        for(i in b){
          
          temp_getMat2 <- getMatrices(mean, covar, weights_set, c(f,i))
          covarF2 <- temp_getMat2$covarF
          covarFB2 <- temp_getMat2$covarFB
          meanF2 <- temp_getMat2$meanF
          wB2 <- temp_getMat2$wB
          temp_compLam <- computeLambda(covarFB2,meanF2,wB2,length(meanF2),
                                         weights_set[i,ncol(weights_set)],covarF2)
          l2 <- temp_compLam$lambda 
          bi <- temp_compLam$bi
          
          if ((length(lambdas)==0 || l2<lambdas[length(lambdas)]) && (l2>l_out)){
            l_out <- l2
            i_out <- i
            temp_getMat_out <- temp_getMat2
          }
        }
      }
      
      if( (l_in <= 0) & (l_out <= 0) ) { # "stop" when at the min var solution! 
        
        # 3) compute minimum variance solution
        lambdas <- c(lambdas,0)
        meanF <- matrix(0,length(meanF))
        temp_compW <- computeW(lambdas,covarFB1,meanF,wB1,covarF1)
        wF <- temp_compW$wF
        g <- temp_compW$gamma
        
        list_temp <- c(list_temp, list(covarF1))
        
        for(i in seq(length(f))){
          w[f[i]]=wF[i]
        }
        weights_set <- cbind(weights_set, w) # store solution
        gammas <- c(gammas, g)
        break
        
      } else {
        
        # 4) decide lambda       
        if(l_in>l_out){ 
          lambdas <- c(lambdas,l_in)
          f <- f[f != i_in] # remove i_in within f
          w[i_in] = bi_in # set value at the correct boundary
          temp_getMat1 <- getMatrices(mean, covar,weights_set, f)
          covarF <- temp_getMat1$covarF
          covarFB <- temp_getMat1$covarFB
          meanF <- temp_getMat1$meanF
          wB <- temp_getMat1$wB
          temp_compW <- computeW(lambdas,covarFB,meanF,wB,covarF)
          wF <- temp_compW$wF
          g <- temp_compW$gamma 
          list_temp <- c(list_temp, list(covarF))
          
        }else{
          lambdas <- c(lambdas,l_out)
          f <- c(f,i_out) # append i_out into f
          
          temp_getMat1 <- temp_getMat_out
          covarF <- temp_getMat1$covarF
          covarFB <- temp_getMat1$covarFB
          meanF <- temp_getMat1$meanF
          wB <- temp_getMat1$wB
          temp_compW <- computeW(lambdas,covarFB,meanF,wB,covarF)
          wF <- temp_compW$wF
          g <- temp_compW$gamma
          list_temp <- c(list_temp, list(covarF))
        }
        
        
        for(i in seq(length(f))){
          w[f[i]]=wF[i]
        }
        
        weights_set <- cbind(weights_set, w) # store solution
        gammas <- c(gammas, g)
        free_indices <- c(free_indices, list(f))
        
      }
      
    } #end While  
    lambdas <- c(NA, lambdas) # The first step has no lambda or gamma, add NA instead.
    gammas <- c(NA, gammas)
    free_indices <- c(free_indices, NA) # The last step stop without weight set, add NA.
    
    # purge 
    i <- purgeNummErr(lB, uB, weights_set, 1e-9, 1e-10) 
    weights_set_purge <- weights_set[,i]
    lambdas_purge <- lambdas[i]
    gammas_purge <- gammas[i]
    free_indices_purge <- free_indices[i]
    
    k <- purgeChull(weights_set_purge, mean, covar)
    weights_set_purge <- weights_set_purge[,k] 
    lambdas_purge <- lambdas_purge[k] 
    gammas_purge <- gammas_purge[k]
    free_indices_purge <- free_indices_purge[k] 

    list(weights_set_purge = weights_set_purge, free_indices_purge = free_indices_purge, 
         gammas_purge = gammas_purge, lambdas_purge = lambdas_purge,
         weights_set = weights_set, free_indices=free_indices, 
         gammas = gammas, lambdas = lambdas, covarF = list_temp)
  }
)