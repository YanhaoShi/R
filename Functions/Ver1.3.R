Env1.3 <- funEnv( 
  info = "lcomputelambda, solve()",
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
  
  # getB1 ---  --------------------------------------------------------------------
  getB = function(mu,f) (1:length(mu))[-f],
  
  # reduceMatrix ---  ------------------------------------------------------------
  reduceMatrix =  function(input, X, Y){
    if  (is.vector(input)){
      return(input[X])
    }
    else { return(input[X,Y,drop=F])}     
  },
  
  
  # getMinVar ---  ---------------------------------------------------------------
  # Get the minimum variance solution
  getMinVar = function(covar,solution_set1){          ## getMinStd !
    n<-ncol(solution_set1)
    var <- rep(1,n)
    for(i in 1:n){
      w<-solution_set1[,i]
      var[i] <- t(w) %*% covar %*% w
    }
    
    list(index_minvar= sqrt(min(var)), weights_minvar=solution_set1[ ,which.min(var)] )
  }, 
  
  # efFrontier -----------------------------------------------------------------
  # Get the efficient frontier
  efFrontier = function(mu, covar, solution_set, dataPoints){ 
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
    mu <- t(mu) %*% weights
    
    sigma <- apply(weights, 2, function(x) sqrt(t(x) %*% covar %*% x))
    
    return(list(mu = mu, weights = weights, sigma = sigma))
  },  
  
  # computeW -----------------------------------------------------------------
  # compute gamma #(Bailey and Lopez de Prado 2013, 11)
  
  computeW = function(lambdas,covarFB,muF,wB, covarF){
    #1) compute gamma
    onesF <- rep(1, length(muF))
    inv <- solve(covarF, cbind(onesF, muF, covarFB %*% wB))
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
  },

  # computeLambda --------------------------------------------------------------
  # (Niedermayer and Niedermayer 2007, 10)
  computeLambda = function(covarFB,muF,wB,i,bi,covarF){
    #1) C
    onesF <- muF
    onesF[] <- 1
    inv <- solve(covarF, cbind(onesF, muF, covarFB %*% wB))
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
  },
  
  #return a list of lambda and b with length(b)
  lcomputeLambda = function(covarFB,muF,wB,i,bi,covarF){
    #1) C
    lambda <- c()
    bi.output <- c()
    onesF <- muF
    onesF[] <- 1
    inv <- solve(covarF, cbind(onesF, muF, covarFB %*% wB))
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
  },
  # getMatrices --- -------------------------------------------------------------
  getMatrices = function(mu, covar, solution_set, f){
    # Slice covarF,covarFB,covarB,muF,muB,wF,wB
    covarF <- covar[f,f]
    muF <- mu[f]
    b <- (1:length(mu))[-f]
    covarFB <- covar[f,b]
    wB <- solution_set[,ncol(solution_set)][b] 
    return(list(covarF=covarF,covarFB=covarFB,muF=muF,wB=wB))
  },
  
  # purgeNumErr ---  -------------------------------------------------------------
  # Purge violations: remove w if sum(w)!=1, w<lB, w>uB
  purgeNummErr = function(lB, uB, solution_set, 
                           tol.s, tol.b){ 
    n <- ncol(solution_set) 
    index.s <- abs(colSums(solution_set)-1) < tol.s
    m_lB <- matrix(rep(lB, n), ncol = n)
    m_uB <- matrix(rep(uB, n), ncol = n) 
    m_logical <- (solution_set - m_lB > -tol.b) | (solution_set - m_uB < tol.b)
    index.b <- apply(m_logical, 2, any)
    index.s & index.b
  }, 
  
  # purgeExcess ---   -------------------------------------------------------------
  # Remove violations of the convex hull, ensure mu2>mu3, ...
  purgeExcess = function(mu, solution_set, tol){ 
    n <- ncol(solution_set)
    Mu <- t(mu) %*% solution_set
    c(which( (Mu[2:(n-1)] - Mu[3:n]) < -tol ), n+1 )  
  },
  
  
  purgeExcess2 = function(mu, solution_set, tol){ 
    n <- ncol(solution_set)
    Mu <- t(mu) %*% solution_set #if mu_n+1 > mu_n, remove mu_n+1
    index.ex <- c(TRUE, TRUE, ((Mu[2:(n-1)] - Mu[3:n]) > -tol ) )
  },
  
  cla.solve = function(cla.input){   
    mu <- cla.input$mu
    covar <- cla.input$covar
    lB <- cla.input$lB
    uB <- cla.input$uB
    # Compute the turning points,free sets and weights
    f <- initAlgo(mu, lB, uB)$index
    w <- as.matrix(initAlgo(mu, lB, uB)$weights)
    weights_set <- w # store solution
    lambdas <-  c()
    gammas <- c()
    free_indices <- list(f) ### 
    list_temp <- list() 
    while ( TRUE ) {
      # 1) case a): Bound one free weight 
      l_in <- 0
      temp_getMat1 <- getMatrices(mu, covar,weights_set, f)
      covarF1 <- temp_getMat1$covarF
      covarFB1 <- temp_getMat1$covarFB
      muF1 <- temp_getMat1$muF
      wB1 <- temp_getMat1$wB
      
      if(length(f) > 1){  
        compl <- lcomputeLambda(covarFB1,muF1,wB1,seq(length(f)),
                                list(lB = lB[f],uB = uB[f]),covarF1)
        lam <- compl$lambda
        bi <- compl$bi
        i_in <- ifelse(max(lam)>0, f[which.max(lam)], 0)
        l_in <- max(lam, 0) 
        bi_in <- ifelse(max(lam)>0, bi[which.max(lam)], 0)
      }  
      
      # 2) case b): Free one bounded weight
      l_out <- 0
      if(length(f) < length(mu)) {
        
        b <- (1:length(mu))[-f]
        for(i in b){
          
          temp_getMat2 <- getMatrices(mu, covar, weights_set, c(f,i))
          covarF2 <- temp_getMat2$covarF
          covarFB2 <- temp_getMat2$covarFB
          muF2 <- temp_getMat2$muF
          wB2 <- temp_getMat2$wB
          temp_compLam <- computeLambda(covarFB2,muF2,wB2,length(muF2),
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
        muF <- matrix(0,length(muF))
        temp_compW <- computeW(lambdas,covarFB1,muF,wB1,covarF1)
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
          temp_getMat1 <- getMatrices(mu, covar,weights_set, f)
          covarF <- temp_getMat1$covarF
          covarFB <- temp_getMat1$covarFB
          muF <- temp_getMat1$muF
          wB <- temp_getMat1$wB
          temp_compW <- computeW(lambdas,covarFB,muF,wB,covarF)
          wF <- temp_compW$wF
          g <- temp_compW$gamma 
          list_temp <- c(list_temp, list(covarF))
          
        }else{
          lambdas <- c(lambdas,l_out)
          f <- c(f,i_out) # append i_out into f
          
          temp_getMat1 <- temp_getMat_out
          covarF <- temp_getMat1$covarF
          covarFB <- temp_getMat1$covarFB
          muF <- temp_getMat1$muF
          wB <- temp_getMat1$wB
          temp_compW <- computeW(lambdas,covarFB,muF,wB,covarF)
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
    
    k <- purgeExcess(mu, weights_set_purge, 1e-12)
    weights_set_purge <- weights_set_purge[,-k] 
    lambdas_purge <- lambdas_purge[-k] 
    gammas_purge <- gammas_purge[-k]
    free_indices_purge <- free_indices_purge[-k] 
    
    list(weights_set_purge = weights_set_purge, free_indices_purge = free_indices_purge, 
         gammas_purge = gammas_purge, lambdas_purge = lambdas_purge,
         weights_set = weights_set, free_indices=free_indices, 
         gammas = gammas, lambdas = lambdas, covarF = list_temp)
    
  }
)