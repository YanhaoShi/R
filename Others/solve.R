funEnv <- function(...) {  ## by Martin Maechler
  e <- list2env(list(...))
  for(n in names(e)) environment(e[[n]]) <- e
  e
}
GetIndex <- function(ind, total_assets){
  list(mu = total_assets$mu[ind, , drop = FALSE], 
       lB = total_assets$lB[ind, , drop = FALSE], 
       uB = total_assets$uB[ind, , drop = FALSE], 
       covar = total_assets$covar[ind, ind])
}
Env2 <- funEnv( 
  info = "purgeChull, colSums != 1", 
  # Initialize the weight -- Find first free weight
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
  
  # computeW -----------------------------------------------------------------
  # compute gamma #(Bailey and Lopez de Prado 2013, 11)
  
  computeW = function(lam,covarFB,muF,wB, covarF){
    #1) compute gamma
    inv <- solve(covarF, cbind(1, muF, covarFB %*% wB, deparse.level = 0))
    g1 <- sum(inv[,2])  
    g2 <- sum(inv[,1])
    if(is.null(wB)){
      g <- as.numeric(-lam*g1/g2+1/g2) 
      w1 <- 0
    } else {
      g3 <- sum(wB)
      w1 <- inv[,3]
      g4 <- sum(w1)
      g <- as.numeric(-lam*g1/g2+(1-g3+g4)/g2)
    }
    #2) compute weights
    w2 <- inv[,1]
    w3 <- inv[,2]
    list(wF = -w1+g*w2+lam*w3, gamma = g)
  },
  
  # computeLambda --------------------------------------------------------------
  # (Niedermayer and Niedermayer 2007, 10)
  # bound to free
  computeLambda = function(covarFB, muF, wB, i, bi, covarF){
    #1) C
    inv <- solve(covarF, cbind(1, muF, covarFB %*% wB))
    c1 <- sum(inv[,1])
    c2i <-  inv[i, 2]
    c3 <- sum(inv[,2])
    c4i <- inv[i, 1]
    Ci <- -c1*c2i + c3*c4i  
    if(Ci == 0){
      return(list(lambda = 0,bi = 0))
    }
    #2) bi
    
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
  
  # free to bound
  #return a list of lambda and b with length(b)
  lcomputeLambda = function(covarFB,muF,wB,i,bi,covarF){
    #1) C
    lambda <- c()
    bi.output <- c()
    inv <- solve(covarF, cbind(1, muF, covarFB %*% wB))
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
      for (tt in seq.int(Ci))
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
  getMatrices = function(mu, covar, w, f){
    # Slice covarF,covarFB,covarB,muF,muB,wF,wB
    covarF <- covar[f,f]
    muF <- mu[f]
    b <- (1:length(mu))[-f]
    covarFB <- covar[f,b]
    wB <- w[b]
    return(list(covarF=covarF,covarFB=covarFB,muF=muF,wB=wB))
  },
  
  # purgeNumErr ---  -------------------------------------------------------------
  # Purge violations: remove w if sum(w)!=1, w<lB, w>uB
  purgeNumErr = function(lB, uB, weights_set, 
                         tol.s, tol.b){ 
    n <- ncol(weights_set) 
    index.s <- abs(colSums(weights_set)-1) < tol.s
    m_lB <- matrix(rep(lB, n), ncol = n)
    m_uB <- matrix(rep(uB, n), ncol = n) 
    m_logical <- (weights_set - m_lB > -tol.b) | (weights_set - m_uB < tol.b)
    index.b <- apply(m_logical, 2, any)
    index.s & index.b
  }, 
  
  ## convex hull
  ## purgeChull
  purgeChull = function(weights_set, mu, covar){
    Sig2 <- colSums(weights_set *(covar %*% weights_set) ) # ms!
    Mu <- t(weights_set) %*% mu
    ch <- chull(cbind(Sig2, Mu))
    # start from largest mu, end at smallest sig2
    range <- c(which.max(Mu[ch]), which.min(Sig2[ch]))
    ch[max(range): min(range)]
  },
  
  cla.solve = function(mu, covar,lB, uB){                                             
    # Compute the turning points,free sets and weights
    f <- initAlgo(mu, lB, uB)$index
    w <- as.matrix(initAlgo(mu, lB, uB)$weights)
    weights_set <- w # store solution
    lambdas <-  c()
    gammas <- NA#########
    free_indices <- list(f) ### 
    list_temp <- list() 
    
    while ( TRUE ) {
      
      
      # 1) case a): Bound one free weight 
      l_in <- 0
      
      
      if(length(f) > 1){  
        get1 <- getMatrices(mu, covar, w, f)
        muF <- get1$muF
        compl <- lcomputeLambda(covarFB = get1$covarFB, muF = get1$muF,
                                wB = get1$wB, i = seq.int(f),
                                bi = list(lB = lB[f],uB = uB[f]), covarF = get1$covarF)
        lam_in<- compl$lambda
        bi <- compl$bi
        l_in <- max(lam_in) 
        
        if(l_in > 0){
          k <- which.max(lam_in)
          i_in <- f[k]
          bi_in <- bi[k]
        } 
      }  
      
      # 2) case b): Free one bounded weight
      
      l_out <- 0
      
      if(length(f) < length(mu)) {
        
        b <- seq_along(mu)[-f]
        
        lam_out <- sapply(b, function(bi) {
          get_i <- getMatrices(mu, covar, w, c(f,bi))
          muF <- get_i$muF
          computeLambda(covarFB = get_i$covarFB, muF = get_i$muF,
                        wB = get_i$wB, i = length(get_i$muF),
                        bi = w[bi], covarF = get_i$covarF)$lambda
        })
        
        if (!is.null(lambdas) && any(sml <- lam_out < lambdas[length(lambdas)])) {
          lam_out <- lam_out[sml]
          b       <- b      [sml]
        }
        
        i_out <- which.max(lam_out)
        l_out <- lam_out[i_out]
        i_out <- b      [i_out] # one only !
        get_out <- getMatrices(mu, covar, w, c(f,i_out))
        
        
      }
      
      if( (l_in <= 0) & (l_out <= 0) ) { # "stop" when at the min var solution! 
        
        # 3) compute minimum variance solution
        lambdas <- c(lambdas,0)
        muF <- matrix(0,length(muF))
        temp_compW <- computeW(lam, covarFB = get1$covarFB, muF ,wB = get1$wB,
                               covarF = get1$covarF)
        wF <- temp_compW$wF
        g <- temp_compW$gamma
        
        list_temp <- c(list_temp, list(get1$covarF))
        
        for(i in seq.int(f)){
          w[f[i]]=wF[i]
        }
        
        #  w[f[seq.int(f)]] <- wF[seq.int(f)]
        
        weights_set <- cbind(weights_set, w) # store solution
        gammas <- c(gammas, g)
        break
        
      } else {
        
        # 4) decide lambda       
        if(l_in>l_out){ 
          lambdas <- c(lambdas,l_in)
          f <- f[f != i_in] # remove i_in within f
          get_in <- getMatrices(mu, covar,w, f) # should move to next line!
          w[i_in] <- bi_in # set value at the correct boundary
          
          muF <- get_in$muF
          
          temp_compW <- computeW(l_in, covarFB = get_in$covarFB,
                                 muF = get_in$muF, wB = get_in$wB, 
                                 covarF = get_in$covarF)
          wF <- temp_compW$wF
          g <- temp_compW$gamma 
          list_temp <- c(list_temp, list(get_in$covarF))
          
        }else{
          lambdas <- c(lambdas,l_out)
          f <- c(f,i_out) # add i_out into f
          muF <- get_out$muF
          temp_compW <- computeW(l_out, covarFB = get_out$covarFB, 
                                 muF = get_out$muF, wB = get_out$wB,
                                 covarF = get_out$covarF)
          wF <- temp_compW$wF
          g <- temp_compW$gamma
          list_temp <- c(list_temp, list(get_out$covarF))
        }
        
        
        for(i in seq.int(f)){
          w[f[i]]=wF[i]
        }
        # w[f] <- as.matrix(wF)[seq.int(f)]
        weights_set <- cbind(weights_set, w) # store solution
        gammas <- c(gammas, g)
        free_indices <- c(free_indices, list(f))
        lam <- lambdas[length(lambdas)]
      }
      
    } #end While  
    lambdas <- c(NA, lambdas) # The first step has no lambda or gamma, add NA instead.
    
    free_indices <- c(free_indices, NA) # The last step stop without weight set, add NA.
    
    # purge 
    i <- purgeNumErr(lB, uB, weights_set, 1e-9, 1e-10) 
    weights_set_purge <- weights_set[,i]
    lambdas_purge <- lambdas[i]
    gammas_purge <- gammas[i]
    free_indices_purge <- free_indices[i]
    
    k <- purgeChull(weights_set_purge, mu, covar)
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

testInv <- function(result, mu, covar){
  free_indices <- result$free_indices
  weights_set <- result$weights_set
  n <- ncol(result$weights_set)
  inv_covarF <- function(j, inv){
    w <- weights_set[, j]
    f <- free_indices[[j]]
    get <- Env2$getMatrices(mu, covar, w, f)
    covarF <- get$covarF
    covarFB <- get$covarFB
    muF <- get$muF
    wB <- get$wB
    if(inv == 1) chol2inv(chol(covarF)) %*% cbind(1, muF, covarFB %*% wB)
    else if(inv == 2) solve(covarF) %*% cbind(1, muF, covarFB %*% wB)
    else if(inv == 3) solve(covarF, cbind(1, muF, covarFB %*% wB))
  }
  micro <- microbenchmark(inv1 = inv1 <- sapply(1:(n-1), function(j) inv_covarF(j, 1)),
                          inv2 = inv2 <- sapply(1:(n-1), function(j) inv_covarF(j, 2)),
                          inv3 = inv3 <- sapply(1:(n-1), function(j) inv_covarF(j, 3)), 
                          times = 10)
  is.allequal <- all.equal(inv1, inv3) & all.equal(inv2, inv3)
  list(micro = micro, is.allequal = is.allequal)
}

solve_time <- list()
for(n in seq(floor(1457/50))){
  assets <- GetIndex(seq(n * 50), assets1457)
  result <- Env2$cla.solve(assets$mu, assets$covar, assets$lB, assets$uB)
  solve_time[[n]] <- testInv(result, mu = assets$mu, covar = assets$covar)$micro
}





