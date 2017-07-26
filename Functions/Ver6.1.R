Env6.1 <- funEnv( 
  info = "new input/output ", 
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
    i <- index.new[i.new] 
    list(index = i, weights = w)         # return the index of first free asset and vector w
  }, 
  
  # computeW -----------------------------------------------------------------
  # compute gamma #(Bailey and Lopez de Prado 2013, 11)
  computeW = function(lam, get, muF = get$muF){
    covarFB <- get$covarFB
    wB <- get$wB
    covarF <- get$covarF
    
    #1) compute gamma
    inv <- solve(covarF, cbind(1, muF, covarFB %*% wB, deparse.level=0L))
    w1 <- inv[,3]
    w2 <- inv[,1]
    w3 <- inv[,2]
    g1 <- sum(w3)  
    g2 <- sum(w2)
    g3 <- sum(wB)
    g4 <- sum(w1)
    g <- as.numeric((-lam*g1+(1-g3+g4))/g2)
    
    #2) compute weights
    list(wF = -w1+g*w2+lam*w3, gamma = g)
  },
  
  # computeLambda --------------------------------------------------------------
  computeLambda = function(get, i, bi.input){
    covarFB <- get$covarFB
    muF <- get$muF
    wB <- get$wB
    covarF <- get$covarF
    
    inv <- solve(covarF, cbind(1, muF, covarFB %*% wB))
    c1 <- sum(inv[,1])
    c2i <-  inv[i, 2]
    c3 <- sum(inv[,2])
    c4i <- inv[i, 1]
    Ci <- -c1*c2i + c3*c4i
    l1 <- sum(wB)
    l3 <- inv[,3]
    l2 <- sum(l3)
    
    if(length(bi.input)==1){                          # 1.bound to free
      if(Ci == 0) 0
      ((1-l1+l2)*c4i-c1*(bi.input+inv[i,3]))/Ci       # return lambda
    } else {                                          # 2.free to bound
      bi.lB <- bi.input[,1]
      bi.uB <- bi.input[,2]
      
      bi          <- bi.lB
      bi[Ci>0]    <- bi.uB[Ci>0]
      bi[Ci == 0] <- 0
      list(lambda = ((1-l1+l2)*c4i-c1*(bi+inv[i,3]))/Ci, bi = bi) 
      # return lambda and boundary
    }
  },
  
  # getMatrices --- -------------------------------------------------------------
  getMatrices = function(mu, covar, w, f){
    # Slice covarF,covarFB,covarB,muF,muB,wF,wB
    covarF <- covar[f,f]
    muF <- mu[f]
    b <- (seq_along(mu))[-f]
    covarFB <- covar[f,b]
    wB <- w[b]
    return(list(covarF = covarF, covarFB = covarFB, muF = muF, wB = wB))
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
  
  ## purgeChull
  purgeChull = function(weights_set, mu, covar){
    Sig2 <- colSums(weights_set *(covar %*% weights_set) ) # ms!
    Mu <- t(weights_set) %*% mu
    ch <- chull(cbind(Sig2, Mu))
    # start from largest mu, end at smallest sig2
    range <- c(which.max(Mu[ch]), which.min(Sig2[ch]))
    ch[max(range): min(range)]
  },
  
  cla.solve = function(cla.input, purge = TRUE){   #tol... = 1e-7  /sqrt(.Machine$double.eps)                                       
    # Compute the turning points, free sets and weights
    mu <- cla.input$mu
    covar <- cla.input$covar
    lB <- cla.input$lB
    uB <- cla.input$uB
    ans <- initAlgo(mu, lB, uB)
    f <- ans$index
    w <- ans$weights
    weights_set <- w  # store solution
    lambdas <- NA  # The first step has no lambda or gamma, add NA instead.
    gammas <- NA
    free_indices <- list(f) 
    lam <- 1 # set non-zero lam
    while ( lam!=0 && length(f) < length(mu)) {
      # 1) case a): Bound one free weight 
      l_in <- 0
      if(length(f) > 1 ){  
        get1 <- getMatrices(mu, covar, w, f)##
        compl <- computeLambda(get = get1, i = seq_along(f),
                                bi.input = cbind(lB[f], uB[f]))
        lam_in <- compl$lambda
        bi     <- compl$bi
        k <- which.max(lam_in)
        i_in  <-      f[k]
        bi_in <-     bi[k]
        l_in  <- lam_in[k]
      }  
      
      # 2) case b): Free one bounded weight
      b <- seq_along(mu)[-f]
      lam_out <- sapply(b, function(bi) {
        get_i <- getMatrices(mu, covar, w, c(f,bi))
        computeLambda(get = get_i, i = length(get_i$muF), bi.input = w[bi])
      })
      
      if (length(lambdas) > 1 && any(!(sml <- lam_out < lam*(1-1e-7)))) {
        lam_out <- lam_out[sml]
        b       <- b      [sml]
      }
      k <- which.max(lam_out)
      i_out <- b      [k] # one only !
      l_out <- lam_out[k]
      
      # 3) decide lambda  
      lam <- max(l_in, l_out, 0)
      if(lam > 0) { # remove i_in from f; or add i_out into f
        if(l_in > l_out ){
          f <- f[f != i_in]
          w[i_in] <- bi_in  # set value at the correct boundary
        } 
        else {
          f <- c(f,i_out)
        }
        getM <- getMatrices(mu, covar, w, f)
        compW <- computeW(lam, get = getM)
      }
      else{ #4) if max(l_in, l_out) < 0, "stop" when at the min var solution! 
        compW <- computeW(lam = lam, get = get1, muF = 0)
      } 
      
      wF <- compW$wF
      g <- compW$gamma
      w[f] <- wF[seq_along(f)] 
      
      lambdas <- c(lambdas, lam)
      weights_set <- cbind(weights_set, w, deparse.level = 0L) # store solution
      gammas <- c(gammas, g)
      free_indices <- c(free_indices, list(sort(f)))
    } #end While  
    
    free_indices[[length(free_indices)]] <- NA #delete?!
    
    if(purge){
      j <- purgeChull(weights_set, mu, covar)
       weights_set <- weights_set[,j] 
       gammas <- gammas[j]
       lambdas <- lambdas[j]
       free_indices <- free_indices[j]
    }
    
    MS.w <- MS(weights_set = weights_set, mu = mu, covar = covar)
    
    list(weights_set = weights_set, 
         free_indices=free_indices, 
         gammas = round(gammas, 3), lambdas = round(lambdas, 3),
         MS_weight = MS.w)
  }
)