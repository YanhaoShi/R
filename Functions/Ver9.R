Env9 <- funEnv( 
  info = "covF_inv: add/remove i from previous covF_inv", 
  # Initialize the weight -- Find first free weight
  initAlgo = function(mu, lB, uB ){
    #New-ordered return, lB, uB with decreasing return  
    w <- c()
    index.new <- order(mu,decreasing = TRUE) # new order with decreasing return
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
  
  computeInv = function(covF.inv.pre, covar, f, i, 
                         add = FALSE, remove = FALSE, mu, w){
    if(add){  # B -> F
      a <- covar[f, i]
      cc <- covF.inv.pre %*% a
      b <- covar[i, i] - sum(a * cc)
      e1 <-  cc%*%t(cc)/b
      m <- cc/b
      covF.inv <- rbind(cbind(covF.inv.pre + e1, -m),
                        cbind(-t(m), 1/b))
      f.new <- c(f, i)
    }
    if(remove){
      ii <- which(f == i)
      
      B <- covF.inv.pre[-ii, -ii]
      b <- covF.inv.pre[ -ii, ii]
      beta <- covF.inv.pre[ii, ii]
      covF.inv <- B - 1/beta * b %*%t(b) 
      f.new <- f[-ii]
    }
    covarF <- covar[f.new, f.new]
    muF <- mu[f.new]
    b <- (seq_along(mu))[-f.new]
    covarFB <- covar[f.new,b]
    wB <- w[b]
    inv <- covF.inv %*% cbind(1, muF, covarFB %*% wB, deparse.level=0L)
    list(inv = inv, covF.inv = covF.inv)
  },
  
  # computeW -----------------------------------------------------------------
  computeW = function(lam, inv, wB){
    # w2 <- inv[,1]; w3 <- inv[,2]; w1 <- inv[,3]
    inv.s <- colSums(inv) # g1 <- inv.s[2]; g2 <- inv.s[1]; g4 <- inv.s[3]
    #1) compute gamma
    g <- (-lam * inv.s[2] + (1- sum(wB) + inv.s[3]))/inv.s[1]
    #2) compute free weights
    list(wF = - inv[,3] + g * inv[,1] + lam * inv[,2], gamma = g)
  },
  
  # computeLambda --------------------------------------------------------------
  computeLambda = function(wB, inv, i, bi.input){
    inv.s <- colSums(inv) 
    # c1 <- inv.s[1]; l2 <- inv.s[3]; c2i <- inv[i,2];
    # c3 <- inv.s[2]; c4i <- inv[i, 1]; l1 <- sum(wB)

    c1 <- inv.s[1]
    if(length(bi.input)==1){                          # 1.bound to free
      c4i <- inv[i, 1]
      Ci <- - c1 * inv[i, 2] + inv.s[2] * c4i
      if(Ci == 0) 0
      ((1- sum(wB) + inv.s[3])* c4i- c1 * (bi.input + inv[i, 3]))/Ci # return lambda
    } else {                                          # 2.free to bound
      c4i <- inv[, 1]
      Ci <- - c1 * inv[, 2] + inv.s[2] * c4i
      bi          <- bi.input[i, 1]         # bi.lB
      bi[Ci > 0]  <- bi.input[i[Ci > 0], 2]  # bi.uB
      bi[Ci == 0] <- 0
      list(lambda = ((1- sum(wB) + inv.s[3]) * c4i- c1 *(bi + inv[, 3]))/Ci, 
           bi = bi) 
      # return lambda and boundary
    }
  },
  
  
  MS = function(weights_set, mu, covar){ 
    Sig2 <- colSums(weights_set *(covar %*% weights_set) )
    cbind(Sig = sqrt(Sig2), Mu = as.vector(t(weights_set) %*% mu))
  },
  
  cla.solve = function(cla.input){   #tol... = 1e-7  /sqrt(.Machine$double.eps) 
    options(digits = 3)
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
    covFinv <- 1/covar[f, f]

    while ( lam > 0 && length(f) < length(mu)) {
      # 1) case a): Bound one free weight  F -> B
      l_in <- 0
      if(length(f) > 1 ){  
        compl <- computeLambda(wB = w[-f], inv = inv,  # inv from last step k (k > 1)
                               i = f, bi.input = cbind(lB, uB))
        lam_in <- compl$lambda
        bi     <- compl$bi
        k <- which.max(lam_in)
        i_in  <-      f[k]
        bi_in <-     bi[k]
        l_in  <- lam_in[k]
      }  
      
      # 2) case b): Free one bounded weight  B -> F
      b <- seq_along(mu)[-f]
      inv_list <- lapply(b, function(bi){
        computeInv(covF.inv.pre = covFinv, covar, f, bi, add = TRUE,
                   mu = mu, w = w)
      })
      
      fi <- length(f) + 1
      lam_out <- sapply(seq_along(b), function(i) {
        computeLambda(wB = w[b[-i]], inv = inv_list[[i]]$inv,
                      i = fi, bi.input = w[b[i]])
      })
      
      if (length(lambdas) > 1 && any(!(sml <- lam_out < lam*(1-1e-7)))) {## tol
        lam_out <- lam_out[sml]
        b       <- b      [sml]
        inv_list <- inv_list[sml]
      }
      k <- which.max(lam_out)
      i_out <- b      [k] # one only !
      l_out <- lam_out[k]
      inv_out <- inv_list[[k]]$inv
      getM <- getMatrices(mu, covar, w, c(f, i_out))
      
      covFinv.out <- inv_list[[k]][[2]]
      
      # 3) decide lambda  
      lam <- max(l_in, l_out, 0)
      if(lam > 0) { # remove i_in from f; or add i_out into f
        if(l_in > l_out ){
          w[i_in] <- bi_in  # set value at the correct boundary
          a <- computeInv(covF.inv.pre = covFinv, covar, 
                          f, i_in, remove = TRUE, mu = mu, w = w)
          f <- f[f != i_in]
          covFinv <- a$covF.inv
          inv <- a$inv
          
        } 
        else {
          f <- c(f,i_out)
          inv <- inv_out
          covFinv <- covFinv.out
        }
        compW <- computeW(lam, inv = inv, wB = w[-f])
      }
      else{ #4) if max(l_in, l_out) < 0, "stop" when at the min var solution! 
        compW <- computeW(lam = lam, inv = inv, wB = w[-f])
        # muF = 0 not necessary, get1 replaced by getM (ie getM from previous step)
      } 
      
      wF <- compW$wF
      g <- compW$gamma
      w[f] <- wF[seq_along(f)] 
      
      lambdas <- c(lambdas, lam)
      weights_set <- cbind(weights_set, w, deparse.level = 0L) # store solution
      gammas <- c(gammas, g)
      free_indices <- c(free_indices, list(sort(f)))
    } #end While  
    
    list(weights_set = weights_set, 
         free_indices = free_indices, 
         gammas = gammas, lambdas = lambdas,
         MS_weight = MS(weights_set = weights_set, mu = mu, covar = covar))
  }
)