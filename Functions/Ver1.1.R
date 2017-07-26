## Previous version--Alexander Norring
Env1.1 <- funEnv(
  info = "AN - adjusted, same as previous",
  initAlgo = function(mu, lB, uB ){             #! mu 
    # Initialize the algo 
    #1) Form data.frame
    a <- data.frame(seq(1:length(mu)),mu)      #! n <- #assets 
    colnames(a) <- c('id','mu')
    temp <- order(a[ ,2])                          #! index.ordered <-...
    #2) Sort data.frame                 #! as in reference, should be decreasing order
    b <- a[temp, ] 
    
    #3) First free weight
    temp <- list(length(mu), lB)     #! delete temp
    i <- temp[[1]] + 1  
    w <- temp[[2]] 
    while(sum(w) < 1) {
      i <- i - 1 
      w[b[i, 1], ] <- uB[b[i, 1], ]  
    }
    w[b[i, 1], ] <- w[b[i, 1], ] + (1 - sum(w))
    
    return (list(b[i, 1], w, i))
  },
  
  
  # getB ---  --------------------------------------------------------------------
  # Deriving B from F
  getB = function(mu,f){
    b <- seq(1:length(mu))[-c(f)] #drop f's
    return (b)
  },
  
  # reduceMatrix ---  ------------------------------------------------------------
  reduceMatrix =function(input_matrix, listX, listY){
    # Reduce a matrix to the provided list of rows and columns
    #if either lists are empty then end function call	
    if(length(listX)==0 || length(listY)==0){ 
      return(as.null(input_matrix))
    }
    if(length(listX)==1){ # handle (1xn) matrix 
      intermediate_matrix <- t(as.matrix(input_matrix[listX,]))
      if(sum(listY) == 0){ #handle (nx1) matrix 
        return_matrix <- as.matrix(intermediate_matrix[,])
      } else {
        return_matrix <- as.matrix(intermediate_matrix[,listY])
      }
      return(t(return_matrix))
    }else{ 
      intermediate_matrix <- as.matrix(input_matrix[listX,])
      if(sum(listY) == 0){ #handle (nx1) matrix
        return_matrix <- as.matrix(intermediate_matrix[,])
      } else {
        return_matrix <- as.matrix(intermediate_matrix[,listY]) 
      }
      return(return_matrix)
    }
  },
  
  # purgeNumErr ---  -------------------------------------------------------------
  # Purge violations of inequality constraints (associated with ill-conditioned covar matrix)
  # tolerance, e.g 10e-10
  purgeNumErr = function(free_indices, gammas, lambdas, lB, uB, weights_set, tol){ 
    i <- 1
    while(TRUE){ 
      flag <- FALSE
      if(i == ncol(weights_set)){
        return(list(weights_set,lambdas,gammas,free_indices))
      }else{
        if((abs(sum(weights_set[i]))-1) > tol){ 
          flag <- TRUE
        } else {
          for(j in seq(length(weights_set[,i]))){
            if((weights_set[j,i] - lB[j]) < (-tol) || (weights_set[j,i] - uB[j]) > (tol)){
              flag <- TRUE
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
      }
    }
    return(list(weights_set,lambdas,gammas,free_indices))
    
  },
  
  # purgeExcess ---   -------------------------------------------------------------
  # Remove violations of the convex hull
  purgeExcess = function(mu, lambdas, gammas, free_indices, weights_set){
    i <- 1
    rep <- FALSE
    while(TRUE){
      if (rep == FALSE){
        i <- i+1
      }else{ #if flag == TRUE 
        i <- i #e.g. repeate again
      }
      if (i == length(weights_set)-1){ #stop i when at second last (not compare last with last)
        return(list(weights_set,lambdas, gammas, free_indices))
      }else{
        w <- weights_set[i]
        mu1 <- (t(w)%*%mu)[1]
        j <- i+1
        rep <- FALSE
        while(TRUE){
          if (j == length(weights_set)){
            break
          }else{
            w <- weights_set[j]
            mu_ <- (t(w)%*%mu)[1]
            if(mu1 < mu_){
              weights_set <- weights_set[,-i] 
              lambdas <- lambdas[-i] 
              gammas <- gammas[-i]
              free_indices <- free_indices[-i] 
              rep <- TRUE
              break
            }else{
              j <- j+1
            }
          }
        }
      }
    }
  },
  
  # getMinVar ---  ---------------------------------------------------------------
  # Get the minimum variance solution
  getMinVar = function(covar,weights_set){
    var <- data.frame()
    for(w in weights_set){
      a <- ((t(w) %*% covar) %*% w)
      var <- rbind(var,a)
    }
    return(list((min(var)**0.5),weights_set[which(var == min(var))]))  #! **? 0.5?
  },
  
  # efFrontier -----------------------------------------------------------------
  # Get the efficient frontier
  efFrontier = function(mu, covar, weights_set, dataPoints){
    mu <- c(); sigma <- c(); weights <- c()
    a <- seq(0,1,length=round(dataPoints/length(weights_set))); a <- a[a != 1] # remove the 1, to avoid duplications
    b <- seq(1,length(weights_set)-1)
    for(i in b){
      w0 <- weights_set[i]; w1 <- weights_set[i+1]
      if(i==b[length(b)]){
        a <- seq(0,1,length=round(dataPoints/length(weights_set))) # include the 1 in the last iteration
      }
      for(j in a){
        w <- as.matrix(w1*j+(1-j)*w0)
        weights <- cbind(weights,c(w))
        mu <- append(mu, (t(w)%*%mu)[1])
        sigma <- append(sigma, ((t(w) %*% covar) %*% w)**0.5)
      }
    }
    return(list(mu, sigma, weights))
  },
  # computeBi ----------------------------------------------------------------
  # needed in computeLambda function (Niedermayer and Niedermayer 2007, 10)
  computeBi =function(c,bi){
    if(c>0){
      bi <- bi[[2]] #ui <- defined bi as bi[[2]]<-ui
    }
    if(c<0){
      bi <- bi[[1]] #li  <- defined bi as bi[[1]]<-li
    }
    return(bi)
  },
  
  # computeW -------------------------------------------------------------------
  # compute gamma #(Bailey and Lopez de Prado 2013, 11)
  computeW = function(lambdas,covarF_inv,covarFB,muF,wB){
    #1) compute gamma
    onesF <- matrix(1,length(muF)) # (kx1)-vector
    g1 <- (t(onesF) %*% covarF_inv) %*% muF
    g2 <- (t(onesF) %*% covarF_inv) %*% onesF
    if(is.null(wB)){
      g <- as.numeric(-tail(lambdas,n=1)*g1/g2+1/g2) 
      w1 <- 0
    } else {
      onesB <- matrix(1,length(wB))
      g3 <- t(onesB) %*% wB
      g4 <- covarF_inv %*% covarFB
      w1 <- g4 %*% wB
      g4 <- t(onesF) %*% w1
      g <- as.numeric(-tail(lambdas,n=1)*g1/g2+(1-g3+g4)/g2)
    }
    
    #2) compute weights
    w2 <- covarF_inv %*% onesF
    w3 <- covarF_inv %*% muF
    
    return(list(-w1+g*w2+tail(lambdas,n=1)*w3,g))
  },
  
  # computeLambda --------------------------------------------------------------
  # (Niedermayer and Niedermayer 2007, 10)
  computeLambda = function(covarF_inv,covarFB,muF,wB,i,bi){
    #1) C
    onesF <- matrix(1,length(muF)) # (kx1)-vector
    c1 <- (t(onesF) %*% covarF_inv) %*% onesF
    c2 <- covarF_inv %*% muF
    c3 <- (t(onesF) %*% covarF_inv) %*% muF
    c4 <- covarF_inv %*% onesF
    c <- (-c1)*c2[i]+c3*c4[i]
    if(c==0){
      return(list(0,0))
    }
    #2) bi
    if(is.list(bi)){
      bi <- computeBi(c,bi)
    }
    #3) lambda
    if(length(wB)==0){
      # All free assets 
      return(list(as.numeric((c4[i]-c1*bi)/c),bi))
    }else{
      onesB <- matrix(1,length(wB))
      l1 <- t(onesB) %*% wB
      l2 <- covarF_inv %*% covarFB
      l3 <- l2 %*% wB
      l2 <- t(onesF) %*% l3
      return(list(as.numeric(((1-l1+l2)*c4[i]-c1*(bi+l3[i]))/c),bi))
    }
  },
  
  # getMatrices --- DONE -------------------------------------------------------------
  getMatrices = function(mu, covar, weights_set, f){
    # Slice covarF,covarFB,covarB,muF,muB,wF,wB
    covarF <- reduceMatrix(covar,f,f) #covar[f,f]
    muF <- reduceMatrix(mu,f,0) #mu[f,]
    b <- getB(mu, f)
    covarFB <- reduceMatrix(covar,f,b) #covar[f,b]
    wB <- reduceMatrix(as.matrix(weights_set[,ncol(weights_set)]),b,0) ##?
    return(list(covarF,covarFB,muF,wB))
  },
  # cla.solve --- DONE -------------------------------------------------------------------
  cla.solve = function(cla.input){ 
    mu <- cla.input$mu
    covar <- cla.input$covar
    lB <- cla.input$lB
    uB <- cla.input$uB
    #initialize data.frames and vectors
    lambdas <- numeric(0) # lambdas 
    gammas <- numeric(0) # gammas
    free_indices <- data.frame() # free weights
    weights_set <- data.frame() # weights
    
    # Compute the turning points,free sets and weights
    temp_initAlgo <- initAlgo(mu, lB, uB)
    f <- temp_initAlgo[[1]]
    w <- temp_initAlgo[[2]]
    weights_set <- rbind(weights_set, w) # store weights
    lambdas <- append(lambdas, NA_integer_) 
    gammas <- append(gammas, NA_integer_) 
    free_indices <- as.data.frame(f)
    
    
    while ( TRUE ) {
      
      # 1) case a): Bound one free weight 
      l_in <- 0 
      if (length(f) > 1) {
        temp_getMat <- getMatrices(mu, covar,weights_set, f)
        covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; muF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
        covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
        j <- 1
        for (i in f) {
          temp_compLam <- computeLambda(covarF_inv,covarFB,muF,wB,j,as.list(c(lB[i],uB[i])))
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
      if(length(f) < length(mu)) {
        b <- getB(mu, f)
        for(i in b){
          temp_getMat <- getMatrices(mu, covar,weights_set,c(f,i))
          covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; muF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
          covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
          temp_compLam <- computeLambda(covarF_inv,covarFB,muF,wB,length(muF),weights_set[,ncol(weights_set)][[i]])
          l <- temp_compLam[[1]]; bi <- temp_compLam[[2]]
          
          if ((is.na(lambdas[length(lambdas)]) || l<lambdas[length(lambdas)]) && (l>l_out)){
            l_out <- l
            i_out <- i
          }
        }
      }
      
      if( (l_in == 0 || l_in < 0) & (l_out == 0 || l_out < 0) ) { # "stop" when at the min var weights! 
        
        # 3) compute minimum variance weights
        lambdas <- append(lambdas,0)
        temp_getMat <- getMatrices(mu, covar,weights_set,f)
        covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; muF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
        covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
        muF <- matrix(0,length(muF))
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
        temp_getMat <- getMatrices(mu, covar,weights_set, f)
        covarF <- temp_getMat[[1]]
        covarFB <- temp_getMat[[2]]
        muF <- temp_getMat[[3]]
        wB <- temp_getMat[[4]]
        covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
      }
      
      # 5) compute solution vector
      temp_compW <- computeW(lambdas,covarF_inv,covarFB,muF,wB)
      wF <- temp_compW[[1]]
      g <- temp_compW[[2]]
      for(i in seq(length(f))){
        w[f[i]]=wF[i]
      }
      weights_set <- cbind(weights_set, w) # store solution
      gammas <- append(gammas, g)
      free_indices <- append(free_indices, as.data.frame(f))
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