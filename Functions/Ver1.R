## Previous version--Alexander Norring
Env1 <- funEnv(
  info = "Alexander Norring's",
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
  purgeNummErr = function(free_weights, gammas, lambdas, lB, uB, solution_set, tol){ 
    i <- 1
    while(TRUE){ 
      flag <- FALSE
      if(i == ncol(solution_set)){
        return(list(solution_set,lambdas,gammas,free_weights))
      }else{
        if((abs(sum(solution_set[i]))-1) > tol){ 
          flag <- TRUE
        } else {
          for(j in seq(length(solution_set[,i]))){
            if((solution_set[i][j,] - lB[j]) < (-tol) || (solution_set[i][j,] - uB[j]) > (tol)){
              flag <- TRUE
            }
          }
        }
        if(flag==TRUE){
          solution_set <- solution_set[,-i] 
          lambdas <- lambdas[-i]
          gammas <- gammas[-i]
          free_weights <- free_weights[-i] 
        }else{
          i <- i+1
        }
      }
    }
    return(list(solution_set,lambdas,gammas,free_weights))
    
  },
  
  # purgeExcess ---   -------------------------------------------------------------
  # Remove violations of the convex hull
  purgeExcess = function(mu, lambdas, gammas, free_weights, solution_set){
    i <- 1
    rep <- FALSE
    while(TRUE){
      if (rep == FALSE){
        i <- i+1
      }else{ #if flag == TRUE 
        i <- i #e.g. repeate again
      }
      if (i == length(solution_set)-1){ #stop i when at second last (not compare last with last)
        return(list(solution_set,lambdas, gammas, free_weights))
      }else{
        w <- solution_set[i]
        mu1 <- (t(w)%*% mu)[1]
        j <- i+1
        rep <- FALSE
        while(TRUE){
          if (j == length(solution_set)){
            break
          }else{
            w <- solution_set[j]
            mu_ <- (t(w)%*%mu)[1]
            if(mu1 < mu_){
              solution_set <- solution_set[,-i] 
              lambdas <- lambdas[-i] 
              gammas <- gammas[-i]
              free_weights <- free_weights[-i] 
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
  getMinVar = function(covar,solution_set){
    var <- data.frame()
    for(w in solution_set){
      a <- ((t(w) %*% covar) %*% w)
      var <- rbind(var,a)
    }
    return(list((min(var)**0.5),solution_set[which(var == min(var))]))  #! **? 0.5?
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
  getMatrices = function(mu, covar, solution_set, f){
    # Slice covarF,covarFB,covarB,muF,muB,wF,wB
    covarF <- reduceMatrix(covar,f,f) #covar[f,f]
    muF <- reduceMatrix(mu,f,0) #mu[f,]
    b <- getB(mu, f)
    covarFB <- reduceMatrix(covar,f,b) #covar[f,b]
    wB <- reduceMatrix(as.matrix(solution_set[,ncol(solution_set)]),b,0) ##?
    return(list(covarF,covarFB,muF,wB))
  },
  # cla.solve --- DONE -------------------------------------------------------------------
  cla.solve = function(cla.input){                                      
    # Compute the turning points, free sets and weights
    mu <- cla.input$mu
    covar <- cla.input$covar
    lB <- cla.input$lB
    uB <- cla.input$uB
    #initialize data.frames and vectors
    lambdas <- numeric(0) # lambdas 
    gammas <- numeric(0) # gammas
    free_weights <- data.frame() # free weights
    solution_set <- data.frame() # solution
    
    # Compute the turning points,free sets and weights
    temp_initAlgo <- initAlgo(mu, lB, uB)
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
        temp_getMat <- getMatrices(mu, covar,solution_set, f)
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
          temp_getMat <- getMatrices(mu, covar,solution_set,c(f,i))
          covarF <- temp_getMat[[1]]; covarFB <- temp_getMat[[2]]; muF <- temp_getMat[[3]]; wB <- temp_getMat[[4]]
          covarF_inv <- chol2inv(chol(covarF)) #solve(covarF)
          temp_compLam <- computeLambda(covarF_inv,covarFB,muF,wB,length(muF),solution_set[,ncol(solution_set)][[i]])
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
        temp_getMat <- getMatrices(mu, covar,solution_set,f)
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
        temp_getMat <- getMatrices(mu, covar,solution_set, f)
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
      solution_set <- cbind(solution_set, w) # store solution
      gammas <- append(gammas, g)
      free_weights <- append(free_weights, as.data.frame(f))
      if(lambdas[length(lambdas)] == 0) {
        break 
      }
    } #end While
    
    
    temp_pNE <- purgeNummErr(free_weights, gammas, lambdas, lB, uB, solution_set, 1e-09) 
    solution_set <- temp_pNE[[1]]
    lambdas <- temp_pNE[[2]]
    gammas <- temp_pNE[[3]]
    free_weights <- temp_pNE[[4]]
    
    temp_purgeEx <- purgeExcess(mu, lambdas, gammas, free_weights, solution_set)
    solution_set <- temp_purgeEx[[1]]
    lambdas <- temp_purgeEx[[2]]
    gammas <- temp_purgeEx[[3]]
    free_weights <- temp_purgeEx[[4]]
    
    ans <- list(solution_set, free_weights, gammas, lambdas)
    return(ans)
  }
)
