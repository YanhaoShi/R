# 1) Portfolio optimization under a given risk--cccp package
#return risk, expected return and weights 
CCCP.solve <- function(cla.input, given.risk){ 
  mu <- cla.input$mu
  covar <- cla.input$covar
  lB <- cla.input$lB
  uB <- cla.input$uB
  N <- nrow(covar)
  sqrt.cov <- sqrm(covar)
  ## Portfolio risk constraint, under a given risk
  soc1 <- socc(F = sqrt.cov, g = rep(0, N), d = rep(0, N), f = given.risk)
  # socc(F, g, d, f),  || Fx + g ||_2 ¡Ü d'x + f.
  #                i.e.|| cov^1/2 w||_2 ¡Ü given risk(Sig)
  ## non-negativity constraint
  nno1 <- nnoc(G =rbind(-diag(N),diag(N)), h = c(lB, uB))
  ## Budget constraint
  A1 <- matrix(rep(1, N), nrow = 1)
  b1 <- 1.0
  ## optimization
  ans <- cccp(q =as.vector(mu)*(-1), A = A1, b = b1, 
              cList = list(nno1, soc1),  
              optctrl = ctrl(trace = FALSE)) 
  # weights sum up to 1, under above constraints
  getx(ans)
}

# 2) Quadratic Programming -- quadprog package
# to minimize 0.5 sqrt(wt cov w)-lambda(w mu), solve w, lambda from 0 to 1-
QP.solve <- function(cla.input, lambda){
  mu <- cla.input$mu
  covar <- cla.input$covar
  lB <- cla.input$lB
  uB <- cla.input$uB
  n <- nrow(mu)
  A <- cbind(-diag(n), diag(n), rep(1,n), rep(-1,n))
  b0 <- c(-uB, lB,1,-1) #uB; lB; sum up to 1
  solve.QP(covar * (1 - lambda), mu * lambda,A,b0)$solution
  # solve.QP(D, d, A, b)
  # D <- covar * (1-lambda)
  # d <- mu * lambda
  # A <- matrix; b0 (bvec) <- linear constraints
  # b -> w to solve
}


  

