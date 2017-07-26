# Compare results

# compare 3 versions of CLA time and results
# cla, qp vs. cccp
w.compare <- function(cla.result, assets, rep = 3){  #as.matrix(lB or uB)!
  mu <- assets$mu
  covar <- assets$covar
  lB <- assets$lB
  uB <- assets$uB
  w.cla <- cla.result$weights_set
  n <- ncol(w.cla)
  lam <- ((0:n)/n)^2
  ms.cla <- MS(w.cla, mu, covar)
  
  t.qp <- cbind(system.time(w.qp <- lapply(lam, 
            function(x) return(tryCatch(QP.solve(assets, x), # skip error, pick "correct" mu
            error=function(e) NULL)))),
    replicate(rep, system.time(lapply(lam, 
                        function(x) return(tryCatch(QP.solve(assets, x), # skip error, pick "correct" mu
                                                  error=function(e) NULL)))) ))
  t.qp <- apply(t(t.qp), 2, mean)[1:3]              
                    
                    
  t.cccp <- cbind(system.time(w.cccp <- sapply(ms.cla[,"Sig"], 
                                      function(x) CCCP.solve (assets, x))),
           replicate(rep, system.time(sapply(ms.cla[,"Sig"], 
                        function(x) CCCP.solve (assets, x)))))
  t.cccp <-  apply(t(t.cccp), 2, mean)[1:3]                      
  
  ind.qp <- which(!sapply(w.qp,is.null))
  w.qp <- sapply(ind.qp, function(x) w.qp[[x]]) # remove NULL terms
  ind.cccp <- which(!is.na(w.cccp[1,]))
  w.cccp <- w.cccp[, ind.cccp]  # remove NaN
  
  list(weights = list(cla = w.cla, qp = w.qp, cccp = w.cccp),
       MS = list( ms.cla = ms.cla,
                  ms.qp = MS(w.qp, mu, covar),
                  ms.cccp = MS(w.cccp, mu, covar)),
       ind.qp = ind.qp, ind.cccp = ind.cccp)# compare micro?
}

compareEF <- function(w.qp, w.cla, Mu, Cov, tol = -16){
  ms.qp <- MS(w.qp, Mu, Cov)
  ms.cla <- MS(w.cla, Mu, Cov)
  if(!anyNA(ms.qp, ms.cla)){ 
    ind <- (ms.qp[, "Mu"] <= max(ms.cla[, "Mu"]))&(ms.qp[, "Mu"] >= min(ms.cla[, "Mu"]))
    Sig2.cla <- sapply(which(ind), function(i){j <- sum(ms.cla[, "Mu"] > ms.qp[i, "Mu"])
    hyperbolic(ms.qp[i, "Mu"], w.cla[,j], w.cla[,j+1], Mu, Cov)})
    Sig2.qp <- ms.qp[ind, "Sig2"]
    cbind(index = which(ind),  Sig2.cla, Sig2.qp,
          diff = Sig2.cla-Sig2.qp,  # Sig2.cla - Sig2.qp
          tol = Sig2.cla-Sig2.qp > -10^tol )
  } else NA
}
