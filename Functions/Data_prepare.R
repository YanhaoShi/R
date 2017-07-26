# 1. Company data
 # Raw company data
if(!file.exists(d.file("company_assets.rds", exists = FALSE))){
  saveRDS(read.table(d.file('03-02_CLA_Data_Tot.csv'), header = TRUE, sep = ','), 
          d.file("company_assets.rds", exists = FALSE))
}
company_assets <- readRDS(d.file("company_assets.rds", exists = FALSE))

 # Total assets: assets 1457 
if(!file.exists(d.file("assets1457.rds", exists = FALSE))){
  transData <- function(d){
    n <- ncol(company_assets)
    assets <- as.matrix(t(d[1:3, 1:n]))
    colnames(assets) <- c("mu", "lB", "uB")
    covar <- as.matrix(d[4:(n+3), 1:n ])
    rownames(covar) <- colnames(covar)
    list(assets = assets, covar = covar)
  }
  # for function input
  getData <- function(d){
    row.names(d$assets) <- NULL
    mu <- as.matrix(d$assets[, "mu"])
    lB <- as.matrix(d$assets[, "lB"])
    uB <- as.matrix(d$assets[, "uB"])
    colnames(d$covar) <- NULL
    rownames(d$covar) <- NULL
    covar <- d$covar
    list(mu = mu, lB = lB, uB = uB, covar = covar)
  }
  
  if(!file.exists(d.file("trans_company_assets.rds", exists = FALSE))){
    saveRDS(transData(company_assets), 
            d.file("trans_company_assets.rds", exists = FALSE))
  }
  saveRDS(getData(transData(company_assets)), d.file("assets1457", exists = FALSE))
}

assets1457 <- readRDS(d.file("assets1457.rds", exists = FALSE))
######################################################################################
# 2. FRAPO Index
# the corresponding assets under a given set of indices from the total assets
Index.ms2 <- function(Index){
  nr <- nrow(Index) 
  nc <- ncol(Index)
  lr <- log(Index[-1, ]/Index[-nr,]) #log return
  loss <- -lr
  gfit <- lapply(loss, garchFit,  formula = ~ garch(1,1),
                 cond.dist = "std", trace = FALSE)
  gpredict <- t(sapply(gfit, function(x) predict(x, n.ahead = 1)))
  mu <- -unlist(gpredict[, "meanForecast"])
  std <- unlist(gpredict[, "standardDeviation"])
  index.stand <- sapply(seq(nc), function(x) lr[,x]/gfit[[x]]@sigma.t)
  diag.std <- diag(std)
  covar <- diag.std %*% cor(index.stand) %*% diag.std
  list(mu = mu, covar = covar)
}

transIndex <- function(Index, lB = as.matrix(rep(0, ncol(Index))), 
                       uB = as.matrix(rep(0.075, ncol(Index)))){
  ms2 <- Index.ms2(Index)
  mu <- as.matrix(unname(ms2$mu))
  covar <- unname(ms2$covar)
  list(mu = mu, covar = covar, lB = lB, uB = uB)
}

GetIndex <- function(ind, total_assets, lB =  total_assets$lB[ind, ],
                     uB = total_assets$uB[ind, ]){
  list(mu = total_assets$mu[ind, , drop = FALSE], 
       covar = total_assets$covar[ind, ind, drop = FALSE],
       lB = as.matrix(lB), 
       uB = as.matrix(uB))
}

if(!file.exists(d.file("trans_sp.rds", exists = FALSE))){
  saveRDS(transIndex(SP500), d.file("trans_sp.rds", exists = FALSE))
}

if(!file.exists(d.file("trans_nasdaq.rds", exists = FALSE))){
  saveRDS(transIndex(NASDAQ), d.file("trans_nasdaq.rds", exists = FALSE))
}

trans.sp <- readRDS(d.file("trans_sp.rds", exists = FALSE))
trans.nasdaq <- readRDS(d.file("trans_nasdaq.rds", exists = FALSE))
#######################################################################################
# 3. readRDS


index.name <- c("company", "sp", "nasdaq")
set.name <- c("Company", "SP", "NASDAQ")
index.n <- c(1457, length(trans.sp$mu), length(trans.nasdaq$mu))
index.total <- c("assets1457", "trans.sp", "trans.nasdaq")

ind.lab <- lapply(1:3, function(i) exp_ind(1.2, 50, index.n[i])$ind)
names(ind.lab) <- index.name
## ind.company, ind.sp and ind.nasdaq
## label.company, label.sp and label.nasdaq
## company.cla.t, sp.cla.t                     #line 111/134/143
## company.cla.mt, sp.cla.mt
## covF.company, covF.sp
## Company55_set1..., SP55_set1,...
## analysis.company, analysis.sp
for(i in 1:3){ 
  # read ind.company, ind.sp and ind.nasdaq
  l <- length(ind.lab[[i]])
  #at least 50 assets, at most n assets, index power of 1.2
  # if 1.1, too much to plot and compute; if 1.3 or larger, too little information.
  if(!file.exists(d.file(paste("ind_", index.name[i], ".rds", sep = ""), exists = FALSE))){ 
    set.seed(2017)
    ind.list <- list()
    for(j in seq(l)){
      ind.list[[j]] <- lapply(rep(ind.lab[[i]][j], 50), function(x) 
        sample.int(n, x, replace = FALSE))
    }
    saveRDS(ind.list, d.file(paste("ind_", index.name[i], ".rds", sep = ""), exists = FALSE))
  }
  ind.list <- readRDS( d.file(paste("ind_", index.name[i], ".rds", sep = ""), exists = FALSE))
  assign(paste0("ind.", index.name[i]), ind.list)
  assign(paste0("label.", index.name[i]), ind.lab[[i]])
} 
### Env6!!! not the newest
for(i in 1:3){ 
  ind.label <- exp_ind(1.2, 50, index.n[i])$ind 
  l <- length(ind.label)
  # read company_cla_t, sp_cla_t, nasdaq_cla_t
  if(!file.exists(d.file(paste0(index.name[i], "_cla_t7.rds"), exists = FALSE))){
    cla.t <- c()
    for(j in seq(l)){
      m <- ind.label[j]
      sys.T <- sapply(1:10, function(k){
        assets <- GetIndex(ind.list[[j]][[k]], get(index.total[i]))
        sys.t.rep <- replicate(9, system.time(Env6$cla.solve(assets$mu,assets$covar, 
                                                             assets$lB, assets$uB)))
        sys.t.rep <- cbind(sys.t.rep, system.time(
          result <- Env6$cla.solve(assets)))
        sys.t <- apply(sys.t.rep[1:3, ], 1, mean)
        saveRDS(result, d.file(paste(set.name[i], m, "_set", k, ".rds", sep = ""), 
                               exists = FALSE))
        sys.t})
      sys.T <- t(sys.T) 
      rownames(sys.T) <- paste0(m, "set", 1:10)
      cla.t <- rbind(cla.t, sys.T)
    }  
    cla.mt <- t(sapply(seq(l), function(m) 
      apply(cla.t[1:10 + (m - 1) * 10, ], 2, mean)))
    row.names(cla.mt) <- ind.label[seq(l)]
    
    cla.t7 <- list()
    cla.t7$time <- cla.t
    cla.t7$time.mean <- cla.mt
    saveRDS(cla.t7, d.file(paste0(index.name[i], "_cla_t7.rds"), exists = FALSE))
  }
  assign(paste0(index.name[i], ".cla.t7"), 
         readRDS(d.file(paste0(index.name[i], "_cla_t7.rds"), exists = FALSE)))
}  

for(i in 1:3){ 
  ind.label <- exp_ind(1.2, 50, index.n[i])$ind 
  l <- length(ind.label)
  # read company_cla_t_pre... using cla.solve form Ver1.R
  if(!file.exists(d.file(paste0(index.name[i], "_cla_t1.rds"), exists = FALSE))){
    cla.t <- c()
    for(j in seq(l)){
      m <- ind.label[j]
      sys.T <- sapply(1:10, function(k){
        assets <- GetIndex(ind.list[[j]][[k]], get(index.total[i]))
        sys.t.rep <- replicate(3, system.time(Env1$cla.solve(assets)))
        sys.t <- apply(sys.t.rep[1:3, ], 1, mean)
        sys.t})
      sys.T <- t(sys.T) 
      rownames(sys.T) <- paste0(m, "set", 1:10)
      cla.t <- rbind(cla.t, sys.T)
    }  
    l <- nrow(cla.t)/10
    cla.mt <- t(sapply(seq(l), function(m) 
      apply(cla.t[1:10 + (m - 1) * 10, ], 2, mean)))
    row.names(cla.mt) <- ind.label[seq(l)]
    
    cla.t1 <- list()
    cla.t1$time <- cla.t
    cla.t1$time.mean <- cla.mt
    saveRDS(cla.t1, d.file(paste0(index.name[i], "_cla_t1.rds"), exists = FALSE))
  }
  assign(paste0(index.name[i], ".cla.t1"), 
         readRDS(d.file(paste0(index.name[i], "_cla_t1.rds"), exists = FALSE)))
} 




for(i in 1:3){ 
  ind.label <- exp_ind(1.2, 50, index.n[i])$ind 
  l <- length(ind.label)
  # read time cost on inversing covF
  if(!file.exists(d.file(paste0("covF_inv_", index.name[i], ".rds"), exists = FALSE))){
    saveRDS(testInv(ind.list = ind.company, assets1457),
            d.file(paste0("covF_inv_", index.name[i], ".rds"), exists = FALSE))
  }
  assign(paste0("covF.", index.name[i]), 
         readRDS(d.file(paste0("covF_inv_", index.name[i], ".rds"), exists = FALSE)))
  # read Index55_set1,...
 # for(k in seq_along(ind.label)){
 #   for(j in 1:10){
 #     if( !file.exists(d.file(paste0(set.name[i], ind.label[k], "_set", j, ".rds", sep = ""), 
 #                             exists = FALSE))){ 
 #       assets <- GetIndex(ind.list[[k]][[j]], get(index.total[i]))
 #       result <- Env6$cla.solve(assets$mu, assets$covar, assets$lB, assets$uB)
 #       saveRDS(result, d.file(paste0(set.name[i], ind.label[k], "_set", j, ".rds", sep = ""), 
 #                              exists = FALSE)) 
 #     }
 #   }
 # }

  # read anlysis.company--------------
  if(!file.exists(d.file(paste0("analysis_", index.name[i],".rds"), exists = FALSE))){ 
    cla <- list()
    t_compare <- list()
    nweights <- matrix(0, ncol = l, nrow = 10)
    nweights_unpurge <- matrix(0, ncol = l, nrow = 10)
    covF.kappa <- matrix(0, ncol = l, nrow = 10)
    covF.rcond <- matrix(0, ncol = l, nrow = 10)
    cov.check <- matrix(0, ncol = l, nrow = 10)
    covF <- list() # 8 sub-list, each with 10 sub-list, with cla-result of covFs
    lambda <- list()
    
    for(k in seq_along(ind.label)){
      re <- lapply(1:10, function(j){
        assets <- GetIndex(ind.list[[k]][[j]], get(index.total[i])) 
        w.compare2(assets$mu, assets$covar, assets$lB, assets$uB)})
      
      cla[[k]] <- lapply(1:10, function(x)re[[x]]$result.cla)
      lambda[[k]] <- lapply(1:10, function(x)re[[x]]$lambda.qp)
      t_compare[[k]] <- sapply(1:10, function(x) {micro <- re[[x]]$micro
      unname(sapply(levels(micro$expr), 
                    function(y){median(micro$time[micro$expr==y])*1e-6}))}) 
      #take the median, miliseconds 
      
      nweights[, k] <- sapply(1:10, function(x) ncol(cla[[k]][[x]]$weights_set_purge))
      nweights_unpurge[,k] <- sapply(1:10, function(x) ncol(cla[[k]][[x]]$weights_set))
      
      covF[[k]] <- lapply(1:10, function(x) lapply(cla[[k]][[x]]$covarF, as.matrix))
      covF.kappa[, k] <- sapply(1:10,function(x) max(sapply(covF[[k]][[x]], kappa)))
      covF.rcond[, k] <- sapply(1:10,function(x) min(sapply(covF[[k]][[x]], rcond)))
      cov.check[, k] <- sapply(1:10, function(x) is.positive.definite(readRDS(d.file(
        paste(toupper(name), length(ind.list[[k]][[x]]), "_data", 1, ".rds", sep = ""), 
        exists = FALSE))$covar) ) #check if the covariance matrix is positive definite
    }
    colnames(nweights) <- ind.label
    colnames(nweights_unpurge) <- ind.label
    colnames(cov.check) <- ind.label
    colnames(covF.kappa) <- ind.label
    colnames(covF.rcond) <- ind.label
    t_mean <- sapply(seq_along(ind.label), 
                     function(x) apply(t_compare[[x]], 1, mean))
    colnames(t_mean) <- ind.label
    rownames(t_mean) <- c("cla", "qp", "cccp")[1:nrow(t_mean)]
    
    analysis <- list(cla_result = cla,
                     t_compare = t_compare,
                     t_mean = t_mean,
                     nweights = nweights,
                     nweights_unpurge = nweights_unpurge,
                     covF.kappa = covF.kappa,
                     covF.rcond = covF.rcond,
                     covF = covF,
                     cov.check = cov.check,
                     lambda.QP = lambda)
    saveRDS(analysis, d.file(paste0("analysis_", index.name[i],".rds"), exists = FALSE))
  }
  assign(paste0("analysis.", index.name[i]),
         readRDS(d.file(paste0("analysis_", index.name[i],".rds"), exists = FALSE)))
  
  
  
}

## read cla results
## read unique_purge_sp. Thesis: testing purgeChull and purgeUnique
if(!file.exists(d.file("purge_unique_sp.rds", exists = FALSE))){
  time <- c()
  diff.num <- c()
  for(i in seq(label.sp)){
    for(j in seq(10)){
      u <- testUnique(GetIndex(ind.sp[[i]][[j]], trans.sp))
      time <- rbind(time, u$time)
      diff.num <- rbind(diff.num, u$diff.num)
    }
  }
  rname <- paste0(rep(label.sp, each = 10), "set", 1:10)
  row.names(time) <- rname
  row.names(diff.num) <- rname
  time.mean <- t(sapply(seq_along(label.company), function(i)
    apply(time[1:10 + (i-1)*10, ], 2, mean)))
  row.names(time.mean) <- label.company
  saveRDS(list(time = time, diff.num = diff.num, time.mean = time.mean), 
          d.file("purge_unique_company.rds", exists = FALSE))
}

