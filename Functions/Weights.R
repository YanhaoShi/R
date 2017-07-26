# "true" free indices 
freeTrue <- function(weights_set){
  lapply(1:ncol(weights_set), function(x)(
    which( !(abs(weights_set[,x]-1e-8)< 1e-10| abs(weights_set[,x]- 0.075)<1e-10))
  ))
}

change_f <- function(free_indices){
  stopifnot((nf <- length(free_indices)) >= 2)
  ii <- seq_len(nf-1)
  lapply(ii, function(x) c(setdiff(free_indices[[x]], free_indices[[x+1]]), 
                           setdiff(free_indices[[x+1]],free_indices[[x]])) )
  
  f2b <- c(NA, sapply(1:(length(free_indices)-1), function(x)
    {d <- setdiff(free_indices[[x]], free_indices[[x+1]])
  ifelse(length(d)!=0, d, 0)}  ))
  b2f <- c(NA, sapply(1:(length(free_indices)-1), function(x) 
    {d <- setdiff(free_indices[[x+1]], free_indices[[x]])
  ifelse(length(d)!=0, d, 0)} ))
  cbind(f2b, b2f)
}

# remove replicated columns
purgeUnique <- function(weights_set){
  n <- ncol(weights_set)
  ind.unique <- which(sapply(1:(n-1), function(x) 
    !isTRUE(all.equal(weights_set[,x],weights_set[,x+1]))))
  c(ind.unique, n)
}

addW <- function(weights_set, mu){
  n <- ncol(weights_set)
  f_ind <- Wf(weights_set)
  dif <- apply(change_f(f_ind) !=0, 2, sum)
  w_add <-c ()
  for(j in which(dif==2)){
    f2b <- setdiff(f_ind[[j]], f_ind[[j+1]]) #row
    w_in <- weights_set[,j]
    w_in[f2b] <- weights_set[f2b, j+1] # set free weight to the boundary
    ind <- intersect(f_ind[[j]], f_ind[[j+1]])
    w_in[ind] <-w_in[ind]/(sum(w_in[ind])) * (1-sum(w_in[-ind]))
    w_add <- cbind(w_add,w_in)
  }
  weights_set <- cbind(weights_set, w_add)
  weights_set[, order(t(weights_set) %*% mu, decreasing = TRUE)]
       
}

plotW <- function(weights_set, mu, covar){
  ind_p <- purgeChull(weights_set, mu, covar)
  pch <- rep(1, ncol(weights_set))
  pch[ind_p] <- 16
  col <- rep("blue", ncol(weights_set))
  col[ind_p] <- "blue"
  if(!is.null(colnames(weights_set))){  
    ind_in <- colnames(weights_set) == "w_in"
    col[ind_in] <- "red"}
  f.label <- change_f(Wf(weights_set))
  
  ms.w <- MS(weights_set, mu = mu, covar = covar)
  plot(ms.w[,"Sig2"], ms.w[,"Mu"], pch = pch, col = col)
  lines(ms.w[ind_p,"Sig2"], ms.w[ind_p,"Mu"], pch = 1, col = adjustcolor("blue", 0.5))
  ind1 <- f.label[1, ]!= 0
  ind2 <- f.label[2, ]!= 0
  text(ms.w[ind1, "Sig2"], ms.w[ind1, "Mu"], pos = 3, paste0("+",f.label[1,ind1]), col = 1)
  text(ms.w[ind2, "Sig2"], ms.w[ind2, "Mu"], pos = 1, paste0("-",f.label[2,ind2]), col = 1)
  legend("bottomright", legend = c("w", "w-in", "w-remove", "w-in-remove"), 
         col = c("blue", "red"), pch = c(16, 16, 1, 1), bty = "n")
}

.freeI2char <- function(f) paste0("(", paste(f, collapse=","), ")")
freeI2char <- function(free_indices) sapply(free_indices, .freeI2char)

Result_table <- function(result, mu, covar){
  options(digits=5)
  ms <- MS(result$weights_set, mu, covar)
  ftrue <- freeTrue(result$weights_set)
  data.frame(ms[,c("Mu", "Sig")], lambda = result$lambdas, 
             gamma = result$gammas, freeI = freeI2char(result$free_indices),
             ind = change_f(result$free_indices), 
             freeI.true = freeI2char(ftrue),
             ind.change = change_f(ftrue))
}

# add f_change!
