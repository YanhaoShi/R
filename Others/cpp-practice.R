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
}




cppFunction('List initAlgo1(NumericVector mu, NumericVector lB, NumericVector uB ){
    int n = mu.size();
  NumericVector m = clone(mu).sort();
  return List::create( mu[3], m(49), mu.sort(), std::max(1.0, m(49)), match(m, mu),
                     match(mu, m),n);
  
            }')

initAlgo1(assets50$mu, assets50$lB, assets50$uB)




cppFunction('List order_(NumericVector x) {
            IntegerVector v=x.size();
             NumericVector sorted = clone(x).sort();
            int k = 0;
            for(int i=0; i<x.size(); i++){
             
             v[i] = i;
            }
            int i2=1;
            while( i2 < 3){ k++; i2++; }
           
            return List::create(v,k);
            }')
order_(c(2,4,1,3,3))
cppFunction('IntegerVector order2(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
            Rf_warning("There are duplicates in x");
            }
            NumericVector sorted = clone(x).sort();
            return match(sorted, x);
            }')




getMatrices = function(mean, covar, solution_set, f){
  # Slice covarF,covarFB,covarB,meanF,meanB,wF,wB
  covarF <- covar[f,f]
  meanF <- mean[f]
  b <- seq_along(mean)[-f]
  covarFB <- covar[f,b]
  wB <- solution_set[,ncol(solution_set)][b] 
  return(list(covarF=covarF,covarFB=covarFB,meanF=meanF,wB=wB))
}
evalCpp("NumericMatrix(Dimension(2, 3))")
t <- getMatrices(s$mu, s$covar, s$solution_set, s$f)
###########################################################################
cppFunction('List get(NumericVector mu, NumericVector covar, 
            NumericMatrix solution_set, NumericVector f){
            int n = f.size();
            NumericMatrix covarF(n);
            NumericVector meanF(n);
            NumericMatrix covarF2(n);
//NumericMatrix::Sub covarF2 = covar(Range(0, n - 1), Range(0, n - 1));
for(int i=0; i < n; i++){
 meanF(i) = mu(f(i)-1);
for(int j =0; j <n; j++){

covarF(i,j) = covar(f(i)-1, f(j)-1);
}
}
return List::create(covarF, meanF, covarF2);
            }')

assets <- GetAssets(1:15, assets1457)
r <-  CLA$M4(assets$mu, assets$covar, assets$lB, assets$uB)
t <- getMatrices(assets$mu, assets$covar, r$weights_set_purge, c(1,2,4))
t2<- get(assets$mu, assets$covar, r$weights_set_purge, c(1,2,4))

all.equal(t2[[1]], t$covarF )
all.equal(t2[[2]], t$muF )



s <- validate(10,20)
microbenchmark(getMatrices(s$mu, s$covar, s$solution_set, c(1,2,4)),
               get(s$mu, s$covar, s$solution_set, c(1,2,4)))

####################################################################################





