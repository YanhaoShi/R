#include<iostream>
using namespace std;
int main(){ 
 
  int mu, lB, uB;
  std::cin>> mu >> lB >> uB;
 /* NumericVector w(n);
  index.new <- sort (mu,decreasing = T)
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
  w[index.new] <- w.new               
  i<-index.new[i.new] 
  return List::create(index=i,weights=w)    */    
  std::cout << mu << std::endl;
  return 0;
}


