#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/


// [[Rcpp::export]]
RObject ex_Matrix_rowcol(NumericMatrix x, int ir, int ic)
{
  NumericMatrix::Row    r = x.row(ir - 1);
  NumericMatrix::Column c = x.column(ic - 1);
  // Row and Column act as vectors
  r[0] = c[0] = 0;
  // They can also be assigned to vectors
  NumericVector rvec = r;
  NumericVector cvec = c;
  return List::create(Named("row") = rvec, Named("col") = cvec);
}

/*** R

x = matrix(as.numeric(1:9), 3, 3)
  ex_Matrix_rowcol(x, 2, 2);
## $row
## [1] 0 5 8
##
## $col
## [1] 0 5 6
x
##      [,1] [,2] [,3]
## [1,]    1    0    7
## [2,]    0    5    8
## [3,]    3    6    9
  
  */

// [[Rcpp::export]]
RObject ex_Matrix_fromsub(NumericMatrix x, int n)
{
  NumericMatrix::Sub sub = x(Range(0, n - 1), Range(0, n - 1));
  NumericMatrix res(sub);
  return res;
}

/*** R

x = matrix(as.numeric(1:9), 3, 3)
  ex_Matrix_fromsub(x, 2);
##      [,1] [,2]
## [1,]    1    4
## [2,]    2    5

*/



//
List get(NumericVector mu, NumericMatrix covar, 
         NumericMatrix solution_set, NumericVector f){
  int n = f.size();
  NumericMatrix covarF(n);
  NumericVector meanF(n);
 // NumericMatrix covarF2(n);
 NumericMatrix::Sub sub = covar(Range(0, n - 1), Range(0, n - 1));
 NumericMatrix covarF2(sub);
  for(int i=0; i < n; i++){
    meanF(i) = mu(f(i)-1);
    for(int j =0; j <n; j++){
      
      covarF(i,j) = covar(f(i)-1, f(j)-1);
    }
  }
  return List::create(covarF, meanF, covarF2);
}
/*** R
 # assets <- GetAssets(1:15, assets1457)
#  r <-  CLA$M4(assets$mu, assets$covar, assets$lB, assets$uB)
#  t <- getMatrices(assets$mu, assets$covar, r$weights_set_purge, c(1,2,4))
 # t2<- get(assets$mu, assets$covar, r$weights_set_purge, c(1,2,4))
  
#  all.equal(t2[[1]], t$covarF )
 # all.equal(t2[[2]], t$muF )
  
*/

