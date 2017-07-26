library(Rcpp)
library(microbenchmark)
cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
            return sum;
            }')


add(1,2,3)
sum(1:3)
add.r <- function(x,y,z) x+y+z
microbenchmark(add(1,2,3), sum(c(1,2,3)), add.r(1,2,3))


cppFunction('int one(){
            return 1;
            }')
one()

# declare the type of each input in the same way we declare the 
# type of the output. 

sumR <- function(x) {
  total <- 0
  for (i in seq_along(x)) {
    total <- total + x[i]
  }
  total
}
#  the cost of loops is much lower in C++.
#  In C++, vector indices start at 0
cppFunction('double sumC(NumericVector x) {
  int n = x.size();
            double total = 0;
            for(int i = 0; i < n; ++i) {
            total += x[i];
            }
            return total;
            }')

x <- runif(1e3)
microbenchmark(
  sum(x),
  sumC(x),
  sumR(x)
)
all.equal(sumR(x), sumC(x))

# Euclidean distance
cppFunction('NumericVector pdistC(double x, NumericVector ys) {
  int n = ys.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
           out[i] = sqrt(pow(ys[i] - x, 2.0));
      }
     return out;
 }')


# NumericMatrix, IntegerMatrix, CharacterMatrix, and LogicalMatrix
# int, double, bool, string
cppFunction('NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
            NumericVector out(nrow);
            
            for (int i = 0; i < nrow; i++) {
            double total = 0;
            for (int j = 0; j < ncol; j++) {
            total += x(i, j);
            }
            out[i] = total;
            }
            return out;
            }')
set.seed(1014)
x <- matrix(sample(100), 10)
rowSums(x)
#subset a matrix with (), not [].
#Use .nrow() and .ncol() get the dimensions of a matrix.
cppFunction('NumericVector attribs() {
  NumericVector out = NumericVector::create(1, 2, 3);
  
  out.names() = CharacterVector::create("a", "b", "c");
  out.attr("my-attr") = "my-value";
  out.attr("class") = "my-class";
  
  return out;
}')

#  ::create create R vector
attribs()

cppFunction('double mpe(List mod) {
  if (!mod.inherits("lm")) stop("Input must be a linear model");
            
            NumericVector resid = as<NumericVector>(mod["residuals"]);
            NumericVector fitted = as<NumericVector>(mod["fitted.values"]);
            
            int n = resid.size();
            double err = 0;
            for(int i = 0; i < n; ++i) {
            err += resid[i] / (fitted[i] + resid[i]);
            }
            return err / n;
            }')

mod <- lm(mpg ~ wt, data = mtcars)
mpe(mod)

cppFunction('List scalar_missings() {
  int int_s = NA_INTEGER;
            String chr_s = NA_STRING;
            bool lgl_s = NA_LOGICAL;
            double num_s = NA_REAL;
            
            return List::create(int_s, chr_s, lgl_s, num_s);
            }')

scalar_missings()
str(scalar_missings)
evalCpp("NAN > 1")

# Vectors
# NA_REAL, NA_INTEGER, NA_LOGICAL, NA_STRING:
cppFunction('List missing_sampler() {
  return List::create(
            NumericVector::create(NA_REAL),
            IntegerVector::create(NA_INTEGER),
            LogicalVector::create(NA_LOGICAL),
            CharacterVector::create(NA_STRING));
            }')
missing_sampler()
str(missing_sampler)

cppFunction('LogicalVector is_naC2(NumericVector x) {
  return is_na(x);
}')
is_naC2(c(1,2,3,NA))

cppFunction('NumericVector pdistC2(double x, NumericVector ys) {
  return sqrt(pow((x - ys), 2));
            }')
pdistC2(1, 1:3)
pdistC(1, 1:3)

cppFunction('bool any_naC(NumericVector x) {
  return is_true(any(is_na(x)));
            }')

cppFunction('double sum3(NumericVector x) {
  double total = 0;
            
            NumericVector::iterator it;
            for(it = x.begin(); it != x.end(); ++it) {
            total += *it;
            }
            return total;
            }')

cppFunction('double sum4(NumericVector x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
            }')
sum4(c(1,2,4))

cppFunction('IntegerVector findInterval2(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());
            
            NumericVector::iterator it, pos;
            IntegerVector::iterator out_it;
            
            for(it = x.begin(), out_it = out.begin(); it != x.end(); 
            ++it, ++out_it) {
            pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
            *out_it = std::distance(breaks.begin(), pos);
            }
            
            return out;
            }')





