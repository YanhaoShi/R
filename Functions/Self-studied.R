# Pointers
f <- vector('list',3)

for(i in 1:length(f)) {
  f[[i]] <- function(){i}
}
# i = 3, function(){3}
for(j in 1:length(f)){
  print(f[[j]]())
}
# print 3,3,3; but not 1,2,3.

# If we want each function to capture the current value of i in its closure
# and not share a common reference:
f <- vector('list',3)

for(i in 1:length(f)) {
  f[[i]] <- function(){i}
  e <- new.env()
  assign("i", i, envir = e)
  environment(f[[i]]) <- e
}

for(j in 1:length(f)){
  print(f[[j]]())
}

############################################################################
m <- 3
f <- function(x){ 
  m <- 5
  f.env <- new.env(); assign("m", m, envir = f.env)
  m <- 2
  c(x + m, x + globalenv()$m, f.env$m)
} 

f(0) 


