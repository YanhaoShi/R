## Set up #####
if(!file.exists("D:/Master Thesis/Codes/R/init.R")) stop("R directory does *not* contain init.R")
source("D:/Master Thesis/Codes/R/init.R")
options(error = NULL)
#######################################################################################
## 1) Packages Required
pa <- c("microbenchmark", "Matrix", "compiler", "png", "fGarch",
        "FRAPO", "ghyp", "Rcpp", "quadprog", "matrixcalc", "xtable")
for (package in pa) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


## 4) check index data in FRAPO package
if(!file.exists(d.file("indexdf", exists = FALSE))){
  index <- list(ESCBFX = ESCBFX, EuroStoxx50 = EuroStoxx50, FTSE100 = FTSE100, 
                INDTRACK1 = INDTRACK1, INDTRACK2 = INDTRACK2, 
                INDTRACK3 = INDTRACK3, INDTRACK4 = INDTRACK4, 
                INDTRACK5 = INDTRACK5, INDTRACK6 = INDTRACK6, 
                MIBTEL = MIBTEL, MultiAsset = MultiAsset, NASDAQ = NASDAQ, 
                SP500 = SP500, StockIndex = StockIndex, 
                StockIndexAdj = StockIndexAdj, StockIndexAdjD = StockIndexAdjD)
  index.df <- data.frame(assets = sapply(index, ncol),
                         observations = sapply(index, nrow),
                         start = sapply(index, function(x) 
                           ifelse(is.null(row.names(x)), 1, row.names(x)[1])),
                         end = sapply(index, function(x) 
                           ifelse(is.null(row.names(x)), nrow(x), row.names(x)[nrow(x)])))
  saveRDS(index.df, d.file("indexdf.rds", exists = FALSE))
}

exp_ind <- function(x, start, end){
  a <- ceiling(log(start)/log(x) )
  b <- floor(log(end + 1)/log(x))
  list(power = a:b, ind = round(x^(a:b)))
}

## 5) Source functions from environment
funEnv <- function(...) {  ## by Martin Maechler
  e <- list2env(list(...))
  for(n in names(e)) environment(e[[n]]) <- e
  e
}

#6) source all files in Function directory

funfile <- paste0("Functions/", list.files("Functions"))
source("Functions/Data_prepare.R") # check if all rds exists
for(filename in seq_along(funfile)[-c(2)]){ # run all other function files
  source(funfile[filename]) 
}

verfile <- grep("Ver", list.files("Functions"), value = TRUE) 
# file names contain "Ver"
verind <- sort(gsub(".R", "", gsub("Ver", "", verfile)))
ver.env <- as.list(verind)
names(ver.env) <- as.list(paste0("Env", verind))
# ordered version index 

info <- sapply(verind, function(vi) get(paste0("Env", vi))$info[[1]])
Ver.info <- data.frame(info = info)

trans.company <- assets1457

result20.list <- lapply(ver.env, function(ind) 
  get(paste0("Env", ind))$cla.solve(GetIndex(1:20, trans.company)))
result50.list <- lapply(ver.env, function(ind) 
  get(paste0("Env", ind))$cla.solve(GetIndex(1:50, trans.company)))


