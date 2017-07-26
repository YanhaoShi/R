################################## Set up ##################################

## Add your user name and working directory path to the statement below
## Original version by Martin Maechler

si    <- Sys.info()
User  <- si[["user"]]
compi <- si[["nodename"]]

d.dirs <- switch(User,
                 'maechler' = c(R = if(compi == "nb-mm4") "~/Dropbox_VS/R_scripts"
                                else "R_dir_at_sfs__TODO",
                                data = "~/Betreute-Arbeiten/data"),
                 'user' = c(R = "D:/Master Thesis/Codes/R", 
                            data = "D:/Master Thesis/Codes/Data"),
                 
                 'shiy' = c(R = "C:/Users/shiy/polybox/Codes/R", 
                            data = "C:/Users/shiy/polybox/Codes/Data"),
                 stop("add your user name ", User, " and directory paths to 'init.R'"))

if(!dir.exists(wd <- d.dirs[["R"]]))    stop("R directory ",    wd, " not valid")
if(!dir.exists(dd <- d.dirs[["data"]])) stop("data directory ", dd, " not valid")

rm(si, User, compi, wd, dd) ## but not 'd.dirs' !!

d.file <- function(file, exists = TRUE) {
    f <- file.path(d.dirs[["data"]], file)
    if(exists && !file.exists(f))
      stop("file ", f, " does not exist (possibly relatively to working dir).")
    f
}

setwd(d.dirs[["R"]]) ## working directory 'R_scripts'


