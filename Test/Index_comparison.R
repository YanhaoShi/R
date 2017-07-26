analysis.company <- readRDS(d.file("analysis_company.rds", exists = FALSE))
analysis.sp <- readRDS(d.file("analysis_SP.rds", exists = FALSE))
analysis.nasdaq <- readRDS(d.file("analysis_NASDAQ.rds", exists = FALSE))
##draft
n.max <- max(1457, ncol(SP500), ncol(NASDAQ))
l.max = max(ncol(analysis.company$nweights), ncol(analysis.sp$nweights),
        ncol(analysis.nasdaq$nweights))
label.max <- exp_ind(1.2, 50, n.max)$ind[1:l.max]

# compare_nweights
nw.mean <- list(company = apply(analysis.company$nweights, 2, mean),
                sp = apply(analysis.sp$nweights, 2, mean),
                nasdaq = apply(analysis.nasdaq$nweights, 2, mean))
plot(nw.mean$company, xlab = "Number of Assets", ylab = "Number", 
        main = "Number of Weights Sets", xaxt = "n", type = "o", pch = 16, 
     xlim = c(0, l.max), ylim = range(nw.mean), col = "blue")
lines(nw.mean$sp, type = "o", pch = 16, col = "green")
lines(nw.mean$nasdaq, type = "o", pch = 16, col = "orange")
axis(1, at = seq(l.max), labels = label.max )
legend("topleft", legend = c("company assets", "SP500", "NASDAQ"),
       col = c("blue", "green", "orange"), lwd = 1, pch = 16)


# compare_kappa
kappa.mean <- list(company = apply(analysis.company$covF.kappa, 2, mean),
                     sp = apply(analysis.sp$covF.kappa, 2, mean),
                     nasdaq = apply(analysis.nasdaq$covF.kappa, 2, mean))
plot(kappa.mean$company, xlab = "Number of Assets", ylab = "kappa", 
     main = "Mean of max.kappa in covF", xaxt = "n", type = "o", pch = 16, 
     xlim = c(0, l.max), ylim = range(kappa.mean), col = "blue")
lines(kappa.mean$sp, type = "o", pch = 16, col = "green")
lines(kappa.mean$nasdaq, type = "o", pch = 16, col = "orange")
axis(1, at = seq(l.max), labels = label.max )
legend("topleft", legend = c("company assets", "SP500", "NASDAQ"),
       col = c("blue", "green", "orange"), lwd = 1, pch = 16)

# compare_rcond
rcond.mean <- list(company = apply(analysis.company$covF.rcond, 2, mean),
                     sp = apply(analysis.sp$covF.rcond, 2, mean),
                     nasdaq = apply(analysis.nasdaq$covF.rcond, 2, mean))
plot(rcond.mean$company, xlab = "Number of Assets", ylab = "rcond", 
     main = "Mean of min.cond in covF", xaxt = "n", type = "o", pch = 16, 
     xlim = c(0, l.max), ylim = range(rcond.mean), col = "blue")
lines(rcond.mean$sp, type = "o", pch = 16, col = "green")
lines(rcond.mean$nasdaq, type = "o", pch = 16, col = "orange")
axis(1, at = seq(l.max), labels = label.max )
legend("topright", legend = c("company assets", "SP500", "NASDAQ"),
       col = c("blue", "green", "orange"), lwd = 1, pch = 16)


# t_compare, mean
l.company <- 11
l.nasdaq <- 16
t_mean <- sapply(seq(l.sp), function(x) apply(analysis.sp$t_compare[[x]], 1, mean))
t.list <- list()
t.list$company <- sapply(seq(l.company),
                            function(x) apply(analysis.company$t_compare[[x]], 1, mean))
t.list$sp = sapply(seq(l.sp), function(x) apply(analysis.sp$t_compare[[x]], 1, mean))
t.list$nasdaq = sapply(seq(l.nasdaq),
                           function(x) apply(analysis.nasdaq$t_compare[[x]], 1, mean))

# compare_t_cla
plot(t.list$company[1, ], xlim = c(0, l.max), type = "o", 
     pch =16, col = "blue", xaxt = "n",
     ylim = range(t.list$company[1, ], t.list$sp[1, ], t.list$nasdaq[1, ]),
     xlab = "Number of Assets", ylab = "milliseconds", main = "CLA Time Comparison")
lines(t.list$sp[1, ], type = "o", pch = 16, col = "green")
lines(t.list$nasdaq[1, ], type = "o", pch = 16, col = "orange")
axis(1, at = seq(l.max), labels = label.max )
legend("topleft", legend = c("company assets", "SP500", "NASDAQ"),
       col = c("blue", "green", "orange"), lwd = 1, pch = 16)

# compare_t_qp
plot(t.list$company[2, ], xlim = c(0, l.max), type = "o", 
     pch =16, col = "blue", xaxt = "n",
     ylim = range(t.list$company[2, ], t.list$sp[2, ], t.list$nasdaq[2, ]),
     xlab = "Number of Assets", ylab = "milliseconds", main = "QP Time Comparison")
lines(t.list$sp[2, ], type = "o", pch = 16, col = "green")
lines(t.list$nasdaq[2, ], type = "o", pch = 16, col = "orange")
axis(1, at = seq(l.max), labels = label.max )
legend("topleft", legend = c("company assets", "SP500", "NASDAQ"),
       col = c("blue", "green", "orange"), lwd = 1, pch = 16)


