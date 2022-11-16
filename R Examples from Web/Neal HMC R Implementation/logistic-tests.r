# TESTS WITH LOGISTIC MODEL.
#
# Radford M. Neal, 2012.

source("mcmc.r")
source("utilities.r")
source("singlevar.r")
source("metropolis.r")
source("basic_hmc.r")
source("pupdate.r")
source("logistic.r")

# Specify the data, which will be conditioned on using "fix".

X <- matrix(c(1,3,7,2,4,9,2,6,
              7,8,2,9,8,0,2,8), 8, 2)
y <-        c(0,0,1,1,0,1,1,1)

colnames(X) <- c("aa","bb")

# Function to look at results.

logistic.look <- function (r) 
{ q <- as.matrix(as.data.frame(r$q))
  print (rbind(mean=colMeans(q),sdev=sd(q)))
  cat("\n")
}

# Sample from the prior (with X fixed as above), the easy way, and the hard way.
# Comparing these results is a check on correctness (see John Geweke's paper,
# "Getting it right: Joint distribution tests of posterior simulators").

cat("prior the easy way:\n"); logistic.look(
rje <- mcmc (lpr_logistic, 
         list (alpha=0, beta=numeric(ncol(X)), log.sigma=0, y=numeric(nrow(X))),
         fix = list (X=X),
         iterations = 5000,
         list (sample_logistic_prior, level=0)))

cat("prior the hard way:\n"); logistic.look(
rjh <- mcmc (lpr_logistic, 
         list (alpha=0, beta=numeric(ncol(X)), log.sigma=0, y=numeric(nrow(X))),
         fix = list (X=X),
         iterations = 40000,
         list (sample_logistic_prior, level=2),
         c("alpha","beta","log.sigma"),
         list (singlevar, list (metropolis, step=0.5, rep=4)),
         list (metropolis, step=0.25, rep=10),
         list (metropolis, step=0.05, rep=10)))

# Sample from the posterior given data (X and y) above.

init <- list (alpha=0, beta=numeric(ncol(X)), log.sigma=0)

iters <- 5000

cat("A:\n"); logistic.look(
r <-   mcmc (lpr_logistic, init, iters, fix = list(X=X,y=y),
             "alpha", list (metropolis, step=0.2, rep=4),
             "beta",  list (metropolis, step=0.2, rep=4),
             "log.sigma",   list (metropolis, step=0.2, rep=4)))

cat("B:\n"); logistic.look(
ri <-  mcmc (ilpr_logistic, init, iters, fix = list(X=X,y=y),
             "alpha", list (metropolis, step=0.2, rep=4),
             "beta",  list (metropolis, step=0.2, rep=4),
             "log.sigma",   list (metropolis, step=0.2, rep=4)))

cat("C:\n"); logistic.look(
rs <-  mcmc (lpr_logistic, init, iters, fix = list(X=X,y=y),
             list (singlevar, list (metropolis, step=0.2, rep=4))))

cat("D:\n"); logistic.look(
rsi <- mcmc (ilpr_logistic, init, iters, fix = list(X=X,y=y),
             list (singlevar, list (metropolis, step=0.2, rep=4))))

cat("E:\n"); logistic.look(
ri1 <- mcmc (ilpr_logistic, init, iters, fix = list(X=X,y=y),
             "alpha", list (metropolis, step=0.2, rep=4),
             "beta",1,list (metropolis, step=0.2, rep=4),
             "beta",2,list (metropolis, step=0.2, rep=4),
             "log.sigma",   list (metropolis, step=0.2, rep=4)))

cat("F:\n"); logistic.look(
rsb <- mcmc (ilpr_logistic, init, iters, fix = list(X=X,y=y),
             "alpha", list (metropolis, step=0.2, rep=4),
             "beta",  list (singlevar, list (metropolis, step=0.2, rep=4)),
             "log.sigma",   list (metropolis, step=0.2, rep=4)))

cat("G:\n"); logistic.look(
rb0 <- mcmc (lpr_logistic, init, iters, fix = list(X=X,y=y), rec=c("acc","step"),
            "alpha", "beta", list (basic_hmc, nsteps=20, rand.step=0.3,
                               step = 0.2),
            "alpha", "beta", list(metropolis, step=0.1),
            "log.sigma", list (metropolis, step=0.2, rep=4)))

cat("H:\n"); logistic.look(
rb1 <- mcmc (lpr_logistic, init, iters, fix = list(X=X,y=y), rec=c("acc","step"),
            "alpha", "beta", list (basic_hmc, nsteps=20, rand.step=0.3,
                               step = function (s) 0.2*exp(min(0,s$log.sigma))),
            "log.sigma", list (metropolis, step=0.2, rep=4)))

cat("I:\n"); logistic.look(
rb2 <- mcmc (lpr_logistic, init, iters, fix = list(X=X,y=y), rec=c("acc","step"),
            "alpha", "beta", list (basic_hmc, nsteps=20, rand.step=0.3,
                               step = function (s) list (alpha=0.2, 
                                beta=rep(0.2*exp(min(0,s$log.sigma)),ncol(X)))),
            "log.sigma", list (metropolis, step=0.2, rep=4)))
