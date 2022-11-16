# TESTS USING OF WINDOWED HMC USING MULTIVARIATE NORMAL DISTRIBUTION.
#
# Radford M. Neal, 2012.

source("mcmc.r")
source("utilities.r")
source("basic_hmc.r")
source("windowed_hmc.r")
source("pupdate.r")
source("mvn.r")


# --------------------------------------------------------------------

# FUNCTIONS FOR SAMPLING A MULTIVARIATE NORMAL WITH BASIC AND WINDOWED HMC.


# SAMPLE USING BASIC HAMILTONIAN MONTE CARLO.

bhmc_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations,
        list(basic_hmc, step=step, nsteps=20),
        ...
  )
}


# SAMPLE USING BASIC HAMILTONIAN MONTE CARLO, BUT USING WINDOWED HMC FUNCTION.

wbhmc_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations,
        list(windowed_hmc, step=step, nsteps=20),
        ...
  )
}


# SAMPLE USING WINDOWED HAMILTONIAN MONTE CARLO.

whmc_mvn <- function (iterations, step, nsteps, window, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations,
        list(windowed_hmc, step=step, nsteps=nsteps, window=window),
        rec=c("acc","ipos","move"), ...
  )
}

# --------------------------------------------------------------------

# FUNCTION TO LOOK AT RESULTS OF THE TESTS BELOW.

mvn.look <- function (r)
{ if (is.list(r$q))
  { q1 <- r$q[[1]]
    q2 <- r$q[[2]]
  }
  else
  { q1 <- r$q[,1]
    q2 <- r$q[,2]
  }
  cat("m1:",round(mean(q1),4)," m2:",round(mean(q2),3),
      " v1:",round(var(q1),4)," v2:",round(var(q2),3),
      " corr:",round(cor(q1,q2),3),
      " act+:",round (2*sum(acf(q1+q2,plot=FALSE,lag.max=100)$acf)-1, 1),
      " act-:",round (2*sum(acf(q1-q2,plot=FALSE,lag.max=100)$acf)-1, 1),
      "\n")
  plot(q1,q2,pch=".")
}


# --------------------------------------------------------------------

# SCRIPT FOR RUNNING THE TESTS.

inv.cov <- solve (matrix (c (0.5, 0.97, 0.97, 2.0), 2, 2))

formals(lpr_mvn)$inv.cov <- inv.cov

# Tests with just four iterations, for comparing results.  Should be
# the same for standard HMC using basic_hmc and windowed_hmc.

cat("SHORT TESTS\n")

cat("\nb:\n");   print ((r4.b <-  bhmc_mvn (4, 0.22)) $ q)
cat("\nwb:\n");  print ((r4.wb <- wbhmc_mvn (4, 0.22)) $ q)

# Tests with 10000 iterations, for comparing with the correct distribution.

cat("\nLONG TESTS\n\n")

cat("b1:\n"); mvn.look(r <- rb1 <- bhmc_mvn (10000, 0.22))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("b2:\n"); mvn.look(r <- rb2 <- bhmc_mvn (10000, 0.23))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("b3:\n"); mvn.look(r <- rb3 <- bhmc_mvn (10000, 0.24))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("wb:\n"); mvn.look(r <- rwb <- wbhmc_mvn(10000, 0.22))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")

cat("w0:\n"); mvn.look(r <- rw0 <- whmc_mvn (10000, 0.22, nsteps=20,window=2))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("w1:\n"); mvn.look(r <- rw1 <- whmc_mvn (10000, 0.23, nsteps=20,window=3))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("w2:\n"); mvn.look(r <- rw2 <- whmc_mvn (10000, 0.23, nsteps=20, 
                           window=structure(3,weights=c(0.5,1,0.5))))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")

cat("w3:\n"); mvn.look(r <- rw3 <- whmc_mvn (10000, 0.23, nsteps=20,window=10))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("w4:\n"); mvn.look(r <- rw4 <- whmc_mvn (10000, 0.23, nsteps=20, 
                           window=structure(10,weights=10:1)))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("w5:\n"); mvn.look(r <- rw5 <- whmc_mvn (10000, 0.23, nsteps=20,window=21))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("w6:\n"); mvn.look(r <- rw6 <- whmc_mvn (10000, 0.23, nsteps=20, 
                                window=structure(21,weights=1:21)))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
