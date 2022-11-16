# TESTS ON MULTIVARIATE NORMAL DISTRIBUTION WITH BOUNDS.
#
# Radford M. Neal, 2012.

source("mcmc.r")
source("utilities.r")
source("metropolis.r")
source("singlevar.r")
source("bounded_hmc.r")
source("pupdate.r")
source("mvn.r")


# --------------------------------------------------------------------

# FUNCTIONS FOR SAMPLING A BOUNDED MULTIVARIATE NORMAL IN VARIOUS WAYS.  


# SAMPLE WITH BOUNDS USING JOINT METROPOLIS UPDATES.

cmet_mvn <- function (iterations, step, lower=NULL, upper=NULL, ...)
{
  mcmc (structure(lpr_mvn,lower=lower,upper=upper), c(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=20),
        ...
  )
}


# SAMPLE WITH BOUNDS USING HAMILTONIAN MONTE CARLO WITH LONG TRAJECTORY.

chmc_mvn <- function (iterations, step, lower=NULL, upper=NULL, ...)
{
  mcmc (structure(lpr_mvn,lower=lower,upper=upper), c(a=0, b=0), iterations, 
        list(bounded_hmc, step=step, nsteps=20),
        ave=c("apr","reflections"), ...
  )
}


# SAMPLE WITH BOUNDS ONE VARIABLE AT A TIME WITH ONE LEAPFROG STEP (LANGEVIN).
# Doesn't use singlevar.

clmc_mvn <- function (iterations, step, lower=NULL, upper=NULL, ...)
{
  mcmc (structure(lpr_mvn,lower=lower,upper=upper), c(a=0, b=0), iterations, 
        1, list(bounded_hmc, step=step),
        2, list(bounded_hmc, step=step),
        ...
  )
}


# SAMPLE WITH BOUNDS ONE VARIABLE AT A TIME WITH ONE LEAPFROG STEP (LANGEVIN).
# Uses singlevar.

sclmc_mvn <- function (iterations, step, lower=NULL, upper=NULL, ...)
{
  mcmc (structure(lpr_mvn,lower=lower,upper=upper), c(a=0, b=0), iterations, 
        list (singlevar, list(bounded_hmc, step=step)),
        ...
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

inv.cov <- solve (matrix (c (0.5, 0.95, 0.95, 2.0), 2, 2))

formals(lpr_mvn)$inv.cov <- inv.cov

lower <- list (NULL, -0.5, NULL, -1.0, c(0,-0.5), -0.2,    -0.2)
upper <- list (NULL, NULL,  1.5,  0.8, c(1,1.5),  c(0,Inf), 0)

it <- 2000

for (i in seq_along(lower))
{ 
  cat("\n----",i,"\n\nBounds:\n")

  print (rbind(lower=lower[[i]], upper=upper[[i]]))

  cat("\ncmet:\n");
  mvn.look(
   assign(paste("rcmet",i,sep=""), cmet_mvn (it, 1.5, lower[[i]], upper[[i]])))

  cat("\nchmc:\n");
  mvn.look(
   assign(paste("rchmc",i,sep=""), chmc_mvn (it, 0.2, lower[[i]], upper[[i]])))

  cat("\nclmc:\n");
  mvn.look(
   assign(paste("rclmc",i,sep=""), clmc_mvn (20*it,0.2,lower[[i]], upper[[i]])))

  cat("\nsclmc:\n");
  mvn.look(
   assign(paste("rsclmc",i,sep=""),sclmc_mvn(20*it,0.2,lower[[i]], upper[[i]])))
}
