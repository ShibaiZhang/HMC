# TESTS USING MULTIVARIATE NORMAL DISTRIBUTION.
#
# Radford M. Neal, 2012.

source("mcmc.r")
source("utilities.r")
source("singlevar.r")
source("metropolis.r")
source("basic_hmc.r")
source("pupdate.r")
source("gibbs.r")
source("mvn.r")


# --------------------------------------------------------------------

# LPR FUNCTIONS FOR TEST PURPOSES.  The lpr_mvn and ilpr_mvn functions
# defined in mvn.r are used along with the two versions below, which 
# take a state that is a list, and which are meant solely for testing 
# the mcmc function.


# LIST VERSION OF LOG DENSITY FUNCTION. 

llpr_mvn <- function (value,grad=FALSE,inv.cov)
{
  v <- c (value$a, value$b)

  g <- - inv.cov %*% v
  r <- t(v) %*% g / 2
  if (grad)
  { g <- as.vector(g)
    attr(r,"grad") <- 
      list (a=g[1:length(value$a)], b=g[(length(value$a)+1):length(v)])
  }

  r
}


# INCREMENTAL LIST VERSION OF LOG DENSITY FUNCTION.  

illpr_mvn <- function (value,grad=FALSE,ch.value,ch.elem,ch.pos,lpr.value,
                       inv.cov)
{
  v <- c (value$a, value$b)

  if (missing(ch.value) || is.null(attr(lpr.value,"g")))
  { g <- - inv.cov %*% v
  }
  else 
  { if (!(ch.elem %in% c("a","b")))
    { stop("Bad ch.elem argument")
    }
    if (missing(ch.pos))
    { if (ch.elem=="a") 
      { ch.pos <- 1:length(value$a)
      }
      else
      { ch.pos <- (length(value$a)+1):length(v)
      }
    }
    else
    { if (ch.elem=="b")
      { ch.pos <- ch.pos + length(value$a)
      }
    }
    g <- attr(lpr.value,"g") - 
           inv.cov[,ch.pos,drop=FALSE] %*% (ch.value-v[ch.pos])
    v[ch.pos] <- ch.value
  }

  r <- t(v) %*% g / 2

  g <- as.vector(g)
  attr(r,"g") <- g
  if (grad) 
  { attr(r,"grad") <- 
      list (a=g[1:length(value$a)], b=g[(length(value$a)+1):length(v)])
  }

  r
}


# --------------------------------------------------------------------

# FUNCTIONS FOR SAMPLING A MULTIVARIATE NORMAL IN VARIOUS WAYS.  These
# are designed as tests for the mcmc function.


# SAMPLE FROM MULTIVARIATE NORMAL USING JOINT METROPOLIS UPDATES.
# Uses non-incremental lpr.

jmet_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=2),
        list (metropolis, step=step/10),
        ...
  )
}


# SAMPLE FROM MULTIVARIATE NORMAL USING JOINT METROPOLIS UPDATES.
# Uses list version of non-incremental lpr.

ljmet_mvn <- function (iterations, step, ...)
{
  mcmc (llpr_mvn, list(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=2),
        list (metropolis, step=lapply(step,function(x)x/10)),
        ...
  )
}


# SAMPLE FROM MULTIVARIATE NORMAL USING JOINT METROPOLIS UPDATES.
# Uses incremental lpr (non-incrementally).

ijmet_mvn <- function (iterations, step, ...)
{
  mcmc (ilpr_mvn, c(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=2),
        list (metropolis, step=step/10),
        ...
  )
}


# SAMPLE FROM MULTIVARIATE NORMAL USING JOINT METROPOLIS UPDATES.
# Uses list version of incremental lpr (non-incrementally).

iljmet_mvn <- function (iterations, step, ...)
{
  mcmc (illpr_mvn, list(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=2),
        list (metropolis, step=lapply(step,function(x)x/10)),
        ...
  )
}


# SAMPLE FROM MULTIVARIATE NORMAL USING SINGLE VARIABLE METROPOLIS UPDATES.
# Uses non-incremental lpr.

smet_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations, 
        list (singlevar, metropolis, step=step, rep=5),
        ...
  )
}


# SAMPLE FROM MULTIVARIATE NORMAL USING SINGLE VARIABLE METROPOLIS UPDATES.
# Uses incremental lpr.

ismet_mvn <- function (iterations, step, ...)
{
  mcmc (ilpr_mvn, c(a=0, b=0), iterations, 
        list (singlevar, metropolis, step=step, rep=5),
        ...
  )
}


# SAMPLE USING SINGLE VARIABLE UPDATES WITHOUT USING SINGLEVAR.
# Uses non-incremental lpr.

omet_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations,
        1, list(metropolis, step=step, rep=5),
        2, list(metropolis, step=step, rep=5),
        ...
  )
}


# SAMPLE USING SINGLE VARIABLE UPDATES WITHOUT USING SINGLEVAR.
# Uses list version of non-incremental lpr.

lomet_mvn <- function (iterations, step, ...)
{
  mcmc (llpr_mvn, list(a=0, b=0), iterations,
        "a", list(metropolis, step=step, rep=5),
        "b", list(metropolis, step=step, rep=5),
        ...
  )
}


# SAMPLE USING SINGLE VARIABLE UPDATES WITHOUT USING SINGLEVAR.
# Uses list version of incremental lpr.

ilomet_mvn <- function (iterations, step, ...)
{
  mcmc (illpr_mvn, list(a=0, b=0), iterations,
        "a", list(metropolis, step=step, rep=5),
        "b", list(metropolis, step=step, rep=5),
        ...
  )
}


# SAMPLE WITH BOUNDS USING JOINT METROPOLIS UPDATES.

cmet_mvn <- function (iterations, step, ...)
{
  mcmc (structure(lpr_mvn,lower=-1,upper=c(Inf,2)), c(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=10),
        ...
  )
}


# SAMPLE WITH BOUNDS USING JOINT METROPOLIS UPDATES.
# Uses list version.

lcmet_mvn <- function (iterations, step, ...)
{
  mcmc (structure(llpr_mvn,lower=-1,upper=c(Inf,2)), list(a=0, b=0), iterations,
        list (metropolis, step=step, rep=10),
        ...
  )
}


# SAMPLE USING BASIC HAMILTONIAN MONTE CARLO WITH LONG TRAJECTORY.

bhmc_mvn <- function (iterations, step, rand.step=0, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations,
        list(basic_hmc, step=step, nsteps=20, rand.step=rand.step),
        ...
  )
}


# SAMPLE USING BASIC HAMILTONIAN MONTE CARLO WITH LONG TRAJECTORY.
# Uses list version of lpr function.

lbhmc_mvn <- function (iterations, step, ...)
{
  mcmc (llpr_mvn, list(a=0, b=0), iterations,
        list(basic_hmc, step=step, nsteps=20),
        ...
  )
}


# SAMPLE USING BASIC HAMILTONIAN MONTE CARLO WITH LONG TRAJECTORY.
# Uses list version of incremental lpr function.

ilbhmc_mvn <- function (iterations, step, ...)
{
  mcmc (illpr_mvn, list(a=0, b=0), iterations,
        list(basic_hmc, step=step, nsteps=20),
        ...
  )
}


# SAMPLE USING BASIC HAMILTONIAN MONTE CARLO WITH LONG TRAJECTORY.
# Keep around momentum.

pbhmc_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations, initial.p = numeric(2),
        pupdate, list(basic_hmc, step=step, nsteps=20),
        ...
  )
}


# GET TRAJECTORIES FROM BASIC HAMILTONIAN UPDATES.

bhmc_traj_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations,
        list(basic_hmc, step=step, nsteps=20, return.traj=TRUE),
        rec = "traj.q", last="traj.H",
        ...
  )
}


# SAMPLE ONE VARIABLE AT A TIME USING HMC WITH ONE LEAPFROG STEP (LANGEVIN).

lhmc_mvn <- function (iterations, step, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations, rep=10,
        1, list(basic_hmc, step=step),
        2, list(basic_hmc, step=step),
        ...
  )
}


# SAMPLE ONE VARIABLE AT A TIME USING HMC WITH ONE LEAPFROG STEP (LANGEVIN).
# Uses list version of incremental lpr function.

illhmc_mvn <- function (iterations, step, ...)
{
  mcmc (illpr_mvn, list(a=0, b=0), iterations, rep=10,
        "a", list(basic_hmc, step=step),
        "b", list(basic_hmc, step=step),
        ...
  )
}


# SAMPLE ONE VARIABLE AT A TIME USING HMC WITH ONE LEAPFROG STEP (LANGEVIN).
# Uses singlevar.

slhmc_mvn <- function (iterations, step, ...)
{
  mcmc (llpr_mvn, list(a=0, b=0), iterations, rep=10,
        list (singlevar, list(basic_hmc, step=step)),
        ...
  )
}


# NORMAL GIBBS SAMPLING. Uses non-incremental lpr.

gibbs_mvn <- function (iterations, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations, 
        list (singlevar, gibbs_normal),
        ...
  )
}


# NORMAL GIBBS SAMPLING. Uses incremental lpr.

igibbs_mvn <- function (iterations, ...)
{
  mcmc (ilpr_mvn, c(a=0, b=0), iterations, 
        list (singlevar, gibbs_normal),
        ...
  )
}


# ADLER OVERRELAXATION. Uses non-incremental lpr.

over_mvn <- function (iterations, a, ...)
{
  mcmc (lpr_mvn, c(a=0, b=0), iterations, 
        list (singlevar, gibbs_normal, a=a),
        ...
  )
}


# ADLER OVERRELAXATION. Uses incremental lpr.

iover_mvn <- function (iterations, a, ...)
{
  mcmc (ilpr_mvn, c(a=0, b=0), iterations, 
        list (singlevar, gibbs_normal, a=a),
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
formals(ilpr_mvn)$inv.cov <- inv.cov
formals(llpr_mvn)$inv.cov <- inv.cov
formals(illpr_mvn)$inv.cov <- inv.cov

# Tests with just four iterations, for comparing results.  Should be
# the same (except for vector vs. list) within each group below.

cat("SHORT TESTS\n\n")

cat("\nj:\n");   print ((r4.j <-   jmet_mvn (4, c(0.2,0.4))) $ q)
cat("\nlj:\n");  print ((r4.lj <-  ljmet_mvn (4, list(a=0.2,b=0.4))) $ q)
cat("\nij:\n");  print ((r4.ij <-  ijmet_mvn (4, c(0.2,0.4))) $ q)
cat("\nilj:\n"); print ((r4.ilj <- iljmet_mvn (4, list(a=0.2,b=0.4))) $ q)

cat("\ns:\n");   print ((r4.s <-   smet_mvn (4, 0.5)) $ q)
cat("\nis:\n");  print ((r4.is <-  ismet_mvn (4, 0.5)) $ q)
cat("\no:\n");   print ((r4.o <-   omet_mvn (4, 0.5)) $ q)
cat("\nlo:\n");  print ((r4.lo <-  lomet_mvn (4, 0.5)) $ q)
cat("\nilo:\n"); print ((r4.ilo <- ilomet_mvn (4, 0.5)) $ q)

cat("\nc:\n");   print ((r4.c <-   cmet_mvn (4, 1.5)) $ q)
cat("\nlc:\n");  print ((r4.lc <-  cmet_mvn (4, 1.5)) $ q)

cat("\nb:\n");   print ((r4.b <-   bhmc_mvn (4, 0.2)) $ q)
cat("\nlb:\n");  print ((r4.lb <-  lbhmc_mvn (4, 0.2)) $ q)
cat("\nilb:\n"); print ((r4.ilb <- ilbhmc_mvn (4, 0.2)) $ q)
cat("\nbt:\n");  print ((r4.bt <-  bhmc_traj_mvn (4, 0.2)) $ q)
cat("\npb:\n");  print ((r4.pb <-  pbhmc_mvn (4, 0.2)) $ q)

cat("\nl:\n");   print ((r4.l <-   lhmc_mvn (4, 0.1)) $ q)
cat("\nill:\n"); print ((r4.ill <- illhmc_mvn (4, 0.1)) $ q)
cat("\nsl:\n");  print ((r4.sl <-  slhmc_mvn (4, 0.1)) $ q)

cat("\ng:\n");   print ((r4.g <-  gibbs_mvn (4)) $ q)
cat("\nig:\n");  print ((r4.ig <- igibbs_mvn (4)) $ q)

cat("\na:\n");   print ((r4.a <-  over_mvn (4, -0.95)) $ q)
cat("\nia:\n");  print ((r4.ia <- iover_mvn (4, -0.95)) $ q)

# Tests with 10000 iterations, for comparing with the correct distribution.

cat("LONG TESTS\n\n")

cat("j:\n");   mvn.look(r <- r10000.j <-   jmet_mvn  (10000, c(0.2,0.4)))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("lj:\n");  mvn.look(r <- r10000.lj <-  ljmet_mvn (10000, list(a=0.2,b=0.4)))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("ij:\n");  mvn.look(r <- r10000.ij <-  ijmet_mvn (10000, c(0.2,0.4)))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("ilj:\n"); mvn.look(r <- r10000.ilj <- iljmet_mvn(10000, list(a=0.2,b=0.4)))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")

cat("s:\n");   mvn.look(r <- r10000.s <-   smet_mvn (10000, 0.5))
  cat("\n")
cat("is:\n");  mvn.look(r <- r10000.is <-  ismet_mvn (10000, 0.5))
  cat("\n")
cat("o:\n");   mvn.look(r <- r10000.o <-   omet_mvn (10000, 0.5))
  cat("\n")
cat("lo:\n");  mvn.look(r <- r10000.lo <-  lomet_mvn (10000, 0.5))
  cat("\n")
cat("ilo:\n"); mvn.look(r <- r10000.ilo <- ilomet_mvn (10000, 0.5))
  cat("\n")

cat("c:\n");   mvn.look(r <- r10000.c <-   cmet_mvn (10000, 1.5))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("lc:\n");  mvn.look(r <- r10000.lc <-  cmet_mvn (10000, 1.5))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")

cat("b:\n");   mvn.look(r <- r10000.b <-   bhmc_mvn (10000, 0.2))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("b2:\n");  mvn.look(r <- r10000.b2<-   bhmc_mvn (10000, 0.28))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("b3:\n");  mvn.look(r <- r10000.b3<-   bhmc_mvn (10000, 0.28,rand.step=0.1))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("lb:\n");  mvn.look(r <- r10000.lb <-  lbhmc_mvn (10000, 0.2))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("ilb:\n"); mvn.look(r <- r10000.ilb <- ilbhmc_mvn (10000, 0.2))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")
cat("bt:\n");  mvn.look(r <- r10000.bt <-  bhmc_traj_mvn (10000, 0.2))
  cat("mean apr:",mean(r$ave[,"apr.1"]),"\n\n")

cat("pb:\n");  mvn.look(r <- r10000.pb <-  pbhmc_mvn (10000, 0.2))
  cat("mean apr:",mean(r$ave[,"apr.2"]),"\n\n")

cat("l:\n");   mvn.look(r <- r10000.l <-   lhmc_mvn (10000, 0.1))
  cat("\n")
cat("ill:\n"); mvn.look(r <- r10000.ill <- illhmc_mvn (10000, 0.1))
  cat("\n")
cat("sl:\n");  mvn.look(r <- r10000.sl <-  slhmc_mvn (10000, 0.1))
  cat("\n")

cat("g:\n");   mvn.look(r <- r10000.g <-   gibbs_mvn (10000))
  cat("\n")
cat("ig:\n");  mvn.look(r <- r10000.ig <-  igibbs_mvn (10000))
  cat("\n")

cat("a:\n");   mvn.look(r <- r10000.a <-   over_mvn (10000, -0.95))
  cat("\n")
cat("ia:\n");  mvn.look(r <- r10000.ia <-  iover_mvn (10000, -0.95))
  cat("\n")

