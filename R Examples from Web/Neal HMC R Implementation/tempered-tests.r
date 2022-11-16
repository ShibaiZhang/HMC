# TESTS USING TEMPERED HMC ON A MIXTURE DISTRIBUTION.
#
# Radford M. Neal, 2012.

source("mcmc.r")
source("utilities.r")
source("metropolis.r")
source("basic_hmc.r")
source("tempered_hmc.r")
source("pupdate.r")
source("mix.r")

# --------------------------------------------------------------------

# FUNCTION TO LOOK AT RESULTS OF THE TESTS BELOW.

mix.look <- function (r, thin=1)
{ if (is.list(r$q))
  { q1 <- r$q[[1]]
    q2 <- r$q[[2]]
  }
  else
  { q1 <- r$q[,1]
    q2 <- r$q[,2]
  }
  q1 <- q1[seq(1,length(q1),by=thin)]
  q2 <- q2[seq(1,length(q2),by=thin)]
  cat("m1:",round(mean(q1),4)," m2:",round(mean(q2),3),
      " v1:",round(var(q1),4)," v2:",round(var(q2),3),
      " act1:",round (2*sum(acf(q1,plot=FALSE,lag.max=200)$acf)-1, 1),
      " act2:",round (2*sum(acf(q2,plot=FALSE,lag.max=200)$acf)-1, 1),
      "\n")
  plot(q1,q2,pch=".")
}


# --------------------------------------------------------------------

# FUNCTIONS FOR SAMPLING THE MIXTURE DISTRIBUTION IN VARIOUS WAYS.  


# SAMPLE MIXTURE USING JOINT METROPOLIS UPDATES.

jmet_mix <- function (iterations, step, ...)
{
  mcmc (lpr_mix, c(a=0, b=0), iterations, 
        list (metropolis, step=step, rep=20),
        ...
  )
}


# SAMPLE MIXTURE USING STANDARD HMC.

hmc_mix <- function (iterations, step, nsteps=20, ...)
{
  mcmc (lpr_mix, c(a=0, b=0), iterations, 
        list (basic_hmc, step=step, nsteps=nsteps),
        ...
  )
}


# SAMPLE MIXTURE USING TEMPERED HMC.

thmc_mix <- function (iterations, step, nsteps=20, temper, ...)
{
  mcmc (lpr_mix, c(a=0, b=0), iterations, 
        list (tempered_hmc, step=step, nsteps=nsteps, temper=temper),
        ...
  )
}


# SAMPLE MIXTURE USING TEMPERED HMC + METROPOLIS.

mthmc_mix <- function (iterations, step, nsteps=20, temper, ...)
{
  mcmc (lpr_mix, c(a=0, b=0), iterations, 
        list (tempered_hmc, step=step, nsteps=nsteps, temper=temper),
        list (metropolis, step=0.3, rep=3),
        ...
  )
}


# --------------------------------------------------------------------

# SCRIPT FOR RUNNING THE TESTS.

formals(lpr_mix)$p1 <- 0.5
formals(lpr_mix)$mu1 <- c(0,0)
formals(lpr_mix)$mu2 <- c(10,10)
formals(lpr_mix)$inv.cov1 <- diag(2)
formals(lpr_mix)$inv.cov2 <- 0.5*diag(2)

it <- 10000

cat("\njmet:\n");
mix.look (rjmet <- jmet_mix (it, 3)); cat("\n")

cat("\nhmc:\n");
mix.look (rhmc <- hmc_mix (it, 1)); cat("\n")

cat("\nthmc1:\n");
mix.look (rthmc1 <- thmc_mix (it, 1, temper=1.1)); cat("\n")

cat("\nthmc2:\n");
mix.look (rthmc2 <- thmc_mix (it, 1, temper=1.2)); cat("\n")

cat("\nthmc3:\n");
mix.look (rthmc3 <- thmc_mix (it, 1, temper=1.4)); cat("\n")

cat("\nthmc4:\n");
mix.look (rthmc4 <- thmc_mix (it, 1, temper=1.6)); cat("\n")

cat("\nthmc5:\n");
mix.look (rthmc5 <- thmc_mix (it, 1, temper=1.9)); cat("\n")

cat("\nthmc6:\n");
mix.look (rthmc6 <- thmc_mix (it, 1, temper=2.1)); cat("\n")

cat("\nthmc7:\n");
mix.look (rthmc7 <- thmc_mix (it, 1, temper=2.5)); cat("\n")

cat("\nmthmc5:\n");
mix.look (rmthmc5 <- mthmc_mix (it, 1, temper=1.9)); cat("\n")

cat("\nmtlmc:\n");
mix.look (rmtlmc <- mthmc_mix (it, 1, nsteps=1, temper=36)); cat("\n")


# --------------------------------------------------------------------

# PRODUCE PLOT AS IN FIGURE 9 OF "MCMC USING HAMILTONIAN DYNAMICS".

par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2,0.7,0))

ss <- 0.3
ns <- 200
temper <- 1.04

t <- tempered_hmc (lpr_mix, 
                   step=ss, nsteps=ns, temper=temper,
                   initial=c(-0.4,-0.9), initial.p=c(0.7,-0.9),
                   return.traj=TRUE)

plot (0:ns, t$traj.H, ylim=c(0,62),
      pch=20, xlab="Leapfrog step number", ylab="Value of Hamiltonian")

plot (t$traj.q[,1], t$traj.q[,2], type="l", xlim=c(-10,25), ylim=c(-10,25),
      xlab="Position coordinate 1", ylab="Position coordinate 2")

t1 <- t

t <- tempered_hmc (lpr_mix, 
                   step=ss, nsteps=ns, temper=temper,
                   initial=c(0.1,1.0), initial.p=c(0.5,0.8),
                   return.traj=TRUE)

print(t$traj.H[ns+1]-t$traj.H[1])

plot (0:ns, t$traj.H, ylim=c(0,62),
      pch=20, xlab="Leapfrog step number", ylab="Value of Hamiltonian")

plot (t$traj.q[,1], t$traj.q[,2], type="l", xlim=c(-10,25), ylim=c(-10,25),
      xlab="Position coordinate 1", ylab="Position coordinate 2")

t2 <- t
