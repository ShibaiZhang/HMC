# HAMILTONIAN MONTE CARLO WITH TEMPERING
#
# Radford M. Neal, 2012.
#
# Arguments:
#
#    lpr          Function returning the log probability of the position part 
#                 of the state, plus an arbitrary constant, with gradient
#                 as an attribute if grad=TRUE is passed.
#    initial      The initial position part of the state (a vector).
#    initial.p    The initial momentum part of the state (a vector), default
#                 is to use momentum variables generated as standard normals.
#    lpr.initial  The value of lpr(initial), maybe with grad.  May be omitted.
#    nsteps       Number of steps in trajectory used to propose a new state.
#                 (Default is 1, giving the "Langevin" method.)
#    step         Stepsize or stepsizes.  May be scalar or a vector of length 
#                 equal to the dimensionality of the state.
#    rand.step    Amount of random jitter for step. The step is multiplied
#                 by exp (runif (length(rand.step), -rand.step, rand.step)).
#                 A scalar or vector of length equal to the dimensionality. 
#                 Default is 0.
#    temper       Tempering factor.  Must be positive, usually slightly greater
#                 than 1.  Default is 1.
#    return.traj  TRUE if values of q and p along the trajectory should be 
#                 returned (default is FALSE).
#
# The value returned is a list containing the following elements:
#
#    final        The new position part of the state
#    final.p      The new momentum part of the state
#    lpr          The value of lpr(final,grad=TRUE)
#    step         Stepsize(s) used (after jittering)
#    acc          1 if the proposal was accepted, 0 if rejected
#    apr          Acceptance probability of the proposal
#    delta        Change in H the acceptance decision was based on
#
# plus the following elements if return.traj was TRUE:
#
#    traj.q     Matrix of nsteps+1 rows with position values along trajectory
#    traj.p     Matrix of nsteps+1 rows with momentum values along trajectory
#    traj.H     Vector of length nsteps+1 with values of H along trajectory

tempered_hmc <- function (lpr, initial, initial.p = rnorm(length(initial)), 
                          lpr.initial = NULL, nsteps = 1, step, rand.step = 0, 
                          temper = 1, return.traj = FALSE)
{
  # Check and process the arguments.

  step <- process_step_arguments (length(initial), step, rand.step)
  process_nsteps_argument (nsteps)
  if (!is.numeric(temper) || length(temper)!=1 || temper<=0)
  { stop("Bad temper argument")
  }

  # Allocate space for the trajectory, if its return is requested.

  if (return.traj)
  { traj.q <- matrix(NA,nsteps+1,length(initial))
    traj.p <- matrix(NA,nsteps+1,length(initial))
    traj.H <- rep(NA,nsteps+1)
  }

  # Evaluate the log probability and gradient at the initial position, if not 
  # already known.

  if (is.null(lpr.initial) || is.null(attr(lpr.initial,"grad")))
  { lpr.initial <- lpr(initial,grad=TRUE)
  }

  # Compute the kinetic energy at the start of the trajectory.

  kinetic.initial <- sum(initial.p^2) / 2

  # Compute the trajectory by the leapfrog method.

  q <- initial
  p <- initial.p
  gr <- attr(lpr.initial,"grad")

  if (return.traj)
  { traj.q[1,] <- q
    traj.p[1,] <- p
    traj.H[1] <- kinetic.initial - lpr.initial
  }

  # Do leapfrog steps with tempering before/after steps.

  sqrt.temper <- sqrt(temper)

  for (i in 1:nsteps)
  { 
    # Before each leapfrog step, multiply in first half, divide in second half.

    if (2*(i-1)<nsteps)
    { p <- p * sqrt.temper
    }
    else
    { p <- p / sqrt.temper
    }

    # Leapfrog step.

    p <- p + (step/2) * gr
    q <- q + step * p     
    lr <- lpr(q,grad=TRUE)
    gr <- attr(lr,"grad")
    p <- p + (step/2) * gr

    # After each leapfrog step, divide in second half, multiply in first half.

    if (2*i>nsteps)
    { p <- p / sqrt.temper
    }
    else
    { p <- p * sqrt.temper
    }

    # Record q, p, and H for the trajectory, if requested.

    if (return.traj)
    { traj.q[i+1,] <- q
      traj.p[i+1,] <- p
      traj.H[i+1] <- sum(p^2)/2 - lr
    }
  }

  # Negate momentum at end of trajectory to make the proposal symmetric.

  p <- -p

  # Look at log probability and kinetic energy at the end of the trajectory.

  lpr.prop <- lr
  kinetic.prop <- sum(p^2) / 2

  # Accept or reject the state at the end of the trajectory.

  H.initial <- -lpr.initial + kinetic.initial
  H.prop <- -lpr.prop + kinetic.prop
  delta <- H.prop - H.initial
  apr <- min(1,exp(-delta))

  if (runif(1) < apr) # accept
  { 
    final.q <- q
    final.p <- p
    lpr.final <- lr
    acc <- 1
  }
  else # reject
  { 
    final.q <- initial
    final.p <- initial.p
    lpr.final <- lpr.initial
    acc <- 0
  }

  # Return new state, its log probability and gradient, plus additional
  # information, including the trajectory, if requested.

  r <- list (final=final.q, final.p=final.p, lpr=lpr.final, step=step, 
             apr=apr, acc=acc, delta=delta)

  if (return.traj) 
  { r$traj.q <- traj.q 
    r$traj.p <- traj.p
    r$traj.H <- traj.H 
  }

  r
}
