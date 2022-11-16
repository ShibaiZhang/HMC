# HAMILTONIAN MONTE CARLO WITH ACCEPT/REJECT WINDOWS.
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
#    window       Size of the accept and reject windows, from 1 to nsteps+1.
#                 May have an attribute "weights" that is a vector of this
#                 length of non-negative weights for window positions, 
#                 not all zero.  The weights need not be normalized to sum 
#                 to one. The weights given are for the reject window; weights
#                 for the accept window are the reverse.  Default is 1, giving
#                 standard HMC.
#    return.traj  TRUE if values of q and p along the trajectory should be 
#                 returned (default is FALSE).
#
# The value returned is a list containing the following elements:
#
#    final        The new position part of the state
#    final.p      The new momentum part of the state
#    lpr          The value of lpr(final,grad=TRUE)
#    step         Stepsize(s) used (after jittering)
#    acc          1 if the accept window was used, 0 if not
#    apr          The probability of using the accept window
#    delta        Change in H that the acceptance decision was based on
#    ipos         Position of initial state within the reject window, from 0
#    move         Position of final state minus position of initial state
#
# plus the following elements if return.traj was TRUE:
#
#    traj.q     Matrix of nsteps+1 rows with position values along trajectory
#    traj.p     Matrix of nsteps+1 rows with momentum values along trajectory
#    traj.H     Vector of length nsteps+1 with values of H along trajectory

windowed_hmc <- function (lpr, initial, initial.p = rnorm(length(initial)), 
                          lpr.initial = NULL, nsteps = 1, step, rand.step = 0,
                          window = 1, return.traj = FALSE)
{
  # Check and process the arguments.

  step <- process_step_arguments (length(initial), step, rand.step)
  process_nsteps_argument (nsteps)

  # Check window argument and set up reject and accept window weights.

  if (!is.numeric(window) || length(window)!=1 || window!=floor(window) 
   || window<1 || window>nsteps+1)
  { stop("Bad window argument")
  }

  rwin <- attr(window,"weights")

  if (is.null(rwin))
  { rwin <- rep(1,window)
  }
  else 
  { if (!is.numeric(rwin) || length(rwin)!=window 
     || any(rwin<0) || all(rwin==0))
    { stop("Bad window weights")
    }
  }

  awin <- rev(rwin)

  one.window <- window==nsteps+1 && all(awin==rwin)

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

  # Compute the Hamiltonian at the start of the trajectory.

  kinetic.initial <- sum(initial.p^2) / 2
  initial.H <- kinetic.initial - lpr.initial

  # Randomly decide on the position of the initial state within the reject 
  # window (from 0).  Don't call 'sample' when window size is one, for speed,
  # and so results will match basic_hmc (useful for testing).

  ipos <- if (window==1) 0 else sample(window,1,prob=rwin) - 1

  # Set up initial data on the accept and reject windows, accounting only
  # for the initial state.

  pos <- ipos

  if (!one.window)
  { n.rej <- 1
    rej.q <- initial
    rej.p <- initial.p
    rej.F <- initial.H - log(rwin[pos+1])
    rej.lpr <- lpr.initial
    rej.pos <- pos
  }

  if (pos <= nsteps-window)
  { n.acc <- 0
  }
  else
  { n.acc <- 1
    acc.q <- initial
    acc.p <- initial.p
    acc.F <- initial.H - log(awin[pos-nsteps+window])
    acc.lpr <- lpr.initial
    acc.pos <- pos
  }

  # Compute the trajectory by the leapfrog method, part backward from the
  # initial state, and part forward from the initial state.  Update states
  # saved as the final state from the reject window, and as the final state
  # from the accept window, and record minus the log of total probability 
  # (minus the log weight) for each window to decide on which to use.

  q <- initial
  p <- initial.p
  gr <- attr(lpr.initial,"grad")

  if (pos==0) 
  { dir <- 1
  }
  else
  { dir <- -1
    step <- -step
  }

  if (return.traj)
  { traj.q[1+pos,] <- q
    traj.p[1+pos,] <- p
    traj.H[1+pos] <- kinetic.initial - lpr.initial
  }

  while (n.acc < window)
  {
    # Do a leapfrog update (forward or backward, depending on step).

    p <- p + (step/2) * gr
    q <- q + step * p     
    lr <- lpr(q,grad=TRUE)
    gr <- attr(lr,"grad")
    p <- p + (step/2) * gr

    # Update position of state.

    pos <- pos + dir

    # Compute Hamiltonian at this state.

    H <- sum(p^2)/2 - lr

    # Record trajectory if asked to.

    if (return.traj)
    { traj.q[1+pos,] <- q
      traj.p[1+pos,] <- p
      traj.H[1+pos] <- H
    }

    # Check for the state being in the reject window.

    if (pos < window && !one.window)
    {
      n.rej <- n.rej+1
      Hw <- H - log(rwin[pos+1])
      rej.F <- min(rej.F,Hw) - log(1+exp(-abs(rej.F-Hw)))

      if (runif(1) < exp(rej.F-Hw))
      { rej.q <- q
        rej.p <- p
        rej.lpr <- lr
        rej.pos <- pos
      }
    }

    # Check for the state being in the accept window.

    if (pos > nsteps-window)
    { 
      n.acc <- n.acc+1
      Hw <- H - log(awin[pos-nsteps+window])

      if (n.acc==1)
      { acc.F <- Hw
        acc.q <- q
        acc.p <- p
        acc.lpr <- lr
        acc.pos <- pos
      }
      else 
      { acc.F <- min(acc.F,Hw) - log(1+exp(-abs(acc.F-Hw)))
        if (runif(1) < exp(acc.F-Hw))
        { acc.q <- q
          acc.p <- p
          acc.lpr <- lr
          acc.pos <- pos
        }
      }
    }

    # Switch to forward simulation if now finished with backward simulation.

    if (pos==0)
    { pos <- ipos
      q <- initial
      p <- initial.p
      gr <- attr(lpr.initial,"grad")
      dir <- 1
      step <- -step
    }
  }

  # Decide whether the final state comes from the accept window or reject
  # window, and set the final state accordingly.

  if (one.window)
  { delta <- 0
    apr <- 1
  }
  else
  { delta <-  acc.F - rej.F
    apr <- min(1,exp(-delta))
  }

  if (one.window || runif(1) < apr) # ACCEPT
  { final.q <- acc.q
    final.p <- - acc.p # negate so reversible
    lpr.final <- acc.lpr
    move <- acc.pos - ipos
    acc <- 1
  }
  else # REJECT
  { final.q <- rej.q
    final.p <- rej.p
    lpr.final <- rej.lpr
    move <- rej.pos - ipos
    acc <- 0
  }

  # Return new state, its log probability with gradient, plus additional
  # information, including the the trajectory, if requested.

  r <- list (final=final.q, final.p=final.p, lpr=lpr.final, step=step, 
             apr=apr, acc=acc, delta=delta, ipos=ipos, move=move)

  if (return.traj) 
  { r$traj.q <- traj.q 
    r$traj.p <- traj.p
    r$traj.H <- traj.H 
  }

  r
}
