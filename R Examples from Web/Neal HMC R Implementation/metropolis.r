# RANDOM WALK METROPOLIS UPDATE WITH NORMAL PROPOSALS.
#
# Radford M. Neal, 2012.
#
# Arguments:
#
#    lpr         Function returning the log probability of a state, plus an
#                arbitrary constant
#    initial     The initial state (a vector)
#    lpr.initial The value of lpr(initial).  May be omitted.
#    step        Stepsize or stepsizes - a scalar or a vector of length equal
#                to the dimensionality of the state.
#    rand.step   Amount of random jitter for step (scalar or vector the
#                lenght of initial)
#    rep         Number of times to repeat the update procedure (default 1).
#
# The value returned is a list containing the following elements:
#
#    final       The new state
#    lpr         The value of lpr(final)
#    step        The stepsize or stepsizes used (after jittering)
#    acc         0 or 1, indicating if the last (or only) update was
#                accepted (1) or rejected (0)
#    apr         Average acceptance probability
#    delta       Change in -lpr the last accept/reject decision was based on

metropolis <- function (lpr, initial, lpr.initial=NULL,
                        step, rand.step=0,  rep=1)
{
  # Check and process the arguments.

  step <- process_step_arguments (length(initial), step, rand.step)
  rep <- process_rep_argument (rep)

  # Evaluate the log probability of the initial state, if not provided.

  if (is.null(lpr.initial))
  { lpr.initial <- lpr(initial)
  }

  # Do 'rep' Metropolis updates.

  apr <- 0
  curr <- initial
  lpr.curr <- lpr.initial

  for (i in 1:rep)
  {
    # Propose a candidate state, and evaluate its log probability.

    proposal <- curr + rnorm (length(initial), 0, step)
    lpr.proposal <- lpr(proposal)
    delta = lpr.curr - lpr.proposal
  
    # Accept or reject the proposed state as the new state.
  
    ap <- min(1,exp(-delta))
    apr <- apr + ap

    if (runif(1) < ap) 
    { curr <- proposal
      lpr.curr <- lpr.proposal
      acc <- 1
    }
    else
    { acc <- 0
    }
  }

  # Return the new state, its log probability, the average acceptance 
  # probability, and the change in -lpr for the last update.

  list (final=curr, lpr=lpr.curr, step=step, acc=acc, apr=apr/rep, delta=delta)
}
