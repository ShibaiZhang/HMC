# MOMENTUM UPDATES FOR HAMILTONIAN MONTE CARLO
#
# Radford M. Neal, 2012.


# MOMENTUM UPDATE FOR HAMILTONIAN MONTE CARLO, FOR GAUSSIAN MOMENTUM
#
# Arguments:
#
#    lpr          Function returning the log probability of the position part 
#                 of the state, plus an arbitrary constant, with gradient
#                 as an attribute if grad=TRUE is passed.  Ignored here.
#    initial      The initial position part of the state (a vector). 
#    initial.p    The initial momentum part of the state (a vector), default
#                 is to use momentum variables generated as standard normals.
#    lpr.initial  The value of lpr(initial), maybe with grad.  May be omitted.
#    persistance  Factor by which to multiply old momentum by. (Default 0.)
#
# The value returned is a list containing the following elements:
#
#    final        The new position part of the state, same as initial
#    final.p      The new momentum part of the state
#    lpr          The value of lpr at final, copied from lpr.initial

pupdate <- function (lpr, initial, initial.p = rnorm(length(initial)), 
                          lpr.initial = NULL, persistence=0)
{
  if (!is.numeric(persistence) || length(persistence) != 1 
        || persistence < -1 || persistence > 1)
  { stop("bad persistence argument")
  }

  if (persistence != 0)
  { final.p <- persistence * initial.p + 
                 rnorm (length(initial.p), 0, sqrt(1-persistence^2))
  }
  else
  { final.p <- rnorm (length(initial.p))
  }

  list (final=initial, final.p=final.p, lpr=lpr.initial)
}
