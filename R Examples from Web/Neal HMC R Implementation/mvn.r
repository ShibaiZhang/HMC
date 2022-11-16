# LOG DENSITY FUNCTION FOR A MULTIVARIATE NORMAL DISTRIBUTION.
#
# Radford M. Neal, 2012.


# LOG DENSITY FUNCTION FOR A MULTIVARIATE NORMAL.  The inv.cov argument
# should be changed before use to have the desired value as its default.
# This might best be done on a copy, as in
#
#     my_lpr <- lpr_mvn; formals(my_lpr) $ inv.cov <- solve(my_cov)
#
# The gradient will be computed if asked for.

lpr_mvn <- function (value, grad=FALSE, inv.cov)
{
  g <- - inv.cov %*% value
  r <- as.vector (t(value) %*% g / 2)
  if (grad)
  { attr(r,"grad") <- as.vector(g)
  }

  r
}


# INCREMENTAL VERSION OF MULTIVARIATE NORMAL LOG DENSITY FUNCTION.  May be 
# faster when the dimensionality is high, and the updates used can exploit 
# incremental computation.

ilpr_mvn <- function (value, grad=FALSE, ch.value, ch.pos, lpr.value, inv.cov)
{
  if (missing(ch.value))
  { g <- - inv.cov %*% value
  }
  else
  { g <- attr(lpr.value,"g") - 
           inv.cov[,ch.pos,drop=FALSE] %*% (ch.value-value[ch.pos])
    value[ch.pos] <- ch.value
  }

  r <- as.vector (t(value) %*% g / 2)

  g <- as.vector(g)
  attr(r,"g") <- g
  if (grad)
  { attr(r,"grad") <- as.vector(g)
  }

  r
}
