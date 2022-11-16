# LOG DENSITY FUNCTION FOR A MIXTURE OF TWO MULTIVARIATE NORMAL DISTRIBUTIONS.
#
# Radford M. Neal, 2012.


# LOG DENSITY FUNCTION FOR MIXTURE OF TWO MULTIVARIATE NORMALS.  The p1, mu1, 
# mu2, inv.cov1, and inv.cov2 arguments should be changed before use to have 
# the desired values as their default.  This might best be done on a copy.
#
# The gradient will be computed if asked for.

lpr_mix <- function (value, grad=FALSE, p1, mu1, mu2, inv.cov1, inv.cov2)
{
  D <- length(value)

  g1 <- - inv.cov1 %*% (value-mu1)
  g2 <- - inv.cov2 %*% (value-mu2)

  lpr1 <- - (D/2) * log(2*pi) + sum(log(diag(chol(inv.cov1)))) + 
            as.vector (t(value-mu1)%*%g1/2)
  lpr2 <- - (D/2) * log(2*pi) + sum(log(diag(chol(inv.cov2)))) + 
            as.vector (t(value-mu2)%*%g2/2)

  r <- log_average_exp (c(lpr1,lpr2), c(p1,1-p1))

  q1 <- p1 * exp(lpr1-r)

  if (grad)
  { attr(r,"grad") <- q1 * as.vector(g1) + (1-q1) * as.vector(g2)
  }

  r
}
