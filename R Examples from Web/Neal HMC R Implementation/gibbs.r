# GIBBS SAMPLING AND RELATED UPDATES.
#
# Radford M. Neal, 2012.


# These functions perform univariate Gibbs sampling for some standard
# distributions, by inferring the parameters of the distribution from
# a few values of the log density.  They are meant to be used as an
# update within singlevar.


# GIBBS SAMPLING FOR A UNIVARIATE NORMAL.  Finds the mean and
# standard deviation by evaluating the log density at initial,
# initial+1, and initial-1.  Standard Adler overrelaxation may be 
# done by specifying a non-zero value for a.  
#
# The density must actually be normal for this to work correctly!

gibbs_normal <- function (lpr, initial, lpr.initial=NULL, a=0)
{
  d0 <- if (is.null(lpr.initial)) lpr(initial) else lpr.initial
  dp <- lpr(initial+1)
  dm <- lpr(initial-1)

  v <- 1 / (2*d0 - dp - dm)
  m <- initial + v * (dp-dm) / 2

  final <- if (a==0) rnorm (1, m, sqrt(v))
           else      rnorm (1, m + a * (initial-m), sqrt((1-a^2)*v))

  list (final = final)
}
