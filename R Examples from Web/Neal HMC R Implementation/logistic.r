# LOG DENSITY FUNCTION FOR POSTERIOR OF A LOGISTIC REGRESSION MODEL.
#
# Radford M. Neal, 2012.


# LOG DENSITY FOR THE POSTERIOR OF A LOGISTIC MODEL.  
#
# The covariates (a matrix X with n rows) and responses (a binary vector y 
# of length n) are provided as elements of the state, which will usually
# be specified with the fix argument of mcmc.  
#
# The regression coefficients are in the element beta, and the constant term
# is in alpha.  Given the value of a hyperparameter log.sigma, the beta 
# parameters are independent in the prior, with each of their priors being a 
# t distribution centred at zero, with scale exp(log.sigma), and degrees of 
# freedom given by the df argument.  To specify df, use an operation such as
#
#     formals(lpr_logistic) $ df <- 4
#
# If this is not done, the default for df is Inf, which gives the equivalent 
# of a normal prior.  The prior for alpha is normal with mean zero and standard 
# deviation alpha.sd, which can be set using "formals" as above; it's value
# if not set is 10.
#
# The log.sigma hyperparameter (the log of the prior standard deviation of 
# the beta parameters when df is Inf) is given by the log.sigma element of 
# the state.  Its prior is normal with mean log.sigma.mean and standard 
# deviation log.sigma.sd, which are parameters of lpr_logistic, which may be 
# set using "formals"; the values are mean 0 and standard deviations 2 if not 
# otherwise set.
#
# This function can compute the gradient with respect to alpha and beta,
# but not with respect to log.sigma.

lpr_logistic <- function (value, grad=FALSE,  # remainder are set with "formals"
                          df=Inf, alpha.sd=10, log.sigma.mean=0, log.sigma.sd=2)
{
  log.prior <- dnorm (value$log.sigma, log.sigma.mean, log.sigma.sd, log=TRUE) +
               dnorm (value$alpha, 0, alpha.sd, log=TRUE) +
               sum (dt (value$beta/exp(value$log.sigma), df, log=TRUE)) - 
                 length(value$beta) * value$log.sigma

  t <- 2*value$y - 1                       # responses in -1 / +1 form
  s <- as.vector (value$X %*% value$beta)  # sums of beta times covariates

  pt <- 1 / (1+exp(-(value$alpha+s)*t))
  log.likelihood <- sum(log(pt))

  r <- log.prior + log.likelihood

  if (grad)
  { d <- t * (1-pt)
    g <- list (alpha = sum(d) - value$alpha/alpha.sd^2, beta = colSums(d*X))
    if (is.finite(df))
    { g$beta <- g$beta - 
          (df+1) * value$beta / (exp(value$log.sigma)^2 * df + value$beta^2)
    }
    else
    { g$beta <- g$beta - value$beta/exp(value$log.sigma)^2
    }
    attr(r,"grad") <- g
  }

  r
}


# INCREMENTAL VERSION OF LOG DENSITY FOR LOGISTIC MODEL POSTERIOR.  The
# quantities saved to help with re-computation are the log of the prior 
# density (log.prior), the sums of regression coefficients times covariates
# (s), and the log of the likelihood (log.likelihood).
#
# This function cannot compute the gradient.

ilpr_logistic <- function (value, ch.value, ch.elem, ch.pos, lpr.value, 
                          df=Inf, alpha.sd=10, log.sigma.mean=0, log.sigma.sd=2)
{
  if (!missing(ch.pos) && length(ch.pos)!=length(ch.value)) 
  { stop("ch.value and ch.pos are incompatible")
  }

  # Initially, nothing has been computed.

  log.likelihood <- NULL
  log.prior <- NULL
  s <- NULL

  # Start with parameters from value argument.

  log.sigma <- value$log.sigma
  alpha <- value$alpha
  beta <- value$beta

  # If only part of the state has changed, recompute some things incrementally,
  # and change log.sigma, alpha, or beta as required.
    
  if (!missing(ch.value))
  { 
    lp <- attr (lpr.value, "log.prior")

    if (ch.elem=="log.sigma")
    { log.likelihood <- attr (lpr.value, "log.likelihood")
      s <- attr (lpr.value, "s")
      log.sigma <- ch.value
    }

    else if (ch.elem=="alpha")
    { log.prior <- lp - dnorm (alpha, 0, alpha.sd, log=TRUE) +
                        dnorm (ch.value, 0, alpha.sd, log=TRUE)
      s <- attr (lpr.value, "s")
      alpha <- ch.value
    }

    else if (ch.elem=="beta")
    { if (!missing(ch.pos))
      { log.prior <- lp - sum (dt (beta[ch.pos]/exp(log.sigma), df, log=TRUE)) +
                          sum (dt (ch.value/exp(log.sigma), df, log=TRUE))
        s <- attr (lpr.value, "s") + 
             value$X[,ch.pos,drop=FALSE] %*% (ch.value - beta[ch.pos])
      }
      beta <- ch.value
    }
  } 

  # Compute whichever of log.prior, s, and log.likelihood haven't been
  # computed incrementally above.

  if (is.null(log.prior))
  { log.prior <- dnorm (log.sigma, log.sigma.mean, log.sigma.sd, log=TRUE) + 
                 dnorm (alpha, 0, alpha.sd, log=TRUE) +
                 sum (dt (beta/exp(log.sigma), df, log=TRUE)) - 
                   length(beta) * log.sigma
  }

  if (is.null(s))
  { s <- as.vector (value$X %*% beta)  # sums of beta times covariates
  }

  if (is.null(log.likelihood))
  { t <- 2*value$y - 1                 # responses in -1 / +1 form
    log.likelihood <- sum (-log(1+exp(-(alpha+s)*t)))
  }

  # Return the log density along with cached values for possible later use.

  structure (log.prior + log.likelihood, 
             log.likelihood = log.likelihood, 
             log.prior = log.prior, 
             s = s)
}


# SPECIAL UPDATE FUNCTION TO SAMPLE FROM THE PRIOR IN A LOGISTIC MODEL.
#
# The level argument indicates what level of the model to start sampling at:
#
#   Level 0 samples everything except covariates (X) from the prior.  
#   Level 1 samples everything except covariates and log.sigma.  
#   Level 2 samples only the responses (y).  
#
# The lpr function passed is not called, but its formals are examined to 
# get alpha.sd, log.sigma.mean, and log.sigma.sd, and its attributes are 
# consulted for bounds.
#
# Bounds (if any) are imposed by repeatedly sampling until they are satisfied.
#
# NOTE:  This update will not leave the correct distribution invariant if
# the level specified samples things that are not part of the mcmc state.

sample_logistic_prior <- structure (special=TRUE, handles.bounds=TRUE,
                         function (lpr, initial, level)
{
  if (level<0 || level>2)
  { stop("invalid level argument")
  }

  formals <- formals(lpr)

  repeat 
  { 
    if (level < 1)
    { initial$log.sigma <- 
        rnorm (1, formals$log.sigma.mean, formals$log.sigma.sd)
    }

    if (level < 2)
    { initial$alpha <- rnorm (1, 0, formals$alpha.sd)
      initial$beta <- rnorm (ncol(initial$X), 0, exp(initial$log.sigma))
    }

    s <- as.vector (initial$X %*% initial$beta)
    p1 <- 1 / (1+exp(-(initial$alpha+s)))
    initial$y <- as.numeric (runif(length(p1)) < p1)

    lower <- attr(lpr,"lower")
    upper <- attr(lpr,"upper")

    if ( is.null(lower) && is.null(upper)
      || (is.null(lower$log.sigma) || initial$log.sigma >= lower$log.sigma)
      && (is.null(upper$log.sigma) || initial$log.sigma <= upper$log.sigma)
      && (is.null(lower$alpha) || initial$alpha >= lower$alpha)
      && (is.null(upper$alpha) || initial$alpha <= upper$alpha)
      && (is.null(lower$beta) || all(initial$beta >= lower$beta))
      && (is.null(upper$beta) || all(initial$beta <= upper$beta))
      && (is.null(lower$y) || all(initial$y >= lower$y))
      && (is.null(upper$y) || all(initial$y <= upper$y)) )
    {
      break
    }
  }

  list (final = initial)
})
