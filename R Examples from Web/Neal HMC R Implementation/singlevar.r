# DO A SEQUENCE OF SINGLE-VARIABLE UPDATES.
#
# Radford M. Neal, 2012.

singlevar <- function (lpr, initial, update, lpr.initial=NULL,
                       seq=1:length(initial), random=FALSE, ...)
{
  # Check that additional arguments are named, and numeric arguments are 
  # either scalars or vectors that give a value for each variable in initial.
  
  args <- list(...)

  if (any(names(args)==""))
  { stop("Arguments to pass on to the update must be named")
  }

  for (a in args)
  { if (is.numeric(a) && length(a)!=1 && length(a)!=length(initial))
    { stop("Invalid argument to pass on to the update")
    }
  }

  # Look at the update, ensure it's a list, and see if it wants lpr.initial.

  if (!is.list(update))
  { update <- list(update)
  }
  if (!is.function(update[[1]]))
  { stop("MCMC update description must start with a function")
  }
  has.lpr.initial <- "lpr.initial" %in% names(formals(update[[1]]))

  # Define the wrapper for lpr that handles just one variable, either by
  # incremental computation or by modifying that variable explicitly.

  if ("ch.pos" %in% names(formals(lpr)))
  { 
    if ("grad" %in% names(formals(lpr)))
    { lpr.single <- function (value, grad=FALSE)
      { if (is.null(lpr.curr))
        { curr[i] <<- value; 
          r <- lpr (curr, grad=grad); 
        }
        else
        { r <- lpr (curr, grad=grad, 
                    ch.value=value, ch.pos=i, lpr.value=lpr.curr)
        }
        if (!is.null(attr(r,"grad")))
        { attr(r,"grad") <- attr(r,"grad")[i]
        }
        r
      }
    }
    else # no grad argument for lpr
    { lpr.single <- function (value)
      { if (is.null(lpr.curr))
        { curr[i] <<- value; 
          lpr (curr); 
        }
        else
        { lpr (curr, ch.value=value, ch.pos=i, lpr.value=lpr.curr)
        }
      }
    }
  }
  else # no ch.pos argument
  { 
    if ("grad" %in% names(formals(lpr)))
    { lpr.single <- function (value, grad=FALSE) 
      { curr[i] <<- value; 
        r <- lpr (curr, grad=grad); 
        if (!is.null(attr(r,"grad")))
        { attr(r,"grad") <- attr(r,"grad")[i]
        }
        r
      }
    }
    else # no grad argument for lpr
    { lpr.single <- function (value) 
      { curr[i] <<- value; 
        lpr(curr); 
      }
    }
  }

  # Randomize sequence order if asked to.

  if (random) 
  { ord <- sample(length(seq))
    seq[ord] <- seq
  }
  else
  { ord <- 1:length(seq)
  }

  # Apply the update to the sequence of variables.

  curr <- initial
  lpr.curr <- lpr.initial
  res <- NULL

  lower <- attr(lpr,"lower")
  upper <- attr(lpr,"upper")

  for (i in seq)
  { 
    attr(lpr.single,"lower") <- if (length(lower)==1) lower else lower[i]
    attr(lpr.single,"upper") <- if (length(upper)==1) upper else upper[i]

    a <- lapply (args, function (x) 
                         if (!is.numeric(x) || length(x)==1) x else x[i])

    if (has.lpr.initial)
    { r <- eval (as.call 
           (c(list (update[[1]], lpr.single, curr[i], lpr.initial = lpr.curr),
              update[-1], a)))
    }
    else
    { r <- eval (as.call 
            (c (list (update[[1]], lpr.single, curr[i]), update[-1], a)))
    }

    if (is.null(r$lpr) && has.lpr.initial)
    { r$lpr <- lpr.single(r$final)
    }

    curr[i] <- r$final
    lpr.curr <- r$lpr
    attr(lpr.curr,"grad") <- NULL

    r$final <- NULL
    r$lpr <- NULL

    res <- if (is.null(res)) r else mapply (base::c, res, r, SIMPLIFY=FALSE)
  }

  c (list (final=curr, lpr=lpr.curr), lapply (res, function (x) x[ord]))
}
