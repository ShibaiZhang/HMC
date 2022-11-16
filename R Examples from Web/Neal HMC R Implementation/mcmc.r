# DO AN MCMC SIMULATION.  
#
# Radford M. Neal, 2012.

mcmc <- function (lpr, initial, iterations, ..., updates=list(), rep=1,
                  fix=NULL, initial.p=NULL, last="lpr", rec="acc", ave="apr",
                  lpr.initial=NULL, initial.arg=initial,
                  seed=1, saved.rand=NULL)
{
  # Check validity of simple arguments.

  if (!is.function(lpr))
  { stop("Log probability density function argument is not a function")
  }

  if (length(names(formals(lpr)))<1 || names(formals(lpr))[1]!="value")
  { stop("Log probability density function must have value as first argument")
  }

  if (!is.numeric(iterations) || length(iterations)!=1 || iterations<1)
  { stop("Invalid iterations argument")
  }

  if (!is.numeric(rep) || length(rep)!=1 || rep<1)
  { stop("Invalid rep argument")
  }

  if (!is.null(last) && !is.character(last))
  { stop("Invalid last argument")
  }

  if (!is.null(rec) && !is.character(rec))
  { stop("Invalid rec argument")
  }

  if (!is.null(ave) && !is.character(ave))
  { stop("Invalid ave argument")
  }

  # Check validity of initial state, and of initial.p and fix arguments.

  validate.mcmc.state(initial)

  if (!is.null(initial.p))
  { validate.mcmc.state (initial.p, initial)
  }

  if (!is.null(fix))
  { if (!is.list(fix))
    { stop("fix argument must be a list")
    }
    if (!is.list(initial))
    { stop("initial state must be a list when part of state is fixed")
    }
    validate.mcmc.state(fix)
    if (any (names(initial) %in% names(fix)))
    { stop("Part of the initial state is also fixed")
    }
  }

  # Get bounds attached to lpr and check their validity.

  lower <- attr(lpr,"lower")
  upper <- attr(lpr,"upper")

  if (!is.null(lower)) 
  { validate.mcmc.state (lower, initial, 
                         allow.scalar=TRUE, allow.infinite=TRUE)
  }

  if (!is.null(upper)) 
  { validate.mcmc.state (upper, initial, 
                         allow.scalar=TRUE, allow.infinite=TRUE)
  }

  # Check that fix argument satisfies bounds, then remove that part from 
  # the bounds.

  if (!is.null(fix))
  { 
    if (value_out_of_bounds (fix, lower, upper))
    { stop("Value given by fix argument is out of bounds")
    }

    lower[names(fix)] <- NULL
    if (length(lower)==0) lower <- NULL

    upper[names(fix)] <- NULL
    if (length(upper)==0) upper <- NULL
  }

  # Find out the characteristics of the lpr function.

  lpr.has.grad.arg <- "grad" %in% names(formals(lpr))
  lpr.has.ch.pos.arg <- "ch.pos" %in% names(formals(lpr))
  lpr.has.ch.elem.arg <- "ch.pos" %in% names(formals(lpr))

  if (lpr.has.ch.pos.arg || lpr.has.ch.elem.arg)
  { if (!all (c("ch.value","lpr.value") %in% names(formals(lpr))))
    { stop("lpr function has ch.elem or ch.pos without ch.value or lpr.initial")
    }
    if (is.list(initial) && lpr.has.ch.pos.arg && !lpr.has.ch.elem.arg)
    { stop("lpr function for list can't have ch.pos argument without ch.elem")
    }
  }

  # Create a wrapper function for lpr called lpr.imp to impose any bounds.

  if (is.null(lower) && is.null(upper)
       || !is.null(attr(lpr,"imposed")) && attr(lpr,"imposed"))
  { 
    lpr.imp <- lpr
  }
  else
  { 
    if (is.list(initial))  # state is a list
    { 
      lpr.imp <- function  # formals replaced below, but something like ...
                   (value, grad, ch.value, ch.elem, ch.pos, lpr.value)
      { 
        if (lpr.has.ch.elem.arg && !missing(ch.elem))
        { if (lpr.has.ch.pos.arg && !missing(ch.pos))
          { out.of.bounds <- value_out_of_bounds (ch.value,
              lower[[ch.elem]][ch.pos], upper[[ch.elem]][ch.pos])
          }
          else
          { out.of.bounds <- value_out_of_bounds (ch.value,
              lower[[ch.elem]], upper[[ch.elem]])
          }
        }
        else
        { out.of.bounds <- value_out_of_bounds (value, lower, upper)
        }

        if (out.of.bounds)
        { r <- -Inf
          if (lpr.has.grad.arg && !missing(grad) && grad)
          { attr(r,"grad") <- numeric(length(value))
          }
          r
        }
        else  # no bounds violation, just call original lpr function
        { if (lpr.has.grad.arg) 
            if (lpr.has.ch.pos.arg)
              lpr (value, grad=grad, ch.value=ch.value, ch.elem=ch.elem, 
                          ch.pos=ch.pos, lpr.value=lpr.value)
            else if (lpr.has.ch.elem.arg)
              lpr (value, grad=grad, ch.value=ch.value, ch.elem=ch.elem, 
                          lpr.value=lpr.value)
            else
              lpr (value, grad=grad)
          else
            if (lpr.has.ch.pos.arg)
              lpr (value, ch.value=ch.value, ch.elem=ch.elem, ch.pos=ch.pos, 
                          lpr.value=lpr.value)
            else if (lpr.has.ch.elem.arg)
              lpr (value, ch.value=ch.value, ch.elem=ch.elem, 
                          lpr.value=lpr.value)
            else
              lpr (value)
        }
      }
    }

    else  # state is a vector
    { 
      lpr.imp <- function  # formals replaced below, but something like ...
                   (value, grad, ch.value, ch.pos, lpr.value)
      { 
        if (lpr.has.ch.pos.arg && !missing(ch.pos))
        { out.of.bounds <- 
            value_out_of_bounds (ch.value, lower[ch.pos], upper[ch.pos])
        }
        else
        { out.of.bounds <-
            value_out_of_bounds (value, lower, upper)
        }

        if (out.of.bounds)
        { r <- -Inf
          if (lpr.has.grad.arg && !missing(grad) && grad)
          { attr(r,"grad") <- numeric(length(value))
          }
          r
        }
        else  # no bounds violation, just call original lpr function
        { if (lpr.has.grad.arg) 
            if (lpr.has.ch.pos.arg)
              lpr (value, grad=grad, ch.value=ch.value, ch.pos=ch.pos,
                          lpr.value=lpr.value)
            else
              lpr (value, grad=grad)
          else
            if (lpr.has.ch.pos.arg)
              lpr (value, ch.value=ch.value, ch.pos=ch.pos, 
                          lpr.value=lpr.value)
            else
              lpr (value)
        }
      }
    }

    formals(lpr.imp) <- formals(lpr)
    attributes(lpr.imp) <- attributes(lpr)
    attr(lpr.imp,"source") <- NULL
    attr(lpr.imp,"imposed") <- TRUE
  }

  # Merge updates from ... and the updates argument.

  if (!is.list(updates))
  { stop("Updates argument is not a list")
  }

  original.updates <- updates <- c(list(...),updates)

  if (length(updates)==0)
  { stop("No MCMC updates specified")
  }

  # Go through the list of MCMC updates / restrictions, checking that they are
  # valid, putting them in a convenient form (updates always a list, indexes 
  # for restrictions sorted and in positive form), recording characteristics
  # of the updates, and for each update creating a wrappers for lpr if 
  # necessary.

  takes.lpr.initial <- rep (FALSE, length(updates))
  takes.initial.p <- rep (FALSE, length(updates))
  handles.bounds <- rep (FALSE, length(updates))
  special <- rep (FALSE, length(updates))
  wrapped.lpr <- rep (list(0), length(updates))
  fun_args <- rep(list(character(0)), length(updates))

  restrict.s <- if (is.list(initial)) names(initial) else NULL
  restrict.i <- NULL

  for (j in 1:length(updates))
  { 
    up <- updates[[j]]

    # Look for a specificaton of a restriction to part of the state.  
    # If found, check it for validity, put it in standard form, and
    # skip to next update / restriction.

    if (is.logical(up) && up || is.character(up) || is.numeric(up))
    {
      if (j==length(updates))
      { stop("list of updates must end with an actual update")
      }

      if (is.logical(up)) # will be TRUE
      { restrict.s <- if (is.list(initial)) names(initial) else NULL
        restrict.i <- NULL
      }
  
      else if (is.character(up))
      { if (j>1 && is.character(updates[[j-1]]))
        { up <- updates[[j]] <- c (updates[[j-1]], up)
        }
        if (!is.list(initial) || length(unique(up))!=length(up) 
                              || !all(up %in% names(initial)))
        { stop("Invalid string restriction")
        }
        restrict.s <- up
        restrict.i <- NULL
      }
  
      else if (is.numeric(up))
      { if (j>1 && is.numeric(updates[[j-1]]))
        { up <- updates[[j]] <- c (updates[[j-1]], up)
        }
        if (is.list(initial) && length(restrict.s)!=1)
        { stop(
          "Integer restriction not allowed with more than one element of state")
        }
        n <- length (if (is.list(initial)) initial[[restrict.s]] else initial)
        up <- updates[[j]] <- unique(sort(as.integer(up)))
        if (all(up<0))
        { if (any(-up>n))
          { stop("Index out of range in restriction")
          }
          up <- setdiff (1:n, -up)
        }
        else
        { if (any(up<1) || any(up>n))
          { stop("Index out of range in restriction")
          }
        }
        restrict.i <- up
      }

      next
    }

    # Check an update for validity, and put it in standard form.

    if (!is.list(up))
    { up <- updates[[j]] <- list(up)
    }
    if (length(up)<1 || !is.function (up[[1]]))
    { stop("MCMC update descriptions must start with a function")
    }
    if (length(formals(up[[1]]))<2 || names(formals(up[[1]]))[1] != "lpr"
                                   || names(formals(up[[1]]))[2] != "initial")
    { stop("update function doesn't have lpr and initial args")
    }

    # Look at named arguments of the update function, changing list arguments
    # to vectors according to restrict.s, with scalar elements expanded as
    # necessary, and subsetting vector arguments according to restrict.i.

    for (e in names(up))
    { if (e != "")
      { if (is.list(up[[e]]))
        { validate.mcmc.state (up[[e]], initial, allow.scalar=TRUE)
          if (!all (restrict.s %in% names(up[[e]])))
          { stop("list argument for update missing some elements being updated")
          }
          up[[e]] <- up[[e]][restrict.s]
          for (f in names(up[[e]]))
          { if (length(up[[e]][[f]])==1)
            { up[[e]][[f]] <- rep (up[[e]][[f]], length(initial[[f]]))
            }
          }
          up[[e]] <- unlist(up[[e]])
          if (!is.null(restrict.i)) up[[e]] <- up[[e]][restrict.i]
        }
        else if (is.numeric(up[[e]]) && length(up[[e]])!=1)
        { if (!is.list(initial) || length(restrict.s)!=1 ||
                length(up[[e]])!=length(unlist(initial[restrict.s])))
          { validate.mcmc.state (up[[e]], initial, allow.scalar=TRUE)
          }
          if (!is.null(restrict.i))
          { up[[e]] <- up[[e]][restrict.i]
          }
        }
        else if (is.function(up[[e]]))
        { fun_args[[j]] <- c(fun_args[[j]],e)
        }
        updates[[j]][[e]] <- up[[e]]
      }
    }

    # Record characteristics of the update.

    a <- attr(up[[1]],"special")
    special[j] <- !is.null(a) && a

    a <- attr(up[[1]],"handles.bounds")
    handles.bounds[j] <- !is.null(a) && a

    takes.lpr.initial[j] <- "lpr.initial" %in% names(formals(up[[1]]))

    takes.initial.p[j] <- 
       !is.null(initial.p) && "initial.p" %in% names(formals(up[[1]]))

    # Check that initial.p (if present) has either all or none of the
    # elements updated by this update.  Set takes.initial.p to FALSE
    # if none of the elements updated are in initial.p.

    if (takes.initial.p[j] && is.list(initial))
    { if (!any (restrict.s %in% names(initial.p)))
      { takes.initial.p[j] <- FALSE
      }
      else if (!all (restrict.s %in% names(initial.p)))
      { stop(
         "initial.p must have all or none of the elements used in an update")
      }
    }

    # Create a suitable wrapper for lpr for use by this update.

    if (!special[j])
    {
      this.lpr <- if (handles.bounds[j]) lpr else lpr.imp

      if (is.list(initial))  # state is a list
      { 
        if (length(restrict.s)==1) # state is a list restricted to one element
        { 
          wrapped.lpr[[j]] <- function  # formals replaced below, but like ...
                                (value, grad, ch.value, ch.pos, lpr.value) 
          { 
            if (lpr.has.ch.pos.arg && !missing(ch.pos))
            { if (is.null(restrict.i)) 
              { curr[[restrict.s]] <<- value
              }
              else
              { curr[[restrict.s]][restrict.i] <<- value
                ch.pos <- restrict.i[ch.pos]
              }
              if (lpr.has.grad.arg) 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               grad=grad, 
                               ch.value=ch.value, ch.elem=restrict.s, 
                               ch.pos=ch.pos, lpr.value=lpr.value)
              else 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               ch.value=ch.value, ch.elem=restrict.s, 
                               ch.pos=ch.pos, lpr.value=lpr.value) 
            }
            else if (lpr.has.ch.pos.arg && !is.null(restrict.i))
            { if (lpr.has.grad.arg) 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               grad=grad, 
                               ch.value=value, ch.elem=restrict.s, 
                               ch.pos=restrict.i, lpr.value=lpr.curr)
              else 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               ch.value=value, ch.elem=restrict.s, 
                               ch.pos=restrict.i, lpr.value=lpr.curr) 
            }
            else if (lpr.has.ch.elem.arg && is.null(restrict.i))
            { if (lpr.has.grad.arg) 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               grad=grad, ch.value=value, ch.elem=restrict.s,
                               lpr.value=lpr.curr)
              else
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               ch.value=value, ch.elem=restrict.s,
                               lpr.value=lpr.curr)
            }
            else
            { if (is.null(restrict.i))curr[[restrict.s]] <<- value
              else                    curr[[restrict.s]][restrict.i] <<- value
              if (lpr.has.grad.arg) 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                               grad=grad) 
              else 
                r <- this.lpr (if (is.null(fix)) curr else c(fix,curr))
            }
            if (!is.null(attr(r,"grad")))
            { attr(r,"grad") <- ( if (is.null(restrict.i))
                                    attr(r,"grad")[[restrict.s]]
                                  else 
                                    attr(r,"grad")[[restrict.s]][restrict.i] )
              if (length(value)!=length(attr(r,"grad")))
              { stop("lpr function didn't compute full gradient needed")
              }
            }
            r
          }
          formals(wrapped.lpr[[j]]) <- formals(this.lpr)
          formals(wrapped.lpr[[j]]) $ ch.elem <- NULL
          if (!lpr.has.ch.pos.arg)
          { formals(wrapped.lpr[[j]]) $ ch.value <- NULL
            formals(wrapped.lpr[[j]]) $ ch.pos <- NULL
            formals(wrapped.lpr[[j]]) $ lpr.value <- NULL
          }
        }
        else  # state is a list, and is not restricted to one element
        { 
          wrapped.lpr[[j]] <-  function # formals replaced below, but like ...
                                 (value, grad)
          { start <- 1
            for (e in restrict.s)
            { end <- start + length(curr[[e]]) - 1
              curr[[e]] <<- value[start:end]
              start <- end + 1
            }
            if (lpr.has.grad.arg) 
              r <- this.lpr (if (is.null(fix)) curr else c(fix,curr), 
                             grad=grad)
            else 
              r <- this.lpr (if (is.null(fix)) curr else c(fix,curr))
            if (!is.null(attr(r,"grad")))
            { attr(r,"grad") <- 
                unlist (attr(r,"grad")[restrict.s], use.names=FALSE)
              if (length(value)!=length(attr(r,"grad")))
              { stop("lpr function didn't compute full gradient needed")
              }
            }
            r
          }
          formals(wrapped.lpr[[j]]) <- formals(this.lpr)
          formals(wrapped.lpr[[j]]) $ ch.value <- NULL
          formals(wrapped.lpr[[j]]) $ ch.elem <- NULL
          formals(wrapped.lpr[[j]]) $ ch.pos <- NULL
          formals(wrapped.lpr[[j]]) $ lpr.value <- NULL
        }

        attributes (wrapped.lpr[[j]]) <- attributes (this.lpr)
        attr (wrapped.lpr[[j]], "source") <- NULL

        for (b in c("lower","upper"))
        { a <- attr (this.lpr, b)
          if (!is.null(a))
          { a <- unlist (a [intersect(names(a),restrict.s)])
            if (!is.null(restrict.i)) a <- if (length(a)==1) a
                                           else a[restrict.i]
            attr (wrapped.lpr[[j]], b) <- a
          }
        }
      }
      else  # state is a vector
      {
        if (is.null(restrict.i))  # state is a vector, and not restricted
        { 
          wrapped.lpr[[j]] <- this.lpr
        }
        else  # state is a vector with a restriction to a subset of indexes
        { 
          wrapped.lpr[[j]] <- function  # formals replaced below, but like ...
                                (value, grad, ch.value, ch.pos, lpr.value) 
          { if (lpr.has.ch.pos.arg)
            { if (missing(ch.pos))
              { if (lpr.has.grad.arg)
                  r <- this.lpr (cv, grad=grad, ch.value=value, 
                                     ch.pos=restrict.i, lpr.value=lpr.curr)
                else
                  r <- this.lpr (cv, ch.value=value, 
                                     ch.pos=restrict.i, lpr.value=lpr.curr)
              }
              else
              { cv[restrict.i] <<- value
                ch.pos <- restrict.i[ch.pos]
                if (lpr.has.grad.arg)
                  r <- this.lpr (cv, grad=grad, ch.value=ch.value, 
                                     ch.pos=ch.pos, lpr.value=lpr.value)
                else
                  r <- this.lpr (cv, ch.value=ch.value, 
                                     ch.pos=ch.pos, lpr.value=lpr.value)
              }
            }
            else
            { cv[restrict.i] <<- value
              if (lpr.has.grad.arg)
                r <- this.lpr (cv, grad=grad) 
              else 
                r <- this.lpr (cv)
            }
            if (!is.null(attr(r,"grad")))
            { attr(r,"grad") <- attr(r,"grad")[restrict.i]
              if (length(value)!=length(attr(r,"grad")))
              { stop("lpr function didn't compute full gradient needed")
              }
            }
            r
          }
          formals(wrapped.lpr[[j]]) <- formals(this.lpr)
          attributes(wrapped.lpr[[j]]) <- attributes(this.lpr)
          attr (wrapped.lpr[[j]], "source") <- NULL

          if (!is.null(restrict.i))
          { for (b in c("lower","upper"))
            { a <- attr (this.lpr, b)
              attr (wrapped.lpr[[j]], b) <- if (length(a)==1) a
                                            else a[restrict.i]
            }
          }
        }
      }
    }
  }

  # Allocate space to store the results in q and (possibly) p.

  if (is.list(initial))
  { q <- list()
    for (e in names(initial))
    { q[[e]] <- matrix (NA, iterations, length(initial[[e]]), 
                        dimnames = list(NULL,names(initial[[e]])))
    }
    if (!is.null(initial.p))
    { p <- q[names(initial.p)]
    }
  }
  else
  { q <- matrix (NA, iterations, length(initial), 
                 dimnames = list(NULL,names(initial)))
    if (!is.null(initial.p))
    { p <- q
    }
  }

  # Set the state of the random number generator.

  if (is.null(saved.rand))
  { set.seed(seed)
  }
  else
  { if (!is.environment(saved.rand) || is.null(saved.rand$saved.state))
    { stop("Invalid saved.rand argument")
    }
    .Random.seed <- saved.rand$saved.state
  }

  # Set up NULL matrices of information recorded.  Replaced once number of
  # columns is known, for those that are used.

  last.mat <- NULL
  rec.mat <- NULL
  ave.mat <- NULL

  # Set the current state to the initial state.

  curr <- initial
  curr.p <- initial.p

  lpr.curr <- lpr.initial
  if (is.null(lpr.curr))
  { lpr.curr <- lpr.imp (if (is.null(fix)) curr else c(fix,curr))
  }
  names(lpr.curr) <- NULL
  dim(lpr.curr) <- NULL

  # Do the specified number of iterations.

  for (i in 1:iterations)
  { 
    # Repeat the sequence of updates the specified number of times.

    for (k in 1:rep)
    { 
      # No restriction of state at the start.

      restrict.s <- if (is.list(initial)) names(initial) else NULL
      restrict.i <- NULL

      # Loop to do updates, or note restrictions.

      this.rec <- NULL
      this.ave.p <- NULL
      u <- 1

      for (j in 1:length(updates))
      { 
        up <- updates[[j]]

        # Look for a specification of a restriction to part of the state.
        # If found, note it and skip to the next update / restriction.

        if (is.logical(up)) # will be TRUE
        { restrict.s <- if (is.list(initial)) names(initial) else NULL
          restrict.i <- NULL
          next
        }
        if (is.character(up))
        { restrict.s <- up
          restrict.i <- NULL
          next
        }
        if (is.numeric(up))
        { restrict.i <- up
          next
        }

        # Forget saved gradient when restriction may have changed.

        if (j==1 || !is.list(updates[[j-1]]))
        {  attr(lpr.curr,"grad") <- NULL
        }

        # Set the lpr function for a full state, depending on whether this
        # update function handles bounds itself.  This is accessed by the
        # wrapped.lpr function, or passed to a special update function.

        this.lpr <- if (handles.bounds[j]) lpr else lpr.imp

        # Set up arguments for the call of a special or general-purpose
        # update function.

        if (special[j]) # special-purpose update gets original state &this.lpr
        { 
          if (takes.initial.p[j])
          { args <- list (this.lpr,
                          initial = if (is.null(fix)) curr else c(fix,curr),
                          initial.p = curr.p)
          }
          else
          { args <- list (this.lpr,
                          initial = if (is.null(fix)) curr else c(fix,curr))
          }

          up.args <- up[-1]
        }
        else  # general-purpose update gets state as a vector & wrapped lpr
        {
          # Extract the part of the state being updated as a vector.

          if (!is.list(curr))
          { cv <- curr
            if (takes.initial.p[j])
            { cv.p <- curr.p
            }
          }
          else if (length(restrict.s)==1)
          { cv <- curr[[restrict.s]]
            if (takes.initial.p[j]) 
            { cv.p <- curr.p[[restrict.s]]
            }
          }
          else
          { cv <- unlist (curr[restrict.s], use.names=FALSE)
            if (takes.initial.p[j]) 
            { cv.p <- unlist (curr.p[restrict.s], use.names=FALSE)
            }
          }

          # Set up first part of argument list for update function.

          if (takes.initial.p[j])
          { if (is.null(restrict.i))
            { args <- list (wrapped.lpr[[j]], initial = cv, 
                                              initial.p = cv.p)
            }
            else
            { args <- list (wrapped.lpr[[j]], initial = cv[restrict.i], 
                                              initial.p = cv.p[restrict.i])
            }
          }
          else
          { if (is.null(restrict.i))
            { args <- list (wrapped.lpr[[j]], initial = cv)
            }
            else
            { args <- list (wrapped.lpr[[j]], initial = cv[restrict.i])
            }
          }

          # Set up addtional arguments for update function, calling functions
          # as necessary.

          up.args <- up[-1]
    
          for (e in fun_args[[j]])
          { if (is.list(initial))
            { curr.other <- curr
              if (is.null(restrict.i))
              { curr.other[restrict.s] <- NULL
              }
              else
              { curr.other[[restrict.s]][restrict.i] <- NA
              }
              up.args[[e]] <- do.call (up.args[[e]],list(c(fix,curr.other)))
              if (is.list(up.args[[e]]))
              { up.args[[e]] <- 
                  unlist (up.args[[e]][restrict.s], use.names=FALSE)
              }
            }
            else  # state is a vector
            { if (is.null(restrict.i))
              { up.args[[e]] <- do.call (up.args[[e]],list())
              }
              else
              { curr.other <- curr
                curr.other[restrict.i] <- NA
                up.args[[e]] <- do.call (up.args[[e]],list(curr.other))
              }
            }
            if (is.list(up.args[[e]]) || is.numeric(up.args[[e]]) &&
                 length(up.args[[e]])!=1 && length(up.args[[e]])!=length(cv))
            { stop("argument for update obtained by function call is invalid")
            }
          }
        }

        if (takes.lpr.initial[j])
        { args$lpr.initial <- lpr.curr
        }

        # Call the update function.

        r <- eval(as.call(c(up[1],args,up.args)))

        if (takes.initial.p[j] && is.null(r$final.p))
        { stop("update taking initial.p didn't return final.p")
        }

        # Change all or part of the state based on the result of the update.
      
        if (special[j])  # special update, with full state
        {  
           curr <- r$final [names(initial)]
           if (takes.initial.p[j])
           { curr.p <- r$final.p [names(initial.p)]
           }
        }
        else  # general-purpose update, with subset of state as vector
        { 
          # Replace part of vector according to integer restriction.

          if (is.null(restrict.i))
          { cv <- r$final
            if (takes.initial.p[j])
            { cv.p <- r$final.p
            }
          }
          else
          { cv[restrict.i] <- r$final
            if (takes.initial.p[j])
            { cv.p[restrict.i] <- r$final.p
            }
          }

          # Replace parts of state according to string restriction.

          if (is.list(initial))  # state is a list
          { 
            if (length(restrict.s)==1)
            { curr[[restrict.s]] <- cv
              if (takes.initial.p[j])
              { curr.p[[restrict.s]] <- cv.p
              }
            }
            else
            { start <- 1
              for (e in restrict.s)
              { end <- start + length(curr[[e]]) - 1
                curr[[e]] <- cv[start:end]
                if (takes.initial.p[j])
                { curr.p[[e]] <- cv.p[start:end]
                }
                start <- end + 1
              }
            }
          }
          else  # state is a vector
          { 
            curr <- cv
            if (takes.initial.p[j])
            { curr.p <- cv.p
            }
          }
        }

        # Compute the lpr value if it was not returned in the update result.

        if (is.null(r$lpr))
        { r$lpr <- lpr.imp (if (is.null(fix)) curr else c(fix,curr))
        }

        names(r$lpr) <- NULL
        dim(r$lpr) <- NULL

        lpr.curr <- r$lpr

        # Record rec and ave information from the update, as requested.

        if (!is.null(rec) && k==rep)
        { info <- unlist(r[intersect(rec,names(r))])
          if (length(info)>0)
          { names(info) <- paste(names(info),".",u,sep="")
            this.rec <- if (is.null(this.rec)) info else c(this.rec,info)
          }
        }

        if (!is.null(ave))
        { info <- unlist(r[intersect(ave,names(r))])
          if (length(info)>0)
          { names(info) <- paste(names(info),".",u,sep="")
            this.ave.p <- if (is.null(this.ave.p)) info 
                          else c(this.ave.p,info)
          }
        }

        u <- u + 1
      }

      # Convert list of ave information to vector, and accumulate averages.

      if (!is.null(ave))
      { if (k==1)
        { this.ave <- this.ave.p
        }
        else
        { if (length(this.ave)!=length(this.ave.p)
           || any(names(this.ave)!=names(this.ave.p)))
          { stop("ave quantities not consistently returned by updates")
          }
          this.ave <- this.ave + this.ave.p
        }
      }
    }

    # Record information for this iteration.

    if (!is.null(last))
    { this.last <- unlist (r[intersect(last,names(r))])
      if (i==1)
      { last.mat <- matrix (NA, iterations, length(this.last))
        colnames(last.mat) <- names(this.last)
      }
      if (ncol(last.mat)!=length(this.last) 
       || any(colnames(last.mat)!=names(this.last)))
      { stop("last quantities not consistently returned by updates")
      }
      last.mat[i,] <- this.last
    }

    if (!is.null(rec))
    { this.rec   <- unlist (this.rec)
      if (i==1)
      { rec.mat <- matrix (NA, iterations, length(this.rec))
        colnames(rec.mat) <- names(this.rec)
      }
      if (ncol(rec.mat)!=length(this.rec) 
       || any(colnames(rec.mat)!=names(this.rec)))
      { stop("rec quantities not consistently returned by updates")
      }
      rec.mat[i,] <- this.rec
    }

    if (!is.null(ave))
    { if (i==1)
      { ave.mat <- matrix (NA, iterations, length(this.ave))
        colnames(ave.mat) <- names(this.ave)
      }
      if (ncol(ave.mat)!=length(this.ave) 
       || any(colnames(ave.mat)!=names(this.ave)))
      { stop("ave quantities not consistently returned by updates")
      }
      ave.mat[i,] <- this.ave / rep
    }

    # Store the new state in q and possibly p.

    if (is.list(initial))
    { for (e in names(initial))
      { q[[e]][i,] <- curr[[e]]
      }
      if (!is.null(initial.p)) 
      { for (e in names(initial.p))
        { p[[e]][i,] <- curr.p[[e]] 
        }
      }
    }
    else
    { q[i,] <- curr 
      if (!is.null(initial.p)) p[i,] <- curr.p
    }
  }

  # Return the results.

  if (is.list(curr))
  { for (e in names(curr))
    { names(curr[[e]]) <- colnames(q[[e]])
    }
  }
  else
  { names(curr) <- colnames(q)
  }

  res <- list (final = curr, q = q, lpr.final = lpr.curr, 
               last = last.mat, rec = rec.mat, ave = ave.mat,
               lpr.arg = lpr, initial.arg = initial.arg, rep.arg = rep, 
               last.arg = last, rec.arg = rec, ave.arg = ave,
               saved.rand = as.environment(list(saved.state=.Random.seed)), 
               seed.arg = seed, fix.arg = fix, updates = original.updates)

  if (!is.null(initial.p))
  { if (is.list(curr))
    { for (e in names(curr.p))
      { names(curr.p[[e]]) <- colnames(p[[e]])
      }
    }
    else
    { names(curr.p) <- colnames(p)
    }
    res$final.p <- curr.p
    res$p <- p
  }

  res
}


# DO MORE MCMC ITERATIONS.  See mcmc.txt for documentation.

more_mcmc <- function (previous, iterations)
{
  new <- mcmc (previous$lpr.arg, 
               previous$final, 
               iterations, 
               updates = previous$updates, 
               rep = previous$rep.arg,
               fix = previous$fix.arg, 
               initial.p = previous$final.p,
               last = previous$last.arg, 
               rec = previous$rec.arg,
               ave = previous$ave.arg, 
               lpr.initial = previous$lpr.final,
               initial.arg = previous$initial.arg,
               seed = previous$seed.arg,
               saved.rand = previous$saved.rand)

  new$q <- rbind (previous$q, new$q)

  if (!is.null(new$p))    new$p    <- rbind (previous$p, new$p)
  if (!is.null(new$last)) new$last <- rbind (previous$last, new$last)
  if (!is.null(new$rec))  new$rec  <- rbind (previous$rec, new$rec)
  if (!is.null(new$ave))  new$ave  <- rbind (previous$ave, new$ave)

  new
}


# CHECK THAT A MCMC STATE (OR RELATED STRUCTURE) IS VALID.  Verifies
# that s is a numeric vector or list of numeric vectors, has suitable
# names, and if m is present, is a subset of m, with matching dimensions
# of elements.  If allow.scalar is TRUE, a scalar element in s can match 
# a vector in m.  If allow.infinite is TRUE, the numerical values in s 
# need not be finite.

validate.mcmc.state <- function (s, m, allow.scalar=FALSE, allow.infinite=FALSE)
{
  lengths.match <- function (s,m)
    length(s)==1 && allow.scalar || length(s)==length(m)

  if (is.numeric(s))
  { if (!missing(m))
    { if (!is.numeric(m) || !lengths.match(s,m))
      { stop("mismatched state values")
      }
    }
  }

  else if (is.list(s))
  { 
    n <- names(s)
    if (is.null(n) || any(n==""))
    { stop("all elements in a list of state values must have names")
    }

    if (!missing(m))
    { if (!is.list(m))
      { stop("mismatched state values")
      }
      if (any(!is.element(n,names(m))))
      { stop("mismatched state values")
      }
    }

    for (e in n)
    { if (!is.numeric(s[[e]]))
      { stop("invalid state value")
      }
      if (!missing(m) && !lengths.match(s[[e]],m[[e]]))
      { stop("mismatched state values")
      }
    }
  }

  else 
  { stop("invalid state value")
  }

  if (!allow.infinite && !all(is.finite(unlist(s))))
  { stop("some state values are not finite")
  }
}
