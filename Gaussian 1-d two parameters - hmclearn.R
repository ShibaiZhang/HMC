# install.packages("hmclearn")
library(hmclearn)


U <- function(theta, X) {
  mu <- theta[1]
  sigma <-  exp(theta[2])
  log_prob <- sum( dnorm(X, mu, sigma, log=TRUE) )
  return(log_prob)
}

# gradient function
# need vector of partial derivatives of U with respect to vector theta
grad_U <- function(theta, X) {
  mu <- theta[1]
  sigma <- exp(theta[2])
  n = length(X)
  d_mu = sum( X - mu )* exp(-2*theta[2])
  d_sigma =  sum( (X - mu)^2 ) * exp(-2*theta[2])- length(X) # log helps with more stable calculation
  return( c(d_mu, d_sigma) ) # negative bc energy is neg-log-prob
}



true_mu = 0
true_sigma = 1
X <- rnorm(100, true_mu, true_sigma)


f <- hmc(N = 100,
         theta.init = c(true_mu+0.1, log(true_sigma)+0.1),
         epsilon = 0.01,
         L = 25,
         logPOSTERIOR = U,
         glogPOSTERIOR = grad_U,
         varnames = c("mu","sigma"),
         param=list(X=X), parallel=FALSE, chains=1)

f$thetaCombined
f$accept
f$N


summary(f)


#diagplots(f, burnin=300)