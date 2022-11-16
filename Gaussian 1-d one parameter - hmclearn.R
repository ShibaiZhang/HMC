# install.packages("hmclearn")
library(hmclearn)

U <- function(theta, X) {
  mu <- theta
  sigma <- 1
  log_prob <- sum( dnorm(X, mu, sigma, log=TRUE) ) 
  return(log_prob)
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
grad_U <- function(theta, X) {
  grad = sum( X - theta ) #+ (5 - q)/1^2
  return( grad ) # negative bc energy is neg-log-prob
}




true_mu = 0
X <- rnorm(100, true_mu)



f <- hmc(N = 100,
         theta.init = true_mu+0.1,
         epsilon = 0.01,
         L = 25,
         logPOSTERIOR = U,
         glogPOSTERIOR = grad_U,
         varnames = "mu",
         param=list(X=X), parallel=FALSE, chains=1)

f$thetaCombined
f$accept
f$N


summary(f)


#diagplots(f, burnin=300)