# install.packages("hmclearn")
library(hmclearn)
library(MASS)

U <- function(theta, X, sigma, det_sigma, inverse_sigma) {
  # U: neg log posterior
  mul <- function(x){
    t(x - theta) %*% inverse_sigma %*% (x-theta)
  }
  
  n = length(X) 
  log_prob <-  -n*log(2*pi) - n/2*log(det_sigma) - 0.5*sum(apply(X, 1, mul))
  return(log_prob)
}


# gradient function
# need vector of partial derivatives of U with respect to vector theta


grad_U <- function(theta, X, sigma, det_sigma, inverse_sigma) {
  calc_d_mu1 <- function(x){
    (x[1] - theta[1])*sigma[2,2] - (x[2] - theta[2]) * sigma[1,2]
  }
  
  calc_d_mu2 <- function(x){
    (x[2] - theta[2])*sigma[1,1] - (x[1] - theta[1]) * sigma[2,1]
  }
  
  d_mu1 = sum(apply(X, 1, calc_d_mu1))
  d_mu2 = sum(apply(X, 1, calc_d_mu2))
  return( c(-d_mu1, -d_mu2) ) 
}



true_mu = c(0,0)
sigma = matrix(c(1, 0.95, 0.95, 1), nrow=2, ncol=2)
inverse_sigma = solve(sigma)
det_sigma = det(sigma)
X = mvrnorm(100, true_mu, sigma)

f <- hmc(N = 1000,
         theta.init = true_mu+0.1,
         epsilon = 0.005,
         L = 20,
         logPOSTERIOR = U,
         glogPOSTERIOR = grad_U,
         varnames = c("mu1","mu2"),
         param=list(X=X, sigma=sigma, det_sigma=det_sigma, inverse_sigma=inverse_sigma), 
         parallel=FALSE, chains=1)

f$thetaCombined
f$accept
f$N


summary(f, burnin=0.1*f$N)


diagplots(f, burnin=100)
