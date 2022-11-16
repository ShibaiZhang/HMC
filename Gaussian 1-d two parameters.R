HMC = function (U, grad_U, epsilon, L, current_q, X)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q, X) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q, X)
    
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q, X) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q, X)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q, X)
  proposed_K = sum(p^2) / 2
  # print('inital p')
  # print(current_p)
  # print('p after leapfrog')
  # print(p)
  # print(paste(current_U,proposed_U,current_K,proposed_K))
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
}



# U = -log(pi(q))
# X ~ normal(mu, sigma)

U <- function(q, X) {
  mu <- q[1]
  sigma <-  exp(q[2])
  log_prob <- sum( dnorm(X, mu, sigma, log=TRUE) )
  return(-log_prob)
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
grad_U <- function(q, X) {
  mu <- q[1]
  sigma <- exp(q[2])
  n = length(X)
  d_mu = sum( X - mu )* exp(-2*q[2])
  d_sigma =  sum( (X - mu)^2 ) * exp(-2*q[2])- length(X) # log helps with more stable calculation
  return( c(-d_mu, -d_sigma) ) # negative bc energy is neg-log-prob
}



true_mu = 10
true_sigma = 10
X <- rnorm(100, true_mu, true_sigma)

n_step = 100
q_history <- matrix(nrow=n_step, ncol = 2)
# use log(sigma) instead of sigma to avoid numeric problems
q_history[1,] <- c(mean(X)+0.1, log(var(X)^0.5)+0.1)
L <- 20 
epsilon <- 0.01 # 0.01 / 0.03

for ( i in 2:n_step ) {
  q_history[i, ] <- HMC( U , grad_U , epsilon , L , q_history[i-1,], X )
}

q_history
# estimated mu
mean(q_history[(0.1*n_step):n_step,1])
# estimated sigma
mean(exp(q_history[(0.1*n_step):n_step,2]))
