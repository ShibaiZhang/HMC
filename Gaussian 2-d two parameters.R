HMC = function (U, grad_U, epsilon, L, current_q, X)
{
  q = current_q
  #print(q)
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
    # print(paste("q =",round(q, 4)))
    # print(paste("p =",round(p, 4)))
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
  # print(current_p)
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


library(MASS) # for mvrnorm()

# U: neg log posterior
mul <- function(x, q, inverse_sigma){
  t(x - q) %*% inverse_sigma %*% (x-q)
}

U <- function(q, X) {
  n = length(X) 
  log_prob <-  -n*log(2*pi) - n/2*log(det_sigma) - 0.5*sum(apply(X, 1, mul, q, inverse_sigma))
  return(-log_prob)
}


# gradient function
# need vector of partial derivatives of U with respect to vector q

calc_d_mu1 <- function(x,q){
  (x[1] - q[1])*sigma[2,2] - (x[2] - q[2]) * sigma[1,2]
}

calc_d_mu2 <- function(x, q){
  (x[2] - q[2])*sigma[1,1] - (x[1] - q[1]) * sigma[2,1]
}

grad_U <- function(q, X) {
  d_mu1 = sum(apply(X, 1, calc_d_mu1, q))
  d_mu2 = sum(apply(X, 1, calc_d_mu2, q))
  return( c(-d_mu1, -d_mu2) ) 
}


true_mu = c(2,5)
sigma = matrix(c(1, 0.95, 0.95, 1), nrow=2, ncol=2)
inverse_sigma = solve(sigma)
det_sigma = det(sigma)
X = mvrnorm(100, true_mu, sigma)

n_step = 1000
q_history <- matrix(nrow=n_step, ncol = 2)
q_history[1,] <- true_mu + 0.01
L <- 25 
epsilon <- 0.01 # 0.01 / 0.03


for ( i in 2:n_step ) {
  q_history[i, ] <- HMC( U , grad_U , epsilon , L , q_history[i-1,], X )
}

q_history
# estimated mu1
mean(q_history[(0.1*n_step):n_step,1])
# estimated mu2
mean(q_history[(0.1*n_step):n_step,2])
