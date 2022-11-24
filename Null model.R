library(geosphere) # for bearing


HMC = function (U, grad_U, epsilon, L, current_q, X_distance, X_bearing_diff)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q, X_distance, X_bearing_diff) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q, X_distance, X_bearing_diff)
    print(p)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q, X_distance, X_bearing_diff) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # print(p)
  # print(q)
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q, X_distance, X_bearing_diff)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q, X_distance, X_bearing_diff)
  proposed_K = sum(p^2) / 2
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

# function transfer locations to distances
distance <- function(X){
  dist_all_step <-dist(X,p=2)
  dist_matrix <- as.matrix(dist_all_step)
  dist_neighbor_step <- diag(dist_matrix[2:nrow(dist_matrix),1:ncol(dist_matrix)-1])
  dist_neighbor_step
}


# function transfer locations to directions
direction <- function(X){
  angle <- (bearing(X)/180)*pi
  angle[-length(angle)]
}



# U: neg log posterior
U <- function(q, X_distance, X_bearing_diff) {
  kappa <- exp(q[1])
  lambda <- exp(q[2])
  rho <- exp(q[3])
  log_prob <- sum(dweibull(X_distance,kappa,lambda,log=TRUE)) + 
    sum(log(sinh(rho))-log(cosh(rho)-cos(X_bearing_diff)))
  
  return(-log_prob)
}


# gradient function
grad_U <- function(q, X_distance, X_bearing_diff) {
  kappa <- exp(q[1])
  lambda <- exp(q[2])
  rho <- exp(q[3])
  n = length(X_distance)
  #d_lambda = -n*kappa/lambda + kappa*sum((X_distance^kappa)/(lambda)^(kappa+1))
  d_lambda = -n*exp(q[1]-q[2]) + kappa*sum(X_distance^kappa)*exp(-q[2]*(kappa+1))
  d_kappa = n/kappa - n*log(lambda)-sum((log(X_distance/lambda))*exp(kappa*(log(X_distance/lambda)))) +sum(log(X_distance)) 
  d_rho = n*cosh(rho)/sinh(rho) - sum(sinh(rho)/(cosh(rho)-cos(X_bearing_diff)))

  return( c(-d_kappa, -d_lambda, -d_rho) )
}


X <- locationAll[,1:2]

X_distance <- distance(X)
X_distance <- X_distance[-1]
X_bearing <- direction(X)
X_bearing_diff <- X_bearing[-1] - X_bearing[-length(X_bearing)]

init_kappa = 1.5 #dunif(1, 0, 5)
init_lambda = 4 #dunif(1, 0, 10)
init_rho = 1.5 #dunif(1, 0, 2)

n_step = 1000
q_history <- matrix(nrow=n_step, ncol = 3)
q_history[1,] <- c(log(init_kappa), log(init_lambda), log(init_rho)) #kappa, lambda, rho
L <- 40 
epsilon <- 0.0005

for ( i in 2:n_step ) {
  q_history[i,] <- HMC(U, grad_U, epsilon, L, q_history[i-1,], X_distance, X_bearing_diff)
}

samples = exp(q_history[(0.1*n_step):n_step, ])
effective_sample = unique(samples)
effective_sample_size = nrow(effective_sample)
accept_rate = effective_sample_size / nrow(samples)
accept_rate

lambda_est = samples[,1]
kappa_est = samples[,2]
rho_est = samples[,3]

# estimated lambda
mean(lambda_est)
# estimated kappa
mean(kappa_est)
# estimated rho
mean(rho_est)

hist(lambda_est, xlab='lambda', main='Histogram of estimated lambda')
hist(kappa_est, xlab='kappa', main='Histogram of estimated kappa')
hist(rho_est, xlab='rho', main='Histogram of estimated rho')