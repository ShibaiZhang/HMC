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


# X ~ normal(mu, 1)


U <- function(q, X) {
  mu <- q
  sigma <- 1
  log_prob <- sum( dnorm(X, mu, sigma, log=TRUE) ) #+ dnorm(q, 5, 1, log=TRUE)
  return(-log_prob)
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
grad_U <- function(q, X) {
  grad = sum( X - q ) #+ (5 - q)/1^2
  return( -grad ) # negative bc energy is neg-log-prob
}



# main
true_mu = 5
X <- rnorm(100, true_mu)

n_step = 10000
q_history <- rep( NA , n_step )
q_history[1] <- 5.5
L <- 20 
epsilon <- 0.01 # 0.01 / 0.03

for ( i in 2:n_step ) {
  q_history[i] <- HMC( U , grad_U , epsilon , L , q_history[i-1], X )
}

q_history[(0.1*n_step):n_step]
mean(q_history[(0.1*n_step):n_step])

par("mar")
par(mar=c(1,1,1,1))

times <- 1:1000000

df <- data.frame(times = times, q_history = q_history)
ggplot(data = df, mapping = aes(x = times, y = q_history)) + geom_line()