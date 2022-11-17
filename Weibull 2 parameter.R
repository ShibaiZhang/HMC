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
    #print(p)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q, X) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # print(p)
  # print(q)
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



#U(q,X)
U <- function(q,X) {
  kappa <- exp(q[1])
  lamda <- exp(q[2])
  log_prob <- sum(dweibull(X,kappa,lamda,log=TRUE))
  
  return(-log_prob)
}


# gradient function
grad_U <- function(q,X) {
  kappa <- exp(q[1])
  lamda <- exp(q[2])
  n = length(X)
  d_lamda = -n*kappa/lamda + kappa*sum((X^kappa)/(lamda)^(kappa+1))
  d_kappa = n/kappa - n*log(lamda)-sum((log(X/lamda))*exp(kappa*(log(X/lamda)))) +sum(log(X)) 
  return( c(-d_kappa, -d_lamda) ) # remember that 
}

#main
true_lamda = 1
true_kappa = 1
X <- rweibull(100, true_lamda, true_kappa)
n_step = 10000
q_history <- matrix(nrow=n_step, ncol = 2)
q_history[1,] <- c(log(true_lamda+0.1), log(true_kappa+0.1))
L <- 20 
epsilon <- 0.01 # 0.01 / 0.03

for ( i in 2:n_step ) {
  q_history[i,] <- HMC( U , grad_U , epsilon , L , q_history[i-1,], X )
}

q_history[(0.1*n_step):n_step]
# estimated lambda
exp(mean(q_history[(0.1*n_step):n_step,1]))
# estimated kappa
exp(mean(q_history[(0.1*n_step):n_step,2]))