library(shape)# for good arrow heads


HMC <- function (U, grad_U, epsilon, L, current_q , ... ) {
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U( q , ... ) / 2
  # Alternate full steps for position and momentum
  qtraj <- matrix(NA,nrow=L+1,ncol=length(q))
  ptraj <- qtraj
  qtraj[1,] <- current_q
  ptraj[1,] <- p
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) {
      p = p - epsilon * grad_U(q,...)
      ptraj[i+1,] <- p
    }
    qtraj[i+1,] <- q
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q,...) / 2
  ptraj[L+1,] <- p
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q,...)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q,...)
  proposed_K = sum(p^2) / 2
  H0 <- current_U + current_K
  H1 <- proposed_U + proposed_K
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  new_q <- q
  accept <- 0
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    new_q <- q  # accept
    accept <- 1
  } else
    new_q <- current_q  # reject
  
  return(
    list(
      q = new_q,
      traj = qtraj,
      ptraj = ptraj,
      accept = accept , 
      dH = H1 - H0 ) )
}


#U(q,X)
U <- function(q,X) {
  kappa <- exp(q[1])
  lambda <- exp(q[2])
  log_prob <- sum(dweibull(X,kappa,lambda,log=TRUE))
  
  return(-log_prob)
}


# gradient function
grad_U <- function(q,X) {
  kappa <- exp(q[1])
  lambda <- exp(q[2])
  n = length(X)
  d_lambda = -n*kappa/lambda + kappa*sum((X^kappa)/(lambda)^(kappa+1))
  d_kappa = n/kappa - n*log(lambda)-sum((log(X/lambda))*exp(kappa*(log(X/lambda)))) +sum(log(X)) 
  return( c(-d_kappa, -d_lambda) ) # remember that 
}

#main
true_kappa = 1.6
true_lambda = 5.5
X <- rweibull(1200, true_kappa, true_lambda)
n_step = 100
q_history <- matrix(nrow=n_step, ncol = 2)
q_history[1,] <- c(log(true_kappa+0.1), log(true_lambda+0.1))
L <- 35 
epsilon <- 0.001 

# ==== plot =====
rad <- 0.1
xr <- c(log(true_kappa)-rad,log(true_kappa)+rad)
yr <- c(log(true_lambda)-rad,log(true_lambda)+rad)
plot( NULL , ylab="log_kappa" , xlab="log_lambda" , xlim=xr , ylim=yr )
xpos <- c(4,2,1,2) # for L=20
# ===============


for ( i in 2:n_step ) {
  Q <- HMC( U , grad_U , epsilon , L , q_history[i-1,], X )
  q_history[i,] <- Q$q
  
  if ( i < 20 ) {
    lines( Q$traj , lwd=2)
    points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )
    Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] , arr.length=0.3 , arr.adj = 0.7 )
    #text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=0.8 , pos=xpos[i] , offset=0.4 )
  }
  points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 16 , 1 ) , 
          col=ifelse( abs(Q$dH)>0.1 , "red" , "black" ) )
  
}

cb <- 0.2
log_kappa_seq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=100) 
log_lambda_seq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=100)
z <- matrix(NA,length(log_kappa_seq),length(log_lambda_seq))
for ( i in 1:length(log_kappa_seq) )
  for ( j in 1:length(log_lambda_seq) )
    z[i,j] <- U( c(log_kappa_seq[i] ,log_lambda_seq[j]), X)
cl <- contourLines( log_kappa_seq , log_lambda_seq , z , nlevels=51 )
for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y ,  lwd=0.5 )
mtext( paste("Leapfrog steps = ",L) )


samples = q_history[(0.1*n_step):n_step]
# estimated lambda
exp(mean(q_history[(0.1*n_step):n_step,1]))
# estimated kappa
exp(mean(q_history[(0.1*n_step):n_step,2]))