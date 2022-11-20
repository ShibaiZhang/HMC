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
init_lambda = 5.3 #dunif(1, 0, 10)
init_rho = 1.2 #dunif(1, 0, 2)

n_step = 100
q_history <- matrix(nrow=n_step, ncol = 3)
q_history[1,] <- c(log(init_kappa), log(init_lambda), log(init_rho)) #kappa, lambda, rho
L <- 20 
epsilon <- 0.001 

# ==== plot =====
rad <- 0.4
xr <- c(log(init_kappa)-rad,log(init_kappa)+rad)
yr <- c(log(init_lambda)-rad,log(init_lambda)+rad)
plot( NULL , ylab="log_kappa" , xlab="log_lambda" , xlim=xr , ylim=yr )
xpos <- c(4,2,1,2) # for L=20
# ===============

for ( i in 2:n_step ) {
  
  Q <- HMC(U, grad_U, epsilon, L, q_history[i-1,], X_distance, X_bearing_diff)
  q_history[i,] <- Q$q
  
  if ( i < 10 ) {
    lines( Q$traj , lwd=2)
    points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )
    Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] , arr.length=0.3 , arr.adj = 0.7 )
    #text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=0.8 , pos=xpos[i] , offset=0.4 )
  }
  points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 16 , 1 ) , 
          col=ifelse( abs(Q$dH)>0.1 , "red" , "black" ) )
  
}

#=== compute contour of log-prob ===

U_weibull <- function(q,X) {
  kappa <- exp(q[1])
  lambda <- exp(q[2])
  log_prob <- sum(dweibull(X,kappa,lambda,log=TRUE))
  return(-log_prob)
}

cb <- 0.2
log_kappa_seq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=100) 
log_lambda_seq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=100)
z <- matrix(NA,length(log_kappa_seq),length(log_lambda_seq))
for ( i in 1:length(log_kappa_seq) )
  for ( j in 1:length(log_lambda_seq) )
    z[i,j] <- U_weibull( c(log_kappa_seq[i] ,log_lambda_seq[j]), X_distance)
cl <- contourLines( log_kappa_seq , log_lambda_seq , z , nlevels=21 )
for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y ,  lwd=0.5 )
mtext( paste("Leapfrog steps = ",L) )
#===================================


samples = exp(q_history[(0.1*n_step):n_step, ])
effective_sample = unique(samples)
effective_sample_size = length(effective_sample)
accept_rate = effective_sample_size / length(samples)
accept_rate

# estimated lambda
mean(samples[,1])
# estimated kappa
mean(samples[,2])
# estimated rho
mean(samples[,3])





