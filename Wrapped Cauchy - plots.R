library(shape)# for good arrow heads
library(geosphere) # for bearing


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



# U: neg log posterior
U <- function(q, X_bearing_diff) {
  rho <- exp(q)
  log_prob <- sum(log(sinh(rho))-log(cosh(rho)-cos(X_bearing_diff)))
  return(-log_prob)
}


# gradient function
grad_U <- function(q, X_bearing_diff) {
  rho <- exp(q)
  n = length(X_bearing_diff)
  d_rho = n*cosh(rho)/sinh(rho) - sum(sinh(rho)/(cosh(rho)-cos(X_bearing_diff)))
  return(-d_rho)
}


# function transfer locations to directions
direction <- function(X){
  angle <- (bearing(X)/180)*pi
  angle[-length(angle)]
}


# main

X <- locationAll[,1:2]

X_bearing <- direction(X)
X_bearing_diff <- X_bearing[-1] - X_bearing[-length(X_bearing)]

init_rho = 0.9 #dunif(1, 0, 2)

n_step = 1000
q_history <- matrix(nrow=n_step, ncol = 1)
q_history[1,] <- log(init_rho)
L <- 35 
epsilon <- 0.001 

# plot
qt <- rep( NA, n_step )
qt[1] <- 1

yr <- 0.5
plot( NULL , xlim=c(1,L*(n_step-1)) , ylim=c(-yr,yr) , xlab="time" , ylab="position" , yaxt="n" )
axis( 2, at=0 , label=0 )
# plot gaussian contour
# compute contour of log-prob
rho_seq <- seq(from=-0.45,to=0.45,length.out=50)
z1d <- matrix(NA,length(rho_seq),length(rho_seq))
for ( i in 1:length(rho_seq) )
  for ( j in 1:length(rho_seq) )
    z1d[i,j] <- U( rho_seq[i] , X_bearing_diff)
cl <- contourLines( rho_seq , 1:length(rho_seq) , z1d , nlevels=18 )
for ( i in 1:length(cl) ) abline( h=cl[[i]]$x[1] , lwd=0.5+2*abs(cl[[i]]$x[1]) )
xt <- 1


for ( i in 2:n_step ) {
  Q <- HMC(U, grad_U, epsilon, L, q_history[i-1,], X_bearing_diff)
  q_history[i,] <- Q$q
  xpts <- (1:(L+1))+xt-1
  #lines( xpts , Q$traj , col="black" , lwd=2 )
  for ( j in 1:length(Q$traj) ) lines( xpts[j:(j+1)] , Q$traj[j:(j+1)] , lwd=2*abs(Q$ptraj[j])+1 )
  #points( xpts[2:L+1] , Q$traj[2:L+1] , pch=16 , col="white" , cex=0.3 )
  points( xpts[L+1] , Q$traj[L+1] , pch=ifelse( Q$accept==1 , 16 , 2 ) , cex=1.2 )
  xt <- xpts[L+1]
  qt[i] <- xt
  
}
points(qt , q_history , col="white" , pch=16 , cex=0.7 )



samples = exp(q_history[(0.1*n_step):n_step, ])
effective_sample = unique(samples)
effective_sample_size = length(effective_sample)
accept_rate = effective_sample_size / length(samples)
accept_rate


# estimated rho
mean(samples)
