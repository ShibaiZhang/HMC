
HMC2 <- function (U, grad_U, epsilon, L, current_q , ... ) {
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

# U needs to return neg-log-probability
# y ~ normal( mu , exp(log_sigma) )
# mu ~ normal( a, b )
# log_sigma ~ normal( k, d )
myU2 <- function( q , a=0 , b=1 , k=0 , d=0.5 ) {
  s <- exp(q[2]) # sigma on log latent scale
  mu <- q[1]
  U <- sum( dnorm(y,mu,s,log=TRUE) ) + dnorm(mu,a,b,log=TRUE) + dnorm(q[2],k,d,log=TRUE)
  return( -U )
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
myU_grad2 <- function( q , a=0 , b=1 , k=0 , d=0.5 ) {
  mu <- q[1]
  s <- exp(q[2])
  G1 <- sum( y - mu ) * exp(-2*q[2]) + (a - mu)/b^2 #dU/dmu
  G2 <- sum( (y - mu)^2 ) * exp(-2*q[2]) - length(y) + (k-q[2])/d^2 #dU/ds
  return( c( -G1 , -G2 ) ) # negative bc energy is neg-log-prob
}

# test data
set.seed(7)
y <- abs(rnorm(50))
y <- c( y , -y ) # ensure mean is zero

###########
# example paths
library(shape) # for good arrow heads
#blank2()

# priors
priors <- list()
priors$a <- 0
priors$b <- 1
priors$k <- 0
priors$d <- 0.3

#ss <- ss + 1
set.seed(42) # seed 9 for examples
# init
n_samples <- 4
Q <- list()
Q$q <- c(-0.4,0.2)
xr <- c(-0.6,0.6)
yr <- c(-0.25,0.4)
plot( NULL , ylab="log_sigma" , xlab="mu" , xlim=xr , ylim=yr )
step <- 0.02
L <- 12 # 0.02/12 okay sampling --- 0.02/20 is good for showing u-turns
xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55

points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
for ( i in 1:n_samples ) {
  Q <- HMC2( myU2 , myU_grad2 , step , L , Q$q , a=priors$a , b=priors$b , k=priors$k , d=priors$d )
  if ( n_samples < 10 ) {
    lines( Q$traj , lwd=2 )
    points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )
    Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] , arr.length=0.3 , arr.adj = 0.7 )
    text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=0.8 , pos=xpos[i] , offset=0.4 )
  }
  points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 16 , 1 ) , 
          col=ifelse( abs(Q$dH)>0.1 , "red" , "black" ) )
}
#points( 0 , log(sd(y)) , pch=3 ) # mark MAP point
# compute contour of log-prob
cb <- 0.2
mu_seq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=50) 
logsigma_seq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=50)
z <- matrix(NA,length(mu_seq),length(logsigma_seq))
for ( i in 1:length(mu_seq) )
  for ( j in 1:length(logsigma_seq) )
    z[i,j] <- myU2( c( mu_seq[i] , logsigma_seq[j] ) , a=priors$a , b=priors$b , k=priors$k , d=priors$d )
cl <- contourLines( mu_seq , logsigma_seq , z , nlevels=21 )
for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y ,  lwd=0.5 )
mtext( paste("Leapfrog steps = ",L) )
#print(ss)
############