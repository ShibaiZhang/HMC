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
      accept = accept ) )
}

# U needs to return neg-log-probability
myU <- function( q , a=0 , b=1 ) {
  U <- sum( dnorm(y,q[1],1,log=TRUE) ) + dnorm(q[1],a,b,log=TRUE)
  return( -U )
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
myU_grad <- function( q , a=0 , b=1 ) {
  G <- sum( y - q[1] ) + (a - q[1])/b^2
  return( -G ) # negative bc energy is neg-log-prob
}

y <- abs(rnorm(50))
y <- c( y , -y )

set.seed(12)
priors <- list()
priors$a <- 0
priors$b <- 1
Q <- list()
ns <- 50
Q$q <- 0
q <- rep( NA , ns )
qt <- rep( NA, ns )
q[1] <- 0
qt[1] <- 1
L <- 20 # 20 for examples
eps <- 0.01 # 0.01 / 0.03
yr <- 0.25
plot( NULL , xlim=c(1,L*(ns-1)) , ylim=c(-yr,yr) , xlab="time" , ylab="position" , yaxt="n" )
axis( 2 , at=0 , label=0 )
# plot gaussian contour
# compute contour of log-prob
mu_seq <- seq(from=-0.45,to=0.45,length.out=50)
z1d <- matrix(NA,length(mu_seq),length(mu_seq))
for ( i in 1:length(mu_seq) )
  for ( j in 1:length(mu_seq) )
    z1d[i,j] <- myU( mu_seq[i] , a=priors$a , b=priors$b )
cl <- contourLines( mu_seq , 1:length(mu_seq) , z1d , nlevels=18 )
for ( i in 1:length(cl) ) abline( h=cl[[i]]$x[1] , lwd=0.5+2*abs(cl[[i]]$x[1]) )
xt <- 1
for ( i in 2:ns ) {
  Q <- HMC2( myU , myU_grad , eps , L , q[i-1] , a=priors$a , b=priors$b )
  if ( Q$accept==1 ) {
    q[i] <- Q$q
  } else {
    q[i] <- q[i-1]
  }
  xpts <- (1:(L+1))+xt-1
  #lines( xpts , Q$traj , col="black" , lwd=2 )
  for ( j in 1:length(Q$traj) ) lines( xpts[j:(j+1)] , Q$traj[j:(j+1)] , lwd=2*abs(Q$ptraj[j])+1 )
  #points( xpts[2:L+1] , Q$traj[2:L+1] , pch=16 , col="white" , cex=0.3 )
  points( xpts[L+1] , Q$traj[L+1] , pch=ifelse( Q$accept==1 , 16 , 2 ) , cex=1.2 )
  xt <- xpts[L+1]
  qt[i] <- xt
}
points( qt , q , col="white" , pch=16 , cex=0.7 )
text( -30 , -0.2 , "south" , xpd=TRUE )
text( -30 , 0.2 , "north" , xpd=TRUE )

q
mean(q)