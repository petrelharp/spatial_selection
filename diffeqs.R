source("sim-patchy-selection-fns.R")
require(deSolve)
require(Matrix)

###
# Solve various equations,
# i.e. the equilibrium frequency
# and the stationary distr'n of the random walk in it

# Default parameters
params <- list(
        mu = 0,           # mutation rate
        r = 0.3,          # reproduction rate
        m = 0.2,         # probability of migration
        N = 1000,         # number of indivs per pop
        range = c(1,1001), # dimensions of species range
        ntypes = 2,       # number of types 
        patchsize = 10,   # patch radius
        sb = .05,
        sm = -.01
    )

# postcompute
if (min(params$range)>1) { # 2D
    params$migrsteps <- c( lapply( 1:5, function (n) { c( 2^(-(n+2)), n, 0 ) } ),  ## 2D
            lapply( 1:5, function (n) { c( 2^(-(n+2)), -n, 0 ) } ),
            lapply( 1:5, function (n) { c( 2^(-(n+2)), floor(sqrt(n/2)), floor(sqrt(n/2)) ) } ), 
            lapply( 1:5, function (n) { c( 2^(-(n+2)), -floor(sqrt(n/2)), floor(sqrt(n/2)) ) } ), 
            lapply( 1:5, function (n) { c( 2^(-(n+2)), floor(sqrt(n/2)), -floor(sqrt(n/2)) ) } ), 
            lapply( 1:5, function (n) { c( 2^(-(n+2)), -floor(sqrt(n/2)), -floor(sqrt(n/2)) ) } ), 
            lapply( 1:5, function (n) { c( 2^(-(n+2)), 0, n ) } ), 
            lapply( 1:5, function (n) { c( 2^(-(n+2)), 0, -n ) } ) )
} else if (params$range[1]==1) {
    params$migrsteps = c( lapply( 1:5, function (n) { c( 2^(-(n+2)), n, 0 ) } ),  ## 1D, y
                          lapply( 1:5, function (n) { c( 2^(-(n+2)), -n, 0 ) } ) )
} else {
    params$migrsteps = c( lapply( 1:5, function (n) { c( 2^(-(n+2)), 0, n ) } ),  ## 1D, x
                          lapply( 1:5, function (n) { c( 2^(-(n+2)), 0, -n ) } ) )
}
params$migrsteps <- { msum <- sum( sapply(params$migrsteps,"[",1) ); lapply( params$migrsteps, function (x) { c(x[1]/msum,x[-1]) } ) }
# define patch
params$s <- with(params, matrix( sm, nrow=range[1], ncol=range[2] ) )
params$s[ (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 <= (params$patchsize/2)^2 ] <- params$sb  # a central circle
params$sigma <- getsigma(params)
# compute distance to center of patch
params$patchdist <- sqrt( (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 ) 

# and compute migration matrix
params$migrmat <- Matrix(0, nrow=prod(params$range), ncol=prod(params$range) )
locs <- matrix( nrow=params$range[1], ncol=params$range[2] )
xlocs <- col(locs)
ylocs <- row(locs)
for (step in params$migrsteps) {
    to.x <- ( xlocs + step[2] )
    to.y <- ( ylocs + step[3] )
    goodones <- ( to.x > 0 ) & ( to.x <= params$range[2] ) & ( to.y > 0 ) & ( to.y <= params$range[1] )
    sources <- (1:prod(params$range))[goodones]
    targets <- ( (to.x-1)*params$range[1] + to.y )[goodones]
    params$migrmat[ cbind(sources,targets) ] <- params$migrmat[ cbind(sources,targets) ] + step[1]
    params$migrmat[ cbind(sources,sources) ] <- params$migrmat[ cbind(sources,sources) ] - step[1]
}
params$migrmat <- params$migrmat * params$m
diag(params$migrmat) <- 1+diag(params$migrmat)
stopifnot(all.equal(rowSums(params$migrmat),rep(1.0,nrow(params$migrmat))))

# the generator for a Moran process
generator.moran <- function (t,y,parms,...) {
    # y is the current set of proportions
    # this function returns the rate of change of each, i.e.
    #    (d/dt) y[i] = (1-y[i]) * sum_j migr[j -> i] * (1+s[j]) * y[j] - y[i] * sum_j migr[j -> i] * (1-y[j])
    # or i.e.
    #    (d/dt) y = (1-y) * ( M %*% (1+s)*y ) - y * ( M %*% (1-y) )
    #             = ( M %*% (1+s)*y ) - y * ( 1 + M %*% s*y )
    #             = ( M %*% y ) + ( M %*% s*y ) - y * ( 1 + M %*% s*y )
    #             = ( (M-I) %*% y ) + (1-y) * ( M %*% s*y )
    ydim <- dim(y)
    y <- as.vector(y)
    y1 <- ( (1-y) * ( parms$migrmat %*% ( (1+as.vector(parms$s)) * y ) ) - y * ( parms$migrmat %*%  (1-y) ) )
    if (!is.null(ydim)) { dim(y1) <- ydim } else { y1 <- as.vector(y1) }
    return(list(y1))
}

# the generator for a Wright-Fisher process
generator.wf <- function (t,y,parms,...) {
    ydim <- dim(y)
    y <- as.vector(y)
    new.offspring <- crossprod(parms$migrmat, as.vector(parms$s)*y)
    y1 <- new.offspring/sum(new.offspring) - y
    if (!is.null(ydim)) { dim(y1) <- ydim } else { y1 <- as.vector(y1) }
    return(list(y1))
}

yinit <- as.numeric( params$s > 0 )
dyinit <- generator(t=1,y=yinit,parms=params)

yint.moran <- ode( y=yinit, times=30*(1:100), func=generator.moran, parms=params, method='lsoda' )
yint.wf <- ode( y=yinit, times=30*(1:100), func=generator.wf, parms=params, method='lsoda' )

yint <- yint.moran

if (FALSE) {
    layout(1:2)
    matplot(yint[,-1],type='l')
    matplot(t(yint[,-1]),type='l')

    plot( yint.moran[nrow(yint.moran),-1], yint.wf[nrow(yint.wf),-1] )
}

# have converged?
niter <- 0
tol <- 1e-4
while( ( ( sqrt(max(colMeans(apply(yint[,-1],2,diff)^2))) > tol ) | ( max(abs(generator(1,y=yint[nrow(yint),-1],parms=params)[[1]])) > tol ) ) & niter < 100 ) {
    cat('.')
    yint <- ode( y=yint[nrow(yint),-1], times=10*(1:100), func=generator, parms=params, method='lsoda' )
    niter <- niter+1
}
stopifnot( all( yint[,-1] >= 0 ) & all( yint[,-1] <= 1 ) )


# equilibrium frequency
equil.p <- yint[nrow(yint),-1]
dim(equil.p) <- params$range


# reverse-time lineage migration *rate* matrix
lineage.migr <- sweep( t(params$migrmat), 2, as.vector(equil.p * (1+params$s)), "*" ) 
diag(lineage.migr) <- (-1) * ( rowSums(lineage.migr) - diag(lineage.migr) )
lineage.prob <- Diagonal(nrow(lineage.migr)) + lineage.migr / max(abs( lineage.migr ))

# stationary distribution, by iteration
iter.pi <- as.vector(equil.p)/sum(equil.p)
for (k in 1:1000) { iter.pi <- crossprod( lineage.prob, iter.pi ) }
iter.pi1 <- iter.pi
for (k in 1:1000) { iter.pi <- crossprod( lineage.prob, iter.pi ) }
stopifnot( all( abs(iter.pi-iter.pi1) < tol ) )

# stationary distribution, by equation
rw.pi <- ( equil.p * (1+params$s) )
rw.pi <- rw.pi / sum(rw.pi)

# it's narrower than p(x)
plot(ep/sum(ep),col='red',type='l',xlim=c(400,600))
lines(as.vector(rw.pi), col='blue',lwd=2)
lines(iter.pi,col='green',lty=2,lwd=2)

# compute probability that has not coalesced by time t
rho <- 2
recomb.rate <- rho * (1-as.vector(equil.p))

killed.generator <- function (t,y,parms,...) {
    # y is the current prob of not having been killed
    # this function returns the rate of change of this, i.e.
    #   (d/dt) y[x] = sum_y m(x<-z) y[z] - y[z] - recomb[z]
    y1 <- as.vector( lineage.migr %*% y - y - recomb.rate*y )
    return(list(y1))
}

dt <- 1/(rho*10)
prob.survival <- ode( y=rep(1,length.out=length(recomb.rate)), times=dt * (1:3000), func=killed.generator, parms=list(), method='lsoda' )
mean.killing <- colSums(prob.survival[,-1]*diff(c(0,prob.survival[,1])))
stopifnot(all(prob.survival[1000,-1]<tol))

plot(mean.killing,col='green',type='l', xlim=c(450,550), ylim=c(0,max(mean.killing)) )
abline(h=1/rho,lty=3)
