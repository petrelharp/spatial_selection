#!/usr/bin/R

source("sim-patchy-selection-fns.R")
source("lineages.R")

# Parameters
params <- list(
        mu = 0,           # mutation rate
        r = 0.3,          # reproduction rate
        m = 0.2,         # probability of migration
        N = 1e8,         # number of indivs per pop
        range = c(1,1001), # dimensions of species range
        ntypes = 2,       # number of types 
        patchsize = 1000,   # patch radius
        sb = .05,
        sm = -.01
    )


params$migrsteps = c( lapply( 1:5, function (n) { c( 2^(-(n+2)), n, 0 ) } ),  ## 1D, y
                      lapply( 1:5, function (n) { c( 2^(-(n+2)), -n, 0 ) } ) )
params$migrsteps <- { msum <- sum( sapply(params$migrsteps,"[",1) ); lapply( params$migrsteps, function (x) { c(x[1]/msum,x[-1]) } ) }
# define patch
params$s <- with(params, matrix( sm, nrow=range[1], ncol=range[2] ) )
params$s[ (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 <= (params$patchsize/2)^2 ] <- params$sb  # a central circle
params$sigma <- getsigma(params)
# compute distance to center of patch
params$patchdist <- sqrt( (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 ) 

initpop <- newpop(params,ntypes=2,nseeds=0)
startpos <- 500
initpop$n[,startpos,1] <- 1e8-1e3
initpop$n[,startpos,2] <- 1e3


# Generate and save runs
pophists <- lapply(1:10, function (k) pophistory( pop=initpop, nsteps=500, step=1 ) )

# compute mean-squared deviation in each generation
dists <- ( (1:params$range[2]) - startpos )
meandists <- sapply(pophists, function (pophist) apply( pophist$pophist[,,2,], 2, function (x) sum(x*dists)/sum(x) ) )
meansqdists <- sapply(pophists, function(pophist) apply( pophist$pophist[,,2,], 2, function (x) sum(x*dists^2)/sum(x) ) )
numbers <- sapply( pophists, function (pophist) colSums( pophist$pophist[,,2,] ) )

range(meandists)

# check estimated sigma?
layout(1:2)
matplot(meansqdists,type='l',main='forward')
abline(0,getsigma(params)^2,col='red')
matplot(nrow(meansqdists):1,(-1)*sweep(meansqdists,2,meansqdists[nrow(meansqdists),],"-"),type='l',main="reversed")
abline(0,getsigma(params)^2,col='red')

# check interpretation of s: GROWTH RATE IS log((1+params$r*(1+params$sb))/(1+params$r))
matplot(numbers,type='l',log='y')
abline(log10(1e3),log10((1+params$r*(1+params$sb))/(1+params$r)),untf=FALSE)
