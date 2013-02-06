#!/usr/bin/R

source("sim-patchy-selection-fns.R")
source("lineages.R")

# Parameters
params <- list(
        mu = 0,           # mutation rate
        r = 0.1,          # reproduction rate
        m = 0.05,         # probability of migration
        N = 1000,         # number of indivs per pop
        range = c(51,51), # dimensions of species range
        ntypes = 2,       # number of types 
        patchsize = 10,   # patch radius
        sb = .1,
        sm = -.02
    )
# get command line modifications
for (x in commandArgs(TRUE)) { eval(parse(text=x)) }
for (x in gsub("^([^ <=]*[ <=])","params$\\1",commandArgs(TRUE))) { eval(parse(text=x)) }
# print(params)
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
    params$migrsteps = c( lapply( 1:5, function (n) { c( 2^(-(n+2)), n, 0 ) } ),  ## 1D
            lapply( 1:5, function (n) { c( 2^(-(n+2)), -n, 0 ) } ) )
} else {
    params$migrsteps = c( lapply( 1:5, function (n) { c( 2^(-(n+2)), 0, n ) } ),  ## 1D
            lapply( 1:5, function (n) { c( 2^(-(n+2)), 0, -n ) } ) )
}
params$migrsteps <- { msum <- sum( sapply(params$migrsteps,"[",1) ); lapply( params$migrsteps, function (x) { c(x[1]/msum,x[-1]) } ) }
# define patch
params$s <- with(params, matrix( sm, nrow=range[1], ncol=range[2] ) )
params$s[ (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 <= (params$patchsize/2)^2 ] <- params$sb  # a central circle
params$sigma <- getsigma(params)
# compute distance to patch
params$patchdist <- sqrt( (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 ) 

initpop <- newpop(params,ntypes=2,nseeds=0)
initpop$n[,,1][params$s>0] <- 0
initpop$n[,,2][params$s>0] <- params$N

run.list <- list.files(".","*-pophistory-run.Rdata")

# Generate and save runs
nsteps <- 10000
stepsize <- 1
run.id <- floor(runif(1)*10000)
filename <- paste(run.id,stepsize,nsteps,"pophistory-run.Rdata",sep="-")
while (filename %in% run.list) {  # make sure don't overwrite something
    run.id <- floor(runif(1)*10000)
    filename <- paste(run.id,stepsize,nsteps,"pophistory-run.Rdata",sep="-")
}
set.seed(run.id)
pophist <- pophistory( pop=initpop, nsteps=nsteps, step=stepsize )

## trace lineages back
longlins <- with(pophist, which( (pophist[,,2,,drop=FALSE]>0) & rep(params$patchdist>8*params$sigma,dim(pophist)[4]), arr.ind=TRUE ) )
longlins <- longlins[sample(1:nrow(longlins),min(nrow(longlins),500)),]
lins <- lineages( pophist, T=longlins[,4], linit=longlins[,1:2], migrsteps=params$migrsteps )

save(pophist,lins,file=filename)

