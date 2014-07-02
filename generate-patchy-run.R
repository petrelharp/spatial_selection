#!/usr/bin/R
# usage:
#   Rscript generate-patchy-run.R 'nsteps=1e3' 'stepsize=1e3' 'range=c(201,201)'
# or, to restart:
#   Rscript generate-patchy-run.R 'nsteps=1e3' 'stepsize=1e3' 'restart="42217-1-10-pophistory-run.Rdata"' 

# SOURCE THESE BELOW:
# source("sim-patchy-selection-fns.R")
# source("lineages.R")

ngens <- 1e0
nsteps <- 100  # record at this many steps
nsteps <- min(ngens,nsteps)
stepsize <- max(1,floor(ngens/nsteps))
nlins <- 500
do.lineages <- FALSE
burnin <- 5e3
restart <- NULL

# Default parameters
default.params <- list(
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

opts <- commandArgs(TRUE)
print(opts)

# get command line modifications: put them in global here (need restart); in params$ below
for (x in opts) { 
    # cat("eval: "); print(x); 
    eval(parse(text=x)) 
}

# load restarting parameters
if (!is.null(restart)) {
    stopifnot(file.exists(restart))
    renv <- new.env()
    load(restart,envir=renv)
    params <- with(renv,pophist$pop$params)
} else {
    params <- default.params
}

# get command line modifications: put them in params$ here; in global above
# EXCLUDE lines with multiple commands in: won't do the right thing.
for (x in gsub("^([^ <=]*[ <=])","params$\\1",grep(";",opts,value=TRUE,invert=TRUE))) { 
    # cat("eval: "); print(x); 
    eval(parse(text=x)) 
}

if (!exists("run.id")) { run.id <- floor(runif(1)*100000) }

filename <- paste(run.id,"-r",paste(params$range,collapse="-"),"-sb",params$sb,"-sm",params$sm,"-pophistory-run.Rdata",sep="")
run.list <- list.files(".","*-pophistory-run.Rdata")
while (filename %in% run.list) {  # make sure don't overwrite something
    run.id <- floor(runif(1)*10000)
    filename <- paste(run.id,"-r",paste(params$range,collapse="-"),"-sb",params$sb,"-sm",params$sm,"-pophistory-run.Rdata",sep="")
}
if (!interactive()) {
    logfile <- gsub(".Rdata",".Rout",filename,fixed=TRUE)
    logcon <- file(logfile,open="wt")
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output")   # send both to log file
}

# parameters all set up now
print(run.id)
paste( paste( names(params), params, sep="="), collapse=", ")
set.seed(run.id)

source("sim-patchy-selection-fns.R")
source("lineages.R")

# postcompute
if (is.null(params$migrsteps)) {
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
}
# define patch
if (is.null(params$s)) {
    params$s <- with(params, matrix( sm, nrow=range[1], ncol=range[2] ) )
    params$s[ (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 <= (params$patchsize/2)^2 ] <- params$sb  # a central circle
}
# compute distance to center of patch (for plotting)
params$patchdist <- sqrt( (row(params$s)-params$range[1]/2)^2 + (col(params$s)-params$range[2]/2)^2 ) 
# and dispersal distance (for plotting)
params$sigma <- getsigma(params)

if (!is.null(restart)) {
    # already loaded this from previous run 
    initpop <- with(renv,pophist$pop)
    initpop$params <- params
} else {
    # if haven't loaded this from previous run already
    initpop <- newpop(params,ntypes=params$ntypes,nseeds=0)
    if (!exists("initdist")) {
        initpop$n[,,1][params$s>0] <- 0
        initpop$n[,,2][params$s>0] <- params$N
    } else {
        # pass in initial distribution of types
        initpop$n <- initdist
    }
}

# Generate and save runs
pophist <- pophistory( pop=initpop, nsteps=nsteps, step=stepsize, progress=max(100,stepsize), burnin=burnin )

if (do.lineages) {
    ## trace lineages back
    # from outer tail:
    demetotals <- rowSums(pophist$pophist[,,2,200:dim(pophist$pophist)[4],drop=FALSE],dims=2)
    fromdemes <- ( demetotals < sort( demetotals )[ 20 ] + 1 )
    # longlins <- arrayInd( sample( (1:length(pophist$pophist[,,2,]))[fromdemes], nlins, replace=TRUE, prob=pophist$pophist[,,2,][fromdemes] ), .dim=dim(pophist$pophist[,,2,,drop=FALSE]) )
    longlins <- arrayInd( sample( (1:length(pophist$pophist[,,2,]))[fromdemes], nlins, replace=TRUE, prob=(pophist$pophist[,,2,][fromdemes]>0) ), .dim=dim(pophist$pophist[,,2,,drop=FALSE]) )
    colnames(longlins) <- c("y","x","type","t")
    # longlins <- with(pophist, which( (pophist[,,2,,drop=FALSE]>0) & rep(params$patchdist>8*params$sigma,dim(pophist)[4]), arr.ind=TRUE ) )
    # longlins <- longlins[sample(1:nrow(longlins),min(nrow(longlins),500)),]
    lins <- lineages( pophist, T=longlins[,4], linit=longlins[,2:1], migrsteps=params$migrsteps )
} else {
    lins <- NULL
}

save(initpop,pophist,lins,restart,run.id,file=filename)

cat("\n",filename,"\n")
