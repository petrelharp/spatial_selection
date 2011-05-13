
source("Spatial_adaptation/lineages.R")
source("ParallelMutation/spatial-spread.R")

s <- .1
mu <- 0
r <- 0.1
m <- 0.05
N <- 1000
npops <- 1000 
ntypes <- 2
migrsteps <- lapply( 1:10, function (n) { c( 2^(-(n+2)), n, 0 ) } )
migrsteps <- c( migrsteps, lapply( 1:10, function (n) { c( 2^(-(n+2)), -n, 0 ) } ) )

wavestart <- floor(npops/2)
initpop <- newpop(range=c(1,npops))
initpop$n[1,,2] <- 0  # put an initial population where we want it
initpop$n[1,wavestart,2] <- N
initpop$n[1,wavestart,1] <- N-initpop$n[1,wavestart,2] 

pophist <- pophistory( pop=initpop, nsteps=100, step=1 )
lins <- lineages( pophist$pophist, nlin=10, migr=migrsteps, T=100, m=m )


