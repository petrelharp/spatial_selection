source("lineages.R")
source("/home/peter/projects/spatial-sweeps/ParallelMutation/spatial-spread.R")

s <- .1
mu <- 0
r <- 0.1
m <- 0.05
N <- 1000
npops <- 100
ntypes <- 2
migrsteps <- lapply( 1:10, function (n) { c( 2^(-(n+2)), n, 0 ) } )
migrsteps <- c( migrsteps, lapply( 1:10, function (n) { c( 2^(-(n+2)), -n, 0 ) } ) )

# number of lineages
nlin <- 1000
# recomb rate
recomb <- 1e-3
# number of time steps:
nsteps <- 1000

# start near the edge
wavestart <- 1
initpop <- newpop(range=c(1,npops))
initpop$n[1,,2] <- 0  # put an initial population where we want it
initpop$n[1,wavestart,2] <- N
initpop$n[1,wavestart,1] <- N-initpop$n[1,wavestart,2] 

# after has run once already:
if (FALSE) {
    initpop$n <- pophist$pophist[,,,100]
    dim(initpop$n) <- dim(pophist$pophist)[1:3]
}

# do the simulations
pophist <- pophistory( pop=initpop, nsteps=nsteps, step=1)
lins <- lineages( pophist$pophist, nlin=nlin, migr=migrsteps, T=nsteps, m=m )

# plot time-by-space in one dimension
zimg <- plotPop( aperm(pophist$pophist[1,,,],c(3,1,2)), xlab="time", ylab="space" )
axis(1); axis(2)
# add lineages
cols <- rainbow(nlin)
for (k in 1:nlin) {
    z <- lins$lineagelocs[,,k]
    lines( dim(z)[2]:1, z["x",], col=cols[k], lwd=2 )
}

# compute functions of the path
lins$recomb <- apply(lins$localsizes, 3, function (x) { cumsum( recomb*(1-x/N) ) })
lins$coal <- apply(lins$localsizes, 3, function (x) { cumsum( 1/(N*x) ) })
# and the probability of recombination
lins$precomb <- apply(1-exp(-lins$recomb), 1, summary)

# initial distance to origin, by lineage
initdist <- apply(lins$lineage[,1,]-c(1,1),2,function(x){sqrt(sum(x^2))})
iord <- order(initdist)  # use this to order by initial distance

#  heterozygosity = P{ two chromosomes differ }
#   = P{ at least one recomb } * global het
global.het <- 0.82
recombvals <- (1:20)/5 * recomb
hetprobs <- sapply(recombvals/recomb, function (relrecomb) { 1-exp( - relrecomb*lins$recomb ) } )
dim(hetprobs) <- c(dim(lins$recomb),length(recombvals))

plot(0,ylim=range(lins$localsizes),xlim=c(0,nsteps),type="n", xlab="time", ylab="local population size")
for (k in 1:nlin) { lines(rev(lins$localsizes[1,,k]), col=cols[k]) }

plot(0,ylim=range(lins$recomb),xlim=c(0,nsteps),type="n")
for (k in 1:nlin) { lines( rev(lins$recomb[,k]), col=cols[k]) }

plot(0,ylim=range(lins$coal),xlim=c(0,nsteps),type="n")
for (k in 1:nlin) { lines( rev(lins$coal[,k]), col=cols[k]) }

# plot the distrn of the prob of recombination
plot(lins$precomb["Mean",], type="l")
lines(lins$precomb["1st Qu.",], lty=2)
lines(lins$precomb["3rd Qu.",], lty=2)
lines(lins$precomb["Min.",], lty=3)
lines(lins$precomb["Max.",], lty=3)
