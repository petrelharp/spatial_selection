source("lineages.R")
source("/home/peter/projects/spatial-sweeps/ParallelMutation/spatial-spread.R")

s <- .1
mu <- 0
r <- 0.1
m <- 0.05
N <- 1000
npops <- 100
ntypes <- 2
migrsteps <- lapply( 1:5, function (n) { c( 2^(-(n+2)), n, 0 ) } )
migrsteps <- c( migrsteps, lapply( 1:5, function (n) { c( 2^(-(n+2)), -n, 0 ) } ) )

# number of lineages
nlin <- 1000
# recomb rate
recomb <- 1e-3
# number of time steps:
nsteps <- 2000

# start near the edge
wavestart <- 1
initpop <- newpop(range=c(1,npops))
initpop$n[1,,2] <- 0  # put an initial population where we want it
initpop$n[1,wavestart,2] <- N
initpop$n[1,wavestart,1] <- N-initpop$n[1,wavestart,2] 

# get good initial values
if (TRUE) {
    pophist <- pophistory( pop=initpop, nsteps=2, step=300 )
    initpop$n <- pophist$pophist[,,,2]
    dim(initpop$n) <- dim(pophist$pophist)[1:3]
}

# do the simulations
pophist <- pophistory( pop=initpop, nsteps=nsteps, step=1)
lins <- lineages( pophist$pophist, nlin=nlin, migr=migrsteps, T=nsteps, m=m )

# save??
save(pophist,initpop,lins,recomb,migrsteps,nlin,nsteps,file="lineage-sims.RData")

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
hetprobs <- sapply(recombvals/recomb, function (relrecomb) { global.het*(1-exp( - relrecomb*lins$recomb[nsteps,] )) } )

# plot time-by-space in one dimension
zimg <- plotPop( aperm(pophist$pophist[1,,,],c(3,1,2)), xlab="time", ylab="space" )
axis(1); axis(2)
# add lineages
# for just lineages do:
plot(0,type="n",xlim=c(0,nsteps),ylim=c(0,npops),xlab="time",ylab="space")
cols <- rainbow(nlin, alpha=min(1,40/nlin))
fullcols <- rainbow(nlin)
for (k in 1:nlin) {
    z <- lins$lineagelocs[,,iord[k]]
    lines( dim(z)[2]:1, z["x",], col=cols[k], lwd=2 )
}

plot(0,ylim=range(lins$localsizes),xlim=c(0,nsteps),type="n", xlab="time", ylab="local population size")
for (k in 1:nlin) { lines(rev(lins$localsizes[1,,k]), col=cols[k]) }

plot(0,ylim=range(lins$recomb),xlim=c(0,nsteps),type="n")
for (k in 1:nlin) { lines( rev(lins$recomb[,k]), col=cols[k]) }

plot(0,ylim=range(lins$coal),xlim=c(0,nsteps),type="n")
for (k in 1:nlin) { lines( rev(lins$coal[,k]), col=fullcols[k]) }

# plot the distrn of the prob of recombination
plot(lins$precomb["Mean",], type="l")
lines(lins$precomb["1st Qu.",], lty=2)
lines(lins$precomb["3rd Qu.",], lty=2)
lines(lins$precomb["Min.",], lty=3)
lines(lins$precomb["Max.",], lty=3)

# plot heterozygosities
# as a function of geographic distance
plot( initdist[iord], hetprobs[iord,1], type="n", ylim=c(0,1), xlab="distance from origin", ylab="heterozygosity" )
legend("topleft", col=c(NA,fullcols), lty=c(NA,1), legend=c("recombination",recombvals) )
for (k in 1:dim(hetprobs)[2]) {
    points(initdist[iord], hetprobs[iord,k], col=cols[k])
    lines( lowess(initdist[iord], hetprobs[iord,k], f=.3), col=fullcols[k] )
}


jpeg(file="het-by-recomb.jpg", width=1000, height=800)
par(cex=1.5)
# ... and as a function of genetic distance
binlocs <- pretty(initdist,n=5)
nbins <- length(binlocs)-1
binmids <- (binlocs[-(nbins+1)]+binlocs[-1])/2
bincols <- rainbow(2*nbins)[1:nbins]
plot( recombvals, hetprobs[1,], type="n", xlim=c(0,max(recombvals)), ylim=c(0,1), xlab="map distance", ylab="heterozygosity" )
legend("topleft", legend=binmids, col=bincols, lty=1, title="Distance to origin", bg="white")
abline(h=global.het, lty=2)
for ( ell in 1:nbins ) {
    iii <- (initdist >= binlocs[ell]) & (initdist < binlocs[ell+1])
    #    lll <- lowess(rep(recombvals,each=sum(iii)), hetprobs[iii,], f=.3)
    #    lines( lll, col=bincols[ell], lwd=2 )
    #    means <- lll$y[match(recombvals,lll$x)]
    rvals <- jitter(recombvals)
    summ <- sapply(1:length(recombvals), function (x) summary(hetprobs[iii,x]) ) 
    lines( recombvals, summ["Mean",], lwd=2, col=bincols[ell] )
    arrows( x0=rvals, y0=summ["1st Qu.",], y1=summ["3rd Qu.",], angle=90, length=0, col=bincols[ell] )
    arrows( x0=rvals, y0=summ["1st Qu.",], y1=summ["3rd Qu.",], angle=90, length=0, col=bincols[ell] )
    for (k in (1:nlin)[iii]) {
        points(jitter(recombvals), hetprobs[k,], col=cols[k], pch=20, cex=0.5)
    }
}
dev.off()
