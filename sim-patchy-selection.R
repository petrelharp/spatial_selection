#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")

run.list <- list.files(".","*-pophistory-run.Rdata")

# load random run
rm(pophist); rm(lins)
run.name <- sample(run.list,1)
run.id <- gsub( "([^-]*)-.*","\\1",run.name )
load(run.name)
print(pophist$pop$params$range)

plotpophist(pophist)

animpophist(pophist)

# 2-D
animpophist(pophist,nsteps=3)
with(pophist$pop, contour( x=1:nrow(params$patchdist), y=1:ncol(params$patchdist), z=params$patchdist/params$sigma, add=TRUE ) )
plotlins(lins,add=TRUE,subset=1:100)

# Exponential decay in mean occupation density
plot( pophist$pop$params$patchdist, rowMeans(pophist$pophist[,,2,200:dim(pophist$pophist)[4],drop=FALSE],dims=2)/pophist$pop$params$N, log='y' )
with(pophist$pop$params, lines( seq(0,max(pophist$pop$params$patchdist),length.out=100), 0.5*exp(-(sqrt(2*abs(sm))/getsigma(pophist$pop$params))*seq(0,max(pophist$pop$params$patchdist),length.out=100)) ) )
with(pophist$pop$params, lines( seq(0,max(pophist$pop$params$patchdist),length.out=100), 0.5*exp(-2*(sqrt(2*abs(sm))/getsigma(pophist$pop$params))*seq(0,max(pophist$pop$params$patchdist),length.out=100)), lty=2 ) )  # fudge factor of 2


# Look at occupation times of farout stuff
plot.times <- (2*pophist$pop$params$patchsize):(4*pophist$pop$params$patchsize)
matplot( pophist$pophist[,plot.times,2,], type='l', col=rainbow(length(plot.times)), lty=1 )
legend("topright",plot.times,lty=1,col=rainbow(length(plot.times)),legend=plot.times)

matplot( pophist$pop$params$patchdist, pophist$pophist[,,2,], type='l', col=adjustcolor(rainbow(dim(pophist$pophist)[4]),.25), lty=1, ylim=c(0,50) )
abline(v=params$patchsize,lwd=2)

######
# plot lineages, 1D
initdist <- apply(lins$lineagelocs["x",,],2,max,na.rm=TRUE)
lins$lineagelocs <- lins$lineagelocs[,,order(initdist),drop=FALSE]
lins$localsizes <- lins$localsizes[,,order(initdist),drop=FALSE]
initdist <- initdist[order(initdist)]
nlin <- dim(lins$lineagelocs)[2]
plotpop( aperm(pophist$pophist[1,,,],c(3,1,2)), xlab="time", ylab="space" )
axis(1); axis(2)
# plot(0,type="n",xlim=c(0,nlin),ylim=c(0,max(dim(lins$lineagelocs)[1:2])),xlab="time",ylab="space")
cols <- rainbow(dim(lins$lineagelocs)[3], 0.1+.7*(initdist/max(initdist)))
matplot( jitter(lins$lineagelocs["x",,]), col=cols, type='l' )


# compute functions of the path
lins$recomb <- apply(lins$localsizes, 3, function (x) { z <- (1-x/pophist$pop$params$N); cumsum( c(z[!is.na(z)],rep(NA,sum(is.na(z)))) ) } )
lins$coal <- apply(lins$localsizes, 3, function (x) { z <- 1/(x*pophist$pop$params$N); cumsum( c(z[!is.na(z)],rep(NA,sum(is.na(z)))) ) } )

