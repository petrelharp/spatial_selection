#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")
source("sim-patchy-selection-fns.R")

run.list <- list.files(".","*-pophistory-run.Rdata")
if (!file.exists("run-info.csv")) {
    rundims <- sapply( run.list, function (x) { tmp <- load(x); if ('pophist' %in% tmp) { dim(pophist$pophist) } else { NA } } )
    write.table(rundims,file="run-info.csv",sep=',')
} else {
    rundims <- read.csv("run-info.csv")
}

# load random run
rm(pophist); rm(lins)
run.name <- sample(run.list,1)
run.id <- gsub( "([^-]*)-.*","\\1",run.name )
load(run.name)
print(c("run",run.id,"range",pophist$pop$params$range))

plotpophist(pophist)

timeslice <- 1000+1:1000
plotpophist(pophist, timeslice=timeslice, maxtimes=length(timeslice))
for (k in 1:dim(lins$lineagelocs)[3]) {
    ltmp <- t(lins$lineagelocs[,timeslice,k])
    usethese <- ( !is.na(ltmp[,2]) &  pophist$pop$params$patchdist[ltmp] > 15 )
    if (any(usethese)) { lines( seq_along(timeslice)[usethese], subset(ltmp,usethese)[,2], col=rainbow(60)[k] ) }
}
points( lins$T - min(timeslice), lins$linit[,2], pch=20, cex=.25, col=adjustcolor("black",.25) )

pdf(file="verylong.pdf",width=100,height=8,pointsize=10)
timeslice <- 1:dim(lins$lineagelocs)[2]
plotpophist(pophist, timeslice=timeslice, maxtimes=length(timeslice))
for (k in 1:dim(lins$lineagelocs)[3]) {
    ltmp <- t(lins$lineagelocs[,timeslice,k])
    usethese <- ( !is.na(ltmp[,2]) &  pophist$pop$params$patchdist[ltmp] > 15 )
    if (any(usethese)) { lines( seq_along(timeslice)[usethese], subset(ltmp,usethese)[,2], col=rainbow(60)[k] ) }
}
points( lins$T - min(timeslice), lins$linit[,2], pch=20, cex=.25, col=adjustcolor("black",.25) )
dev.off()

# 1-D
animpophist(pophist)

# 2-D
animpophist(pophist,nsteps=3)
with(pophist$pop, contour( x=1:nrow(params$patchdist), y=1:ncol(params$patchdist), z=params$patchdist/params$sigma, add=TRUE ) )
plotlins(lins,add=TRUE,subset=1:100)

# Exponential decay in mean occupation density
plot( pophist$pop$params$patchdist, rowMeans(pophist$pophist[,,2,200:dim(pophist$pophist)[4],drop=FALSE],dims=2)/pophist$pop$params$N, log='y', col=adjustcolor("black",.25), pch=20 )
with(pophist$pop$params, 
    lines( seq(0,max(pophist$pop$params$patchdist),length.out=100), 
        1.5*exp(-(sqrt(2*abs(sm))/getsigma(pophist$pop$params))*seq(0,max(pophist$pop$params$patchdist),length.out=100)) ) )
with(pophist$pop$params, 
    lines( seq(0,max(pophist$pop$params$patchdist),length.out=100), 
        1.5*exp(-2*(sqrt(2*abs(sm))/getsigma(pophist$pop$params))*seq(0,max(pophist$pop$params$patchdist),length.out=100)), lty=2 ) )  # fudge factor of 2 

# Look at occupation times of farout stuff
plot.times <- (4*pophist$pop$params$patchsize):(5*pophist$pop$params$patchsize)
matplot( t(pophist$pophist[,plot.times,2,]), type='l', col=rainbow(length(plot.times)), lty=1, xlab='time', ylab='density' )
legend("topright",plot.times,lty=1,col=rainbow(length(plot.times)),legend=plot.times,title='space')

matplot( pophist$pop$params$patchdist, pophist$pophist[,,2,], type='l', col=adjustcolor(rainbow(dim(pophist$pophist)[4]),.25), lty=1, ylim=c(0,50) )
with(pophist$pop, abline(v=params$patchsize,lwd=2) )

######
# plot lineages, 1D
initdist <- apply(lins$lineagelocs["x",,],2,max,na.rm=TRUE)
lins$lineagelocs <- lins$lineagelocs[,,order(initdist),drop=FALSE]
lins$localsizes <- lins$localsizes[,,order(initdist),drop=FALSE]
initdist <- initdist[order(initdist)]
nlin <- dim(lins$lineagelocs)[2]
plotpop( aperm(pophist$pophist[1,,,],c(3,1,2)), params=pophist$pop$params, xlab="time", ylab="space" )
axis(1); axis(2)
# plot(0,type="n",xlim=c(0,nlin),ylim=c(0,max(dim(lins$lineagelocs)[1:2])),xlab="time",ylab="space")
cols <- rainbow(dim(lins$lineagelocs)[3], 0.1+.7*(initdist/max(initdist)))
matlines( jitter(lins$lineagelocs["x",,]), col=cols, type='l' )


#######
# compute functions of the path
lins$recomb <- apply(lins$localsizes, 3, function (x) { z <- (1-x/pophist$pop$params$N); cumsum( c(z[!is.na(z)],rep(NA,sum(is.na(z)))) ) } )
lins$coal <- apply(lins$localsizes, 3, function (x) { z <- 1/(x*pophist$pop$params$N); cumsum( c(z[!is.na(z)],rep(NA,sum(is.na(z)))) ) } )

