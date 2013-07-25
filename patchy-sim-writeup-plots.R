#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")
source("sim-patchy-selection-fns.R")

run.list <- list.files(".","*-pophistory-run.Rdata")
rundims <- read.csv("run-info.csv")

run.names <- c( '3975-1-10000-pophistory-run.Rdata', '5473-1-10000-pophistory-run.Rdata' )

run.id <- gsub( "([^-]*)-.*","\\1","3975-1-10000-pophistory-run.Rdata" )
load(run.name)
print(c("run",run.id,"range",pophist$pop$params$range))

# time slices of the process
pdf(file='example-equilibrium.pdf', width=6, height=6, pointsize=10)
par(mar=c(5,4,1,1)+.1)
layout(1:2)
timeslice <- floor( seq( 2000, 10000, length.out=50 ) )

matplot( pophist$pophist[,,2,timeslice]/pophist$pop$params$N, type='l', xlab='deme number (space)', ylab='allele frequency' )
abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
text( mean(which(as.vector(pophist$pop$params$s)>0)), .100, labels=as.expression(substitute(s==sb,list(sb=max(pophist$pop$params$s)))) )
text( c(.05,.95)*dim(pophist$pophist)[2], .900, labels=as.expression(substitute(s==sb,list(sb=min(pophist$pop$params$s)))) )
lines( rowMeans(pophist$pophist[,,2,2000:10000])/pophist$pop$params$N, lwd=2 )

matplot( pophist$pophist[,,2,timeslice]/pophist$pop$params$N, type='l', xlab='deme number (space)', ylab='allele frequency', log='y' )
abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
lines( rowMeans(pophist$pophist[,,2,2000:10000])/pophist$pop$params$N, lwd=2 )

dev.off()

# the exp'l decay
patchloc <- pmax( min(which(pophist$pop$params$s>0))-seq_along(pophist$pop$params$s), seq_along(pophist$pop$params$s)-max(which(pophist$pop$params$s>0)) )

plot( sort(unique(patchloc)), tapply(rowMeans(pophist$pophist[,,2,2000:10000])/pophist$pop$params$N,patchloc,mean), log='y' )
lines( sort(unique(patchloc)), (1/2) * exp( - sort(unique(patchloc)) * sqrt(2*abs(min(pophist$pop$params$s))) / (pophist$pop$params$sigma) ) )


#######
# get equlibrium solution by iterating f() until convergence
# with boundary conditions fixed to zero

xx <- getequilibrium(pophist,keep.convergence=TRUE)
matplot(xx[1,,floor(seq(1,nreps,length.out=100))],type='l',col=rainbow(1.2*100))
