#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")


run.names <- c( '3975-1-10000-pophistory-run.Rdata', '5473-1-10000-pophistory-run.Rdata' )

run.id <- gsub( "([^-]*)-.*","\\1",run.names[1] )
load(run.name)
print(c("run",run.id,"range",pophist$pop$params$range))

# time slices of the process
pdf(file='example-equilibrium.pdf', width=6, height=6, pointsize=10)
par(mar=c(5,4,1,1)+.1)
layout(1:2)
timeslice <- floor( seq( .05*dim(pophist$pophist)[4], dim(pophist$pophist)[4], length.out=50 ) )

matplot( pophist$pophist[,,2,timeslice]/pophist$pop$params$N, type='l', xlab='deme number (space)', ylab='allele frequency' )
abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
text( mean(which(as.vector(pophist$pop$params$s)>0)), .100, labels=as.expression(substitute(s==sb,list(sb=max(pophist$pop$params$s)))) )
text( c(.05,.95)*dim(pophist$pophist)[2], .900, labels=as.expression(substitute(s==sb,list(sb=min(pophist$pop$params$s)))) )
lines( rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N, lwd=2 )

matplot( pophist$pophist[,,2,timeslice]/pophist$pop$params$N, type='l', xlab='deme number (space)', ylab='allele frequency', log='y' )
abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
lines( rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N, lwd=2 )

dev.off()

# the exp'l decay
patchloc <- pmax( min(which(pophist$pop$params$s>0))-seq_along(pophist$pop$params$s), seq_along(pophist$pop$params$s)-max(which(pophist$pop$params$s>0)) )


plot( sort(unique(patchloc)), tapply(pophist$occupation[,,2]/(pophist$pop$gen*pophist$pop$params$N),patchloc,mean), log='y' )
points( sort(unique(patchloc)), tapply(rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N,patchloc,mean), pch=20, col='red' )
lines( sort(unique(patchloc)), (1/2) * exp( - sort(unique(patchloc)) * sqrt(2*abs(min(pophist$pop$params$s))) / (pophist$pop$params$sigma) ) )

with(pophist$pop$params, lines( sort(unique(patchloc)), (1/2) * exp( - sort(unique(patchloc)) * sqrt(2*abs(getgrowth(pophist$pop$params)$gm)) / (sigma) ), col='green' ) )


#######

s <- pophist$pop$params$s
migrsteps <- pophist$pop$params$migrsteps
x <- matrix( 0.5, nrow=dim(pophist$pophist)[1], ncol=dim(pophist$pophist)[2] )
f <- function (x) {
    upm <- downm <- (x * 0)
    for (migr in migrsteps) {
        upmigrants <- migr[1] * x * (1+s)
        downmigrants <- migr[1] * (1-x)
        dx <- migr[2]
        if (length(migr)==2) { dy <- 0 } else { dy <- migr[3] }
        nc <- dim(x)[2] # x-extent of populations
        nr <- dim(x)[1] # y-extent of populations
        if (dx>0) {
            upmigrants[,-(1:dx)] <- upmigrants[,-((nc-dx+1):nc),drop=FALSE]
            downmigrants[,-(1:dx)] <- downmigrants[,-((nc-dx+1):nc),drop=FALSE]
        } else if (dx<0) {
            upmigrants[,-((nc+dx+1):nc)] <- upmigrants[,-1:dx,drop=FALSE]
            downmigrants[,-((nc+dx+1):nc)] <- downmigrants[,-1:dx,drop=FALSE]
        }
        if (dy>0) {
            upmigrants[-(1:dy),] <- upmigrants[-((nr-dy+1):nr),,drop=FALSE]
            downmigrants[-(1:dy),] <- downmigrants[-((nr-dy+1):nr),,drop=FALSE]
        } else if (dy<0) {
            upmigrants[-((nc+dy+1):nc),] <- upmigrants[-1:dy,,drop=FALSE]
            downmigrants[-((nc+dy+1):nc),] <- downmigrants[-1:dy,,drop=FALSE]
        }
        upm <- upm + upmigrants
        downm <- downm + upmigrants
    }
    pmigr <- sum( sapply(migrsteps,'[[',1) )
    x <- ( (1-pmigr)*x + upm ) / ( (1-pmigr)*x + upm + downm )
    return(x)
}
xx <- array( 0, dim=c(dim(x),200) )
xx[,,1] <- x
for (k in 2:dim(xx)[3]) {
    tmpx <- xx[,,k-1]
    dim(tmpx) <- dim(xx)[1:2]
    xx[,,k] <- f(tmpx)
}
