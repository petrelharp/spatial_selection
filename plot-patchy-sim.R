#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")

run.name <- commandArgs(TRUE)[1]

run.id <- gsub( "([^-]*)-.*","\\1",run.name)
load(run.name)
print(c("run",run.id,"range",pophist$pop$params$range))

dimension <- sum(dim(pophist$pophist)[1:2]>1)
theory.decay <- sqrt(2*abs(getgrowth(pophist$pop$params)$gm) / getsigma(pophist$pop$params) )

# time slices of the process
timeslice <- floor( seq( .05*dim(pophist$pophist)[4], dim(pophist$pophist)[4], length.out=50 ) )

# the exp'l decay
patchloc <- with(pophist$pop$params, (-1)^(s>0) * sqrt( abs( ( row(s) - range[1]/2 )^2 + ( col(s) - range[2]/2 )^2 - (patchsize/2)^2 ) ) )
obs.freqs <- pophist$occupation[,,2]/(pophist$pop$params$N*(pophist$pop$gen-pophist$burnin))
pred.freqs <-  ( if (dimension==2) {patchloc^(-1/2)} else {1} )  * exp( - patchloc * theory.decay )
fit.const <- exp( mean( log(obs.freqs/pred.freqs)[obs.freqs<.01], na.rm=TRUE ) )
tmpdists <- seq(min(patchloc),max(patchloc),length.out=27)
tmplocs <- tmpdists[-1] - diff(tmpdists)/2
patchdist <- cut( as.vector(patchloc), breaks=tmpdists, include.lowest=TRUE, ordered_result=TRUE )
# equil <- getequilibrium(pophist,init=rowMeans(pophist$pophist[,,2,2000:10000,drop=FALSE])/pophist$pop$params$N) #,keep.convergence=TRUE)
if (is.null(pophist$occupation)) { pophist$occupation <- rowSums(pophist$pophist,dim=3); pophist$burnin <- 0 }

if (dimension==1) {
    sliced <- pophist$pophist[,,2,timeslice]/pophist$pop$params$N
    # if (length(dim(sliced))>2) { dim(sliced) <- c( prod(dim(sliced)[1:2]), dim(sliced)[3] ) }

    pdf(file=paste(run.id,'example-equilibrium.pdf',sep='-'), width=6, height=6, pointsize=10)
    par(mar=c(5,4,1,1)+.1)
    layout(1:2)

    matplot( sliced, type='l', xlab='deme number (space)', ylab='allele frequency' )
    abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
    text( mean(which(as.vector(pophist$pop$params$s)>0)), .100, labels=as.expression(substitute(s==sb,list(sb=max(pophist$pop$params$s)))) )
    text( c(.05,.95)*dim(pophist$pophist)[2], .900, labels=as.expression(substitute(s==sb,list(sb=min(pophist$pop$params$s)))) )
    lines( rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N, lwd=2 )

    matplot( sliced, type='l', xlab='deme number (space)', ylab='allele frequency', log='y' )
    abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
    lines( rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N, lwd=2 )

    dev.off()
}

pdf(file=paste(run.id,'expl-decay.pdf',sep='-'), width=6, height=6, pointsize=10)

plot( patchloc, pophist$occupation[,,2]/(pophist$pop$params$N*(pophist$pop$gen-pophist$burnin)), log='y', xlab='deme number (space)', ylab='allele frequency', pch=20, cex=.5, col=grey(.7) )
points( tmplocs, tapply(pophist$occupation[,,2]/(pophist$pop$gen*pophist$pop$params$N),patchdist,mean), lwd=2 )
lines( tmplocs, fit.const * ( if (dimension==2) {tmplocs^(-1/2)} else {1} ) * exp( - tmplocs * theory.decay ), lwd=2 )
abline(v=0,lty=2)

dev.off()
