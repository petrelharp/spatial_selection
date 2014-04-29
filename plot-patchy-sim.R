#!/usr/bin/R
source("sim-patchy-selection-fns.R")
require(gsl)

# compute analytic equilibrium? (takes a minute)
do.equilibrium <- FALSE

run.name <- commandArgs(TRUE)[1]

run.id <- gsub( "([^-]*)-.*","\\1",run.name)
load(run.name)
print(c("run",run.id,"range",pophist$pop$params$range))

if (is.null(pophist$occupation) | max(pophist$occupation)==0 ) { stop("occupation densities not recorded") }

dimension <- sum(dim(pophist$pophist)[1:2]>1)
theory.decay <- sqrt(2*abs(getgrowth(pophist$pop$params)$gm)) / getsigma(pophist$pop$params) 
# yes, the parentheses differ in the exponents for x:
expl.approx <- function (x,d) { sqrt(pi/2) * x^((1-d)/2) * exp(-x) }
bessel.approx <- function (x,d) { x^(1-(d/2)) * bessel_Knu( nu=1-d/2, x=x ) } 

# the exp'l decay
patchloc <- with(pophist$pop$params, (-1)^(s>0) * sqrt( abs( ( row(s) - range[1]/2 )^2 + ( col(s) - range[2]/2 )^2 - (patchsize/2)^2 ) ) )
obs.freqs <- pophist$occupation[,,2]/(pophist$pop$params$N*(pophist$pop$gen-pophist$burnin))
expl.freqs <-  expl.approx( (patchloc*theory.decay), d=dimension )
# and the bessel function
bessel.freqs <-  bessel.approx( (patchloc*theory.decay), d=dimension )
expl.const <- exp( mean( log(obs.freqs/expl.freqs)[(log(obs.freqs)<quantile(log(obs.freqs),.5))&(log(obs.freqs)>quantile(log(obs.freqs),.2))], na.rm=TRUE ) )
bessel.const <- exp( mean( log(obs.freqs/bessel.freqs)[(log(obs.freqs)<quantile(log(obs.freqs),.9))&(log(obs.freqs)>quantile(log(obs.freqs),.4))], na.rm=TRUE ) )
tmpdists <- seq(min(patchloc),max(patchloc),length.out=27)
tmplocs <- tmpdists[-1] - diff(tmpdists)/2
patchdist <- cut( as.vector(patchloc), breaks=tmpdists, include.lowest=TRUE, ordered_result=TRUE )

# equilibrium solution to discrete equations
# note this can converge verrry slowly and might not be correct?
if (do.equilibrium) {
    equil <- getequilibrium(pophist$pop$params, keep.convergence=FALSE, nreps=1000, init=obs.freqs)
# equil <- getequilibrium( pophist$pop$params, keep.convergence=TRUE, nreps=100, init=(equil[,,dim(equil)[3]]+equil[dim(equil)[1]:1,dim(equil)[2]:1,dim(equil)[3]])/2 )
}

# distance vs freq
pdf(file=paste(run.id,'expl-decay.pdf',sep='-'), width=6, height=6, pointsize=10)

plot( patchloc, obs.freqs, log='y', xlab='deme number (space)', ylab='allele frequency', pch=20, cex=.5, col=grey(.7) )
abline(v=0,lty=2)
points( tmplocs, tapply(obs.freqs,patchdist,mean), lwd=2 )
lines( tmplocs, bessel.const * expl.approx( tmplocs*theory.decay, d=dimension ), lwd=2, col='green' )
lines( tmplocs, bessel.const * bessel.approx( tmplocs*theory.decay, d=dimension ), lwd=2 )
# points( patchloc, equil, col='red' )
if (do.equilibrium) { lines( tmplocs, tapply(equil,patchdist,mean), col='red' ) }

dev.off()

if (dimension==1) {

    # time slices of the process
    timeslice <- floor( seq( .05*dim(pophist$pophist)[4], dim(pophist$pophist)[4], length.out=50 ) )
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

