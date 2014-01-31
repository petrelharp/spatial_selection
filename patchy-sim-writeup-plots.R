#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")

run.list <- list.files(".","*-pophistory-run.Rdata")
rundims <- read.csv("run-info.csv")

# run.name <- "3975-1-10000-pophistory-run.Rdata"

for (run.name in run.list) {

    run.id <- gsub( "([^-]*)-.*","\\1",run.name)
    load(run.name)
    print(c("run",run.id,"range",pophist$pop$params$range))

    dimension <- sum(dim(pophist$pophist)[1:2]>1)
    theory.decay <- dimension * sqrt(2*abs(getgrowth(pophist$pop$params)$gm) / getsigma(pophist$pop$params) )

    # time slices of the process
    timeslice <- floor( seq( .05*dim(pophist$pophist)[4], dim(pophist$pophist)[4], length.out=50 ) )

    # the exp'l decay
    patchloc <- with(pophist$pop$params, (-1)^(s>0) * sqrt( abs( ( row(s) - range[1]/2 )^2 + ( col(s) - range[2]/2 )^2 - (patchsize/2)^2 ) ) )
    tmpdists <- seq(min(patchloc),max(patchloc),length.out=27)
    tmplocs <- tmpdists[-1] - diff(tmpdists)/2
    patchdist <- cut( as.vector(patchloc), breaks=tmpdists, include.lowest=TRUE, ordered_result=TRUE )
    # equil <- getequilibrium(pophist,init=rowMeans(pophist$pophist[,,2,2000:10000,drop=FALSE])/pophist$pop$params$N) #,keep.convergence=TRUE)
    if (is.null(pophist$occupation)) { pophist$occupation <- rowSums(pophist$pophist,dim=3) }

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

    plot( patchloc, rowMeans(pophist$pophist[,,2,2000:10000,drop=FALSE],dim=2)/pophist$pop$params$N, log='y', xlab='deme number (space)', ylab='allele frequency', pch=20, cex=.5, col=grey(.7) )
    points( tmplocs, tapply(pophist$occupation[,,2]/(pophist$pop$gen*pophist$pop$params$N),patchdist,mean), lwd=2 )
    lines( tmplocs, (1/2) * exp( - tmplocs * theory.decay ), lwd=2 )
    abline(v=0,lty=2)

    dev.off()
}

### snapshots
run.name <- "3975-1-10000-pophistory-run.Rdata"
load(run.name)
dimension <- sum(dim(pophist$pophist)[1:2]>1)
theory.decay <- dimension * sqrt(2*abs(getgrowth(pophist$pop$params)$gm) / getsigma(pophist$pop$params) )

# lineages

pdf(file="sim-snapshots.pdf",width=6,height=3,pointsize=10)
layout(t(1:2))
par(mar=c(5,4,1,1)+.1)

timeslice <- floor( seq( .05*dim(pophist$pophist)[4], dim(pophist$pophist)[4], length.out=50 ) )
sliced <- pophist$pophist[,,2,timeslice]/pophist$pop$params$N

matplot( sliced, type='l', xlab='deme number (space)', ylab='allele frequency' )
abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
text( mean(which(as.vector(pophist$pop$params$s)>0)), .05, labels=as.expression(substitute(s==sb,list(sb=max(pophist$pop$params$s)))) )
text( c(.15,.85)*dim(pophist$pophist)[2], .900, labels=as.expression(substitute(s==sb,list(sb=min(pophist$pop$params$s)))) )
lines( rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N, lwd=2 )

maxouttime <- 8250
maxoutloc <- 130
outind <- with(pophist, which(
        ( abs( row( pophist[,,2,] ) - maxoutloc ) < 15 ) &
        ( abs( col( pophist[,,2,] ) - maxouttime ) < 150 ) &
        ( pophist[,,2,] > 0 ), arr.ind=TRUE ) )
sum( pophist$pophist[,,2,][ outind ] )
outind <- cbind( x=outind[,1], y=1, t=outind[,2] )
stopifnot( all(pophist$pophist[cbind(outind[,2:1],2,outind[,3])] > 0 ) )
lins <- lineages( pophist, linit=outind[sample(1:nrow(outind),50),] )

# just lineages
timeslice <- seq(7000,8600,length.out=200)
subph <- (pophist$pophist[,,2,timeslice])/pophist$pop$params$N
cols <- c(list(c("#FFFFFF")),list(colorRampTrans('red',n=32)[c(1,5:32)]))
LL <- plotlins(lins,plotit=FALSE)
Ltimes <- as.numeric(dimnames(LL)[[2]])
useLtimes <- (Ltimes>min(timeslice)) & (Ltimes<max(timeslice))
LL <- LL[,useLtimes,] 
# plotpophist(pophist,timeslice=timeslice,maxtimes=300,plotlegend=FALSE, cols=cols,transposed=TRUE)
image( x=1:dim(pophist$pophist)[2], y=seq(min(timeslice),max(timeslice),length.out=200), subph, col=colorRampTrans('red',n=32)[c(1,7:32)], xlab='space (demes)', ylab='time (generations)' )
image( x=1:dim(pophist$pophist)[2], y=seq(min(timeslice),max(timeslice),length.out=200), subph>0, col=c(NA,adjustcolor("red",.1)), add=TRUE )
contour( x=1:dim(pophist$pophist)[2], y=seq(min(timeslice),max(timeslice),length.out=200), subph, levels=c(.5,.05), add=TRUE, col=grey(.4) )
invisible( lapply( 1:dim(LL)[3], function (k) { lines( LL[2,,k], Ltimes[useLtimes] ) } ) )

dev.off()

### occupation frequencies
pdf(file="sim-occupation-freqs.pdf",width=6,height=3,pointsize=10)
layout(t(1:2))
par(mar=c(5,4,3,1)+.1)
for (run.name in c("5038-1-10000-pophistory-run.Rdata","352-1-10000-pophistory-run.Rdata")) {
    load(run.name)
    dimension <- sum(dim(pophist$pophist)[1:2]>1)
    theory.decay <- dimension * sqrt(2*abs(getgrowth(pophist$pop$params)$gm) / getsigma(pophist$pop$params) )
    # the exp'l decay
    patchloc <- with(pophist$pop$params, (-1)^(s>0) * sqrt( abs( ( row(s) - range[1]/2 )^2 + ( col(s) - range[2]/2 )^2 - (patchsize/2)^2 ) ) )
    tmpdists <- seq(min(patchloc),max(patchloc),length.out=27)
    tmplocs <- tmpdists[-1] - diff(tmpdists)/2
    patchdist <- cut( as.vector(patchloc), breaks=tmpdists, include.lowest=TRUE, ordered_result=TRUE )
    if (is.null(pophist$occupation)) { pophist$occupation <- rowSums(pophist$pophist,dim=3) }
    plot( patchloc, rowMeans(pophist$pophist[,,2,2000:10000,drop=FALSE],dim=2)/pophist$pop$params$N, log='y', xlab='deme number (space)', ylab='allele frequency', pch=20, cex=.5, col=grey(.7), main=substitute(d==thisdim,list(thisdim=dimension)), ylim=c(1e-4,1) )
    points( tmplocs, tapply(pophist$occupation[,,2]/(pophist$pop$gen*pophist$pop$params$N),patchdist,mean), lwd=2 )
    lines( tmplocs, (1/2) * exp( - tmplocs * theory.decay ), lwd=2 )
    abline(v=0,lty=2)
}
dev.off()

if (FALSE) {
    # trippy contours for 2D sims:
    f <- function (k,...) contour(pophist$pophist[,,2,k]/pophist$pop$params$N,levels=c(2/pophist$pop$params$N,.125,.5,.75),...)
    tbase <- 3200; f(tbase); for (k in 1:50) { f(tbase+k,add=TRUE,col=rainbow(60)[k]) }

    pdf(file="contours.pdf",width=10,height=10,pointsize=10)
    layout(matrix(1:4,nrow=2))
    par(mar=c(2,2,0,0)+.1)
    tbase <- 3000; f(tbase); for (k in 1:50) { f(tbase+k,add=TRUE,col=rainbow(60)[k]) }
    text(.9,.95,labels=paste("t =",tbase))
    tbase <- 3050; f(tbase); for (k in 1:50) { f(tbase+k,add=TRUE,col=rainbow(60)[k]) }
    text(.9,.95,labels=paste("t =",tbase))
    tbase <- 3100; f(tbase); for (k in 1:50) { f(tbase+k,add=TRUE,col=rainbow(60)[k]) }
    text(.9,.95,labels=paste("t =",tbase))
    tbase <- 3150; f(tbase); for (k in 1:50) { f(tbase+k,add=TRUE,col=rainbow(60)[k]) }
    text(.9,.95,labels=paste("t =",tbase))
    dev.off()


    ########
    ## more exploration with picking out lineages
    run.name <- "3975-1-10000-pophistory-run.Rdata"
    load(run.name)
    # get lineages
    lins <- lineages( pophist, nlin=20 )
    plotpophist(pophist)
    plotlins(lins)

    # identify a good outlying family
    occupation <- rowSums(pophist$pophist[,,2,])
    outloc <- which( ( occupation < min( occupation[occupation>50] ) ) & ( occupation > 0 ) )
    outtimes <- sapply( outloc, function (x) { which( pophist$pophist[,x,2,] > 0 ) } )
    maxout <- which.max( sapply(outtimes,length) )
    maxoutloc <- outloc[maxout]
    maxouttime <- outtimes[[maxout]][1]

    maxouttime <- 8250
    maxoutloc <- 130
    outind <- with(pophist, which(
            ( abs( row( pophist[,,2,] ) - maxoutloc ) < 15 ) &
            ( abs( col( pophist[,,2,] ) - maxouttime ) < 150 ) &
            ( pophist[,,2,] > 0 ), arr.ind=TRUE ) )
    sum( pophist$pophist[,,2,][ outind ] )
    outind <- cbind( x=outind[,1], y=1, t=outind[,2] )
    stopifnot( all(pophist$pophist[cbind(outind[,2:1],2,outind[,3])] > 0 ) )
    lins <- lineages( pophist, linit=outind[sample(1:nrow(outind),50),] )

    # just lineages
    timeslice <- seq(7000,8600,length.out=200)
    subph <- t(pophist$pophist[,,2,timeslice])/pophist$pop$params$N
    cols <- c(list(c("#FFFFFF")),list(colorRampTrans('red',n=32)[c(1,5:32)]))
    plotpophist(pophist,timeslice=timeslice,maxtimes=300,plotlegend=FALSE,ylim=c(50,150),cols=cols)
    contour( x=seq(min(timeslice),max(timeslice),length.out=200), y=1:dim(pophist$pophist)[2], subph, levels=c(.5,.05), add=TRUE )
    plotlins(lins,lwd=2,col='black',linit=FALSE,coalevents=FALSE)

    # contours etc
    timeslice <- seq(7000,8600,length.out=200)
    subph <- t(pophist$pophist[,,2,timeslice])/pophist$pop$params$N
    plotpophist(pophist,timeslice=timeslice,maxtimes=200,plotlegend=FALSE)
    image( x=seq(min(timeslice),max(timeslice),length.out=200), y=1:dim(pophist$pophist)[2], subph < 1e-3, col=adjustcolor(c(NA,'green'),.5), add=TRUE )
    image( x=seq(min(timeslice),max(timeslice),length.out=200), y=1:dim(pophist$pophist)[2], subph < 1e-2, col=adjustcolor(c(NA,'green'),.5), add=TRUE )
    contour( x=seq(min(timeslice),max(timeslice),length.out=200), y=1:dim(pophist$pophist)[2], subph, levels=c(c(1,10,100)/pophist$pop$params$N,.2,.5), add=TRUE )


}