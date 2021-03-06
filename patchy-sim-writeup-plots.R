#!/usr/bin/R
source("sim-patchy-selection-fns.R")
source("lineages.R")
require(gsl)

plotfns <- c(
             pdf=function (...) { pdf(...,pointsize=10) },
             tiff=function (...) { tiff(...,pointsize=10,res=288,units="in") }
         )

# yes, the parentheses differ in the exponents for x:
expl.approx <- function (x,d) { sqrt(pi/2) * x^((1-d)/2) * exp(-x) }
bessel.approx <- function (x,d) { x^(1-(d/2)) * bessel_Knu( nu=1-d/2, x=x ) } 

run.list <- list.files(".","*-pophistory-run.Rdata")
rundims <- read.csv("run-info.csv")

# run.name <- "3975-1-10000-pophistory-run.Rdata"

#######
### occupation frequencies
plotdecay <- function (pophist, shift=FALSE) {
    # return averaged freqs and theoretical prediction
    dimension <- sum(dim(pophist$pophist)[1:2]>1)
    theory.decay <- sqrt(2*abs(getgrowth(pophist$pop$params)$gm)) / getsigma(pophist$pop$params) 
    # the exp'l decay
    patchloc <- with(pophist$pop$params, (-1)^(s>0) * sqrt( abs( ( row(s) - range[1]/2 )^2 + ( col(s) - range[2]/2 )^2 - (patchsize/2)^2 ) ) )
    obs.freqs <- pophist$occupation[,,2]/(pophist$pop$params$N*(pophist$pop$gen-pophist$burnin))
    expl.freqs <-  expl.approx( (patchloc*theory.decay), d=dimension )
    # also estimate the shift?
    if (shift) {
        f <- function (x0) { 
            goodones <- (obs.freqs>3e-3) & (patchloc > -x0)
            bfreqs <- (bessel.approx( (patchloc+x0)*theory.decay, d=dimension ))
            mean( ((log(obs.freqs)-log(bfreqs))[goodones] - mean(log(obs.freqs)[goodones]) + mean(log(bfreqs)[goodones]) )^2 ) }
        bessel.shift <- optimize( f, interval=c(0,50) )
        bessel.shift.const <- exp( mean( log(obs.freqs/bessel.approx( (patchloc+bessel.shift$minimum)*theory.decay, d=dimension ))[(log(obs.freqs)<quantile(log(obs.freqs),.9))&(log(obs.freqs)>quantile(log(obs.freqs),.4))], na.rm=TRUE ) )
    } else {
        bessel.shift.const <- bessel.shift <- NA
    }
    # and the bessel function
    bessel.freqs <-  bessel.approx( (patchloc*theory.decay), d=dimension )
    expl.const <- exp( mean( log(obs.freqs/expl.freqs)[(log(obs.freqs)<quantile(log(obs.freqs),.5))&(log(obs.freqs)>quantile(log(obs.freqs),.2))], na.rm=TRUE ) )
    bessel.const <- exp( mean( log(obs.freqs/bessel.freqs)[(log(obs.freqs)<quantile(log(obs.freqs),.9))&(log(obs.freqs)>quantile(log(obs.freqs),.4))], na.rm=TRUE ) )
    plotdists <- seq(min(patchloc),max(patchloc),length.out=27)
    plotlocs <- plotdists[-1] - diff(plotdists)/2
    patchdist <- cut( as.vector(patchloc), breaks=plotdists, include.lowest=TRUE, ordered_result=TRUE )
    return( list( patchloc=patchloc, obs.freqs=obs.freqs, plotlocs=plotlocs, patchdist=patchdist, bessel.const=bessel.const, bessel.shift.const=bessel.shift.const, bessel.shift=bessel.shift, theory.decay=theory.decay, dimension=dimension) )
}

plotinfo <- lapply( c("3994-1000-1000-pophistory-run.Rdata","7705124-r101-101-sb0.05-sm-0.005-pophistory-run.Rdata"), function (run.name) {
    load(run.name); plotdecay(pophist) } )

for (k in seq_along(plotfns)) {
    outfile <- paste("sim-occupation-freqs.",names(plotfns)[k],sep='')
    plotfns[[k]](file=outfile, width=6,height=3)
    # pdf(file="sim-occupation-freqs.pdf",width=6,height=3,pointsize=10)
    # postscript(file="sim-occupation-freqs.eps",width=6,height=3,pointsize=10,paper="special", horizontal=FALSE)
    layout(t(1:2))
    par(mar=c(5,4,1,1)+.1)
    # averaged
    plot( 1, 1, xlim=range(unlist(lapply(plotinfo,"[[","patchloc"))), ylim=range(unlist(lapply(plotinfo,"[[","obs.freqs"))), 
            log='y', xlab='distance from patch', ylab='allele frequency', type='n' )
    mtext("A",side=3,line=0.1,cex=1.5,adj=-0.3,padj=1)
    abline(v=0,lty=2)
    lapply( plotinfo, function (x) { with(x, {
            lines( patchloc, obs.freqs, pch=20, cex=.5, col=grey(.7) )
            lines( plotlocs, bessel.const * bessel.approx( plotlocs*theory.decay, d=dimension ), col=c('red','green')[dimension] )
            points( plotlocs, tapply(obs.freqs,patchdist,mean), pch=20, cex=.5 )
            if (FALSE) {  lines( plotlocs, bessel.shift.const * bessel.approx( (plotlocs+bessel.shift$minimum)*theory.decay, d=dimension ), lwd=2, lty=2 ) }
    } ) } )
    legend('topright',legend=paste("d=",c(1,2)),lty=1,col=c('red','green'))
    # and, snapshots
    run.name <- "3975-1-10000-pophistory-run.Rdata"
    load(run.name)
    dimension <- sum(dim(pophist$pophist)[1:2]>1)
    theory.decay <- sqrt(2*abs(getgrowth(pophist$pop$params)$gm)) / getsigma(pophist$pop$params)
    timeslice <- floor( seq( .05*dim(pophist$pophist)[4], dim(pophist$pophist)[4], length.out=50 ) )
    sliced <- pophist$pophist[,,2,timeslice]/pophist$pop$params$N
    # plot
    matplot( sliced, type='l', xlab='deme number (space)', ylab='allele frequency', col=adjustcolor(2:10,.5) )
    mtext("B",side=3,line=0.1,cex=1.5,adj=-0.3,padj=1)
    abline(v=0.5 + which(diff(as.vector(pophist$pop$params$s))!=0), lty=2 )
    text( mean(which(as.vector(pophist$pop$params$s)>0)), .05, labels=as.expression(substitute(s[p]==sb,list(sb=max(pophist$pop$params$s)))) )
    text( c(.15,.85)*dim(pophist$pophist)[2], .900, labels=as.expression(substitute(s[m]==sm,list(sm=min(pophist$pop$params$s)))) )
    lines( rowMeans(pophist$pophist[,,2,min(timeslice):max(timeslice)])/pophist$pop$params$N, lwd=2 )
    # done
    dev.off()
}


########
### snapshots
run.name <- "3975-1-10000-pophistory-run.Rdata"
load(run.name)
dimension <- sum(dim(pophist$pophist)[1:2]>1)
theory.decay <- sqrt(2*abs(getgrowth(pophist$pop$params)$gm)) / getsigma(pophist$pop$params)

# lineages

for (k in seq_along(plotfns)) {
    outfile <- paste("sim-snapshots.",names(plotfns)[k],sep='')
    plotfns[[k]](file=outfile,width=5,height=3)
    # pdf(file="sim-snapshots.pdf",width=5,height=3,pointsize=10)
    par(mar=c(5,4,1,1)+.1)
    # countours
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
}


###
# watch a lineage cross over
run.name <- "67890-r1-301-sb0.05-sm-0.006-pophistory-run.Rdata"
load(run.name)
dimension <- sum(dim(pophist$pophist)[1:2]>1)


for (k in seq_along(plotfns)) {
    outfile <- paste("sim-transit.",names(plotfns)[k],sep='')
    plotfns[[k]](file=outfile, width=5,height=3)
    # pdf(file="sim-transit.pdf",width=5,height=3,pointsize=10)
    layout(t(1:2),width=c(3,1))
    par(mar=c(5,4,1,1)+.1)
    # make a copy of the slice we want
    tmphist <- pophist
    tmphist$pophist <- tmphist$pophist[,,,(5.5e4:6.3e4),drop=FALSE]
    ntimes <- dim(tmphist$pophist)[4]
    # plotpophist(tmphist,transposed=TRUE,plotlegend=FALSE)
    timeslice <- seq(1,ntimes,length.out=200)
    subph <- (tmphist$pophist[,,2,timeslice])/tmphist$pop$params$N
    image( x=1:dim(tmphist$pophist)[2], y=seq(min(timeslice),max(timeslice),length.out=200), subph, 
        col=colorRampTrans('red',n=32)[c(1,7:32)], xlab='space (demes)', ylab='time (generations)' )
    image( x=1:dim(tmphist$pophist)[2], y=seq(min(timeslice),max(timeslice),length.out=200), subph>0, col=c(NA,adjustcolor("red",.1)), add=TRUE )
    contour( x=1:dim(tmphist$pophist)[2], y=seq(min(timeslice),max(timeslice),length.out=200), subph, levels=c(.5,.05), add=TRUE, col=grey(.4) )
    # lineages
    linit <- matrix( c(
            1, 1, ntimes,
            # 10, 1, ntimes,
            # 20, 1, ntimes,
            # 30, 1, ntimes,
            # 200, 1, ntimes,
            # 190, 1, ntimes ,
            # 180, 1, ntimes ,
            170, 1, ntimes ),
            ncol=3, byrow=TRUE )
    lins <- lineages( tmphist, linit=linit )
    LL <- plotlins(lins,plotit=FALSE)
    ltimes <- seq(from=1,to=dim(LL)[2],length.out=250)
    lcols <- rep("black",length(linit)) # terrain_hcl(length(linit))
    invisible( lapply( 1:dim(LL)[3], function (k) { lines( LL[2,ltimes,k], ltimes, col=lcols[k]) } ) )
    # simulate shared haplotype length
    localprops <- 1-lins$localsizes[,,,drop=TRUE] / tmphist$pop$params$N
    recombs <- list( matrix( rexp(prod(dim(localprops))), nrow=nrow(localprops) ),
            matrix( rexp(prod(dim(localprops))), nrow=nrow(localprops) ) )
    coaltime <- lins$coalevents$t[2]
    visible <- list( matrix( rbinom(prod(dim(localprops)),size=1, prob=ifelse(seq_along(localprops)>coaltime,localprops[,1],0) ), nrow=nrow(localprops) ),
            matrix( rbinom(prod(dim(localprops)),size=1, prob=ifelse(seq_along(localprops)>coaltime,localprops[,2],0) ), nrow=nrow(localprops) ) )
    for (k in 1:2) {
        recombs[[k]][!visible[[k]]] <- Inf
        recombs[[k]] <- cummin( recombs[[k]] )
    }
    par(mar=c(5,1,1,1)+.1)
    changepoints <- lapply( recombs, function (x) { z <- sort(unique(c(1,outer(which(diff(x)!=0),c(1,0),"+")),length(x))); z[z>0] } )
    x1 <- (-1)*recombs[[1]][changepoints[[1]]]
    x2 <- rev(recombs[[2]][changepoints[[2]]])
    xx <- c( -1, x1, max(x1), min(x2), x2, 1 )
    y1 <- changepoints[[1]]
    y2 <- rev(changepoints[[2]])
    yy <- c(1, y1, length(recombs[[1]]), length(recombs[[1]]), y2 , 1)
    plot( xx, yy, type='l', xlim=c(-1,1)*1e-3, ylim=range(timeslice), yaxt='n', ylab='', xlab="shared haplotype" )
    # lines( (-1)*recombs[[2]], seq_along(recombs[[1]]), type='l', xlim=c(-1,1)*1e-3 )
    polygon( x=xx, y=yy, col='red' )
    dev.off()
}

# same sim, longer plot
# lineages
linit <- matrix( c(
        1, 1, max(timeslice),
        10, 1, max(timeslice),
        20, 1, max(timeslice),
        30, 1, max(timeslice),
        200, 1, max(timeslice),
        190, 1, max(timeslice) ,
        180, 1, max(timeslice) ,
        170, 1, max(timeslice) ),
        ncol=3, byrow=TRUE )
lins <- lineages( pophist, linit=linit )
LL <- plotlins(lins,plotit=FALSE)
ltimes <- seq(from=1,to=dim(LL)[2],length.out=length(timeslice))

# plot
for (k in seq_along(plotfns)) {
    outfile <- paste("sim-stationary.",names(plotfns)[k],sep='')
    plotfns[[k]](file=outfile, width=5,height=8)
    # pdf(file="sim-stationary.pdf",width=5,height=8,pointsize=10)
    par(mar=c(5,4,2,2)+.1)
    timeslice <- seq(1,dim(pophist$pophist)[4],length.out=1000)
    subph <- (pophist$pophist[,,2,timeslice])/pophist$pop$params$N
    image( x=1:dim(pophist$pophist)[2], y=timeslice, subph, 
        col=colorRampTrans('red',n=32)[c(1,7:32)], xlab='space (demes)', ylab='time (generations)' )
    image( x=1:dim(pophist$pophist)[2], y=timeslice, subph>0, col=c(NA,adjustcolor("red",.1)), add=TRUE )
    contour( x=1:dim(pophist$pophist)[2], y=timeslice, subph, levels=c(.5,.05), add=TRUE, col=grey(.4) )
    lcols <- diverge_hcl(length(linit), h = c(255, 330), l = c(40, 40))
    invisible( lapply( 1:dim(LL)[3], function (k) { 
                lines( LL[2,ltimes,k], ltimes, col=lcols[k]) 
            } ) )
    # coalescent events
    with( lins$coalevents, points( x, t, cex=4 ) )
    dev.off()
}


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
