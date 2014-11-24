#!/usr/bin/Rscript

source("sim-patchy-selection-fns.R")

maxtime <- 25000  # number of gens these were run for

mutsims <- read.table("mutation-sims-info.tsv",sep='\t',header=TRUE)

names(mutsims)[names(mutsims)=="N"] <- "N.name"
mutsims$N.name <- factor( mutsims$N.name, levels=levels(mutsims$N.name)[order(as.numeric(gsub("N-","",levels(mutsims$N.name))))])
names(mutsims)[names(mutsims)=="N.1"] <- "N"
names(mutsims)[names(mutsims)=="sm"] <- "sm.name"
mutsims$sm.name <- factor( mutsims$sm.name, levels=levels(mutsims$sm.name)[order(as.numeric(gsub("sm-","",levels(mutsims$sm.name))))])
names(mutsims)[names(mutsims)=="sm.1"] <- "sm"
names(mutsims)[names(mutsims)=="mu"] <- "mu.name"
mutsims$mu.name <- factor( mutsims$mu.name, levels=levels(mutsims$mu.name)[order(as.numeric(gsub("mu","",levels(mutsims$mu.name))))])
names(mutsims)[names(mutsims)=="mu.1"] <- "mu"
mutsims$gb <- getgrowth(mutsims)$gb

# create factor for parameter combination
mutsims$paramstring <- with(mutsims, factor( paste(mu.name,sm.name,N.name,dims,sep='_') ) )

mutsims$adjust <- 0  # with( mutsims, log(N)/gb + size1/2 / ( sigma * sqrt(gb) ) )
mutsims$muttime <- with( mutsims, 1/(2 * mu * gb * size1 * N) + adjust )
mutsims$adapted <- ( 
    ( mutsims$dims == "1D" ) &
    ( mutsims$final1 > 0.5 ) &
    ( mutsims$time1 >= mutsims$hit100.1 ) )

mutsims$adapttime <- mutsims$hit100.1

# look at residuals
with( subset(mutsims,adapted), {
        x <- log(N)
        plot( x, adapttime-muttime, log='x' ) 
        points( sort(unique(x)), tapply(adapttime-muttime,factor(x,levels=sort(unique(x))),mean,na.rm=TRUE), col='red', cex=2 )
        resid.lm <- lm( adapttime-muttime ~ x )
        abline(coef(resid.lm),untf=TRUE)
        coef(resid.lm)
    } )

muttime.lm <- with( subset(mutsims,adapted), lm( adapttime ~ muttime ) )
muttime.lm

pdf(file="mutation-times-predicted.pdf", width=5, height=5, pointsize=10)
with( subset(mutsims,adapted), 
        plot( muttime, adapttime, log='xy', xlab="mean time to adaptation by mutation", ylab="time to hit 100 in patch",
                col=sm.name, pch=as.numeric(N.name) )
        )
abline(0,1,untf=TRUE)
# abline(coef(muttime.lm),col='blue',lwd=2,untf=TRUE)
for ( smval in 1:nlevels(mutsims$sm.name) ) {
    with( subset(mutsims,adapted & as.numeric(sm.name)==smval), 
        points( 
            sort(unique(muttime)), 
            tapply(adapttime,factor(muttime,levels=sort(unique(muttime))),mean,na.rm=TRUE), 
            cex=3, pch=20 )
    )
}
legend("bottomright", legend=c(levels(mutsims$sm.name),levels(mutsims$N.name)), pch=c(rep(1,nlevels(mutsims$sm.name)),1:nlevels(mutsims$N.name)), col=c(1:nlevels(mutsims$sm.name),rep(1,nlevels(mutsims$N.name))) )
dev.off()


pairs( mutsims[c("time1", "hit100.1", "final1", "muttime" )  ], 
    upper.panel=function (x,y,...) { points(x,y,...,col=mutsims$sm.name) }, 
    lower.panel=function (x,y,...) { points(x,y,...,col=mutsims$N.name) } )

layout(t(1:2))
with( mutsims, plot( final1, adapttime, col=sm.name, pch=1+!adapted ) )
with( subset(mutsims,adapted&N>100), abline(coef(lm( adapttime ~ final1 ) ) ) )
legend("bottomleft",legend=levels(mutsims$sm.name),col=1:nlevels(mutsims$sm.name),pch=1)
with( mutsims, plot( final1, adapttime, col=N.name, pch=1+!adapted ) )
with( subset(mutsims,adapted&N>100), abline(coef(lm( adapttime ~ final1 ) ) ) )
legend("bottomleft",legend=levels(mutsims$N.name),col=1:nlevels(mutsims$N.name),pch=1)


################
# migration

migsims <- read.table("migration-sims-info.tsv",sep='\t',header=TRUE)

names(migsims)[names(migsims)=="N"] <- "N.name"
migsims$N.name <- factor( migsims$N.name, levels=levels(migsims$N.name)[order(as.numeric(gsub("N-","",levels(migsims$N.name))))])
names(migsims)[names(migsims)=="N.1"] <- "N"
names(migsims)[names(migsims)=="sm"] <- "sm.name"
migsims$sm.name <- factor( migsims$sm.name, levels=levels(migsims$sm.name)[order(as.numeric(gsub("sm-","",levels(migsims$sm.name))))])
names(migsims)[names(migsims)=="sm.1"] <- "sm"
names(migsims)[names(migsims)=="R"] <- "R.name"
migsims$R.name <- factor( migsims$R.name, levels=levels(migsims$R.name)[order(as.numeric(gsub("R-","",levels(migsims$R.name))))])
migsims$R <- as.numeric(gsub("R-","",migsims$R.name))
migsims$gb <- getgrowth(migsims)$gb
migsims$gm <- getgrowth(migsims)$gm

migsims$paramstring <- with(migsims, factor( paste(sm.name,R.name,N.name,dims,sep='_') ) )

migsims$adjust <- 0  # with( migsims, log(N)/gb + size1/2 / ( sigma * sqrt(gb) ) )
migsims$migtime <- with( migsims, 0.5/( 4 * N * gb * abs(sm) * exp(- R * sqrt(2*abs(sm))/sigma ) ) + adjust )
migsims$adapted <- ( 
    ( migsims$dims == "1D" ) &
    ( migsims$hit100.2 < 24000 ) &
    ( migsims$time2 >= migsims$hit100.2 ) )

pdf(file="migration-time-predicted.pdf",width=9, height=5, pointsize=10)
layout(t(1:2))
with( subset(migsims,adapted), {
        plot( migtime, hit100.2, col=sm.name, pch=as.numeric(N.name), log='xy', xlim=c(10,1e8), ylab="time to hit 100 in second patch", xlab="mean migration time") 
        # points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, col=adjustcolor("red",.8), cex=2 )
        points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), median, na.rm=TRUE ), pch=20, cex=2 )
        abline(0,1,untf=TRUE)
        legend("bottomright",legend=c(levels(sm.name),levels(N.name)), col=c(1:nlevels(sm.name),rep("black",nlevels(N.name))), pch=c(rep(1,nlevels(sm.name)),1:nlevels(N.name)) )
    } )
with( subset(migsims,adapted), {
        plot( jitter(R), hit100.2, col=sm.name, pch=as.numeric(N.name), log='xy', ylab="time to hit 100 in second patch", xlab="distance between patches")
        points( tapply(R,factor(migtime,levels=sort(unique(migtime))),mean),  tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
        legend("bottomright",legend=c(levels(sm.name),levels(N.name)), col=c(1:nlevels(sm.name),rep("black",nlevels(N.name))), pch=c(rep(1,nlevels(sm.name)),1:nlevels(N.name)) )
    } )
dev.off()


layout(matrix(1:4,nrow=2))
with( subset(migsims,adapted), {
        plot( migtime, hit100.2, col=R.name, log='xy', xlim=c(10,1e8)) 
        points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
        abline(0,1,untf=TRUE)
        legend("bottomright",pch=1,legend=levels(R.name), col=1:nlevels(R.name) )
        plot( migtime, hit100.2, col=N.name, log='xy', xlim=c(10,1e8))
        points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
        abline(0,1,untf=TRUE)
        legend("bottomright",pch=1,legend=levels(N.name), col=1:nlevels(N.name) )
        plot( migtime, hit100.2, col=sm.name, log='xy', xlim=c(10,1e8)) 
        points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
        abline(0,1,untf=TRUE)
        legend("bottomright",pch=1,legend=levels(sm.name), col=1:nlevels(sm.name) )
    } )

pairs( migsims[c("N", "gb", "gm", "R", "sigma", "migtime")],
    upper.panel=function (x,y,...) { points(x,y,...,col=migsims$sm.name) }, 
    lower.panel=function (x,y,...) { points(x,y,...,col=migsims$N.name) } )

pairs( migsims[c("final1", "time1", "hit100.1", "final2", "time2", "hit100.2","migtime")],
    upper.panel=function (x,y,...) { points(x,y,...,col=migsims$sm.name) }, 
    lower.panel=function (x,y,...) { points(x,y,...,col=migsims$N.name) } )

####
# for paper

require(colorspace)

pdf( file="times-predicted-observed.pdf", width=6.5, height=3.3, pointsize=10 )
layout(t(1:2))
par(mar=c(4,4,1,1)+.1)
# MUTATION
    # remove simulation configurations where a substantial portion haven't yet adapted
    usethese <- with( mutsims, adapted & ( muttime < maxtime/10 ) ) 
    with( droplevels(subset(mutsims,usethese)), {
            sm.vals <- (levels(droplevels(sm.name)))
            sm.col.pal <- sequential_hcl(length(sm.vals))
            sm.cols <- sm.col.pal[ match( levels(sm.name)[tapply(sm.name,paramstring,"[",1)], sm.vals ) ]
            N.vals <- (levels(droplevels(N.name)))
            N.col.pal <- terrain_hcl(length(N.vals))
            N.cols <- N.col.pal[ match( levels(N.name)[tapply(N.name,paramstring,"[",1)], N.vals ) ]
            mu.vals <- (levels(droplevels(mu.name)))
            mu.pch <- match( levels(mu.name)[tapply(mu.name,paramstring,"[",1)], mu.vals )
            xx <- exp(jitter(log(tapply(muttime,paramstring,mean,na.rm=TRUE))))
            plot( 0, type='n', 
                    log='xy', xlim=c(100,10000),ylim=c(100,10000),
                    xlab="mean time to adaptation by mutation", 
                    ylab="time to hit 100 in patch",
            ) 
            segments( x0=xx,
                    y0=tapply(adapttime,paramstring,quantile,.25,na.rm=TRUE),
                    y1=tapply(adapttime,paramstring,quantile,.75,na.rm=TRUE),
                    col=sm.cols
                    )
            points( xx,
                    tapply(adapttime,paramstring,median,na.rm=TRUE), 
                    pch=20+mu.pch, bg=sm.cols, lwd=2, col=N.cols
            ) 
            abline(0,1,untf=TRUE)
            legend("bottomright", cex=0.5, pt.cex=1,
                legend=c(levels(sm.name),levels(N.name),levels(mu.name)), 
                pch=c(rep(21,nlevels(sm.name)),rep(1,nlevels(N.name)),20+1:nlevels(mu.name)), 
                col=c("black",N.col.pal)[c(rep(1,nlevels(sm.name)),1+1:nlevels(N.name),rep(1,nlevels(mu.name)))],
                pt.bg=c(NA,sm.col.pal)[c(1+1:nlevels(sm.name), rep(1,nlevels(N.name)), rep(1,nlevels(mu.name)))] 
            )
        } )
# MIGRATION
    usethese <- with( migsims, adapted  & ( migtime < maxtime ) ) 
    with( droplevels(subset(migsims,usethese)), {
            sm.vals <- (levels(droplevels(sm.name)))
            sm.col.pal <- sequential_hcl(length(sm.vals))
            sm.cols <- sm.col.pal[ match( levels(sm.name)[tapply(sm.name,paramstring,"[",1)], sm.vals ) ]
            R.vals <- (levels(droplevels(R.name)))
            R.col.pal <- terrain_hcl(length(R.vals))
            R.cols <- R.col.pal[ match( levels(R.name)[tapply(R.name,paramstring,"[",1)], R.vals ) ]
            N.vals <- (levels(droplevels(N.name)))
            N.pch <- match( levels(N.name)[tapply(N.name,paramstring,"[",1)], N.vals )
            plot( 0, type='n', log='xy',
                    xlim=c(10,3.5e4), ylim=c(10,3.5e4), 
                    ylab="time to hit 100 in second patch", xlab="mean migration time") 
            segments( x0=tapply(migtime,paramstring,mean,na.rm=TRUE),
                    y0=tapply(hit100.2,paramstring,quantile,.25,na.rm=TRUE),
                    y1=tapply(hit100.2,paramstring,quantile,.75,na.rm=TRUE),
                    col=sm.cols
                    )
            points( tapply(migtime,paramstring,mean,na.rm=TRUE),
                    tapply(hit100.2,paramstring,median,na.rm=TRUE),
                    pch=20+N.pch, col=R.cols, bg=sm.cols, lwd=2 )
            abline(0,1,untf=TRUE)
            legend("topleft", cex=0.5, pt.cex=1,
                legend=c(levels(sm.name),levels(N.name),levels(R.name)), 
                pch=c(rep(21,nlevels(sm.name)),20+(1:nlevels(N.name)),rep(21,nlevels(R.name))), 
                col=c("black",R.col.pal)[c(rep(1,nlevels(sm.name)),rep(1,nlevels(N.name)), 1+1:nlevels(R.name))],
                pt.bg=c(NA,sm.col.pal)[c(1+1:nlevels(sm.name), rep(1,nlevels(N.name)), rep(1,nlevels(R.name)))] 
            )
    } )
dev.off()

##
# phase space picture of probability of migration versus mutation

pmut <- function(sm,mu,R,A,sigma,C=1/2) { ( A * mu ) / ( A*mu + 2*C*sm*exp(-R*sqrt(2*sm)/sigma^2) ) }
paramgrid <- expand.grid( N=seq(25,1600,length.out=100), sm=10^(seq(-4,-1,length.out=100)), mu=1e-6, R=seq(20,160,length.out=100), A=99, sigma=0.953770108230929 )

paramgrid$pmut <- with( paramgrid, pmut(sm,mu,R,A,sigma) )

# helper functions
con <- function (xvar, yvar, zvar, ... ) {
    xvals <- sort(unique(xvar))
    yvals <- sort(unique(yvar))
    zvals <- matrix( NA, nrow=length(xvals), ncol=length(yvals) )
    zvals[ cbind( match(xvar,xvals), match(yvar,yvals) ) ] <- zvar
    contour( x=xvals, y=yvals, z=zvals, ... )
}
fcon <- function (xvar, yvar, zvar, levels=pretty(range(zvar), nlevels), nlevels=20, col=grey(seq(0,1,length.out=nlevels)) ) {
    xvals <- sort(unique(xvar))
    yvals <- sort(unique(yvar))
    zvals <- matrix( NA, nrow=length(xvals), ncol=length(yvals) )
    zvals[ cbind( match(xvar,xvals), match(yvar,yvals) ) ] <- zvar
    .filled.contour( x=xvals, y=yvals, z=zvals, levels=levels, col=col )
}

layout(t(1:2))
# sm against N
usethese <- with(paramgrid, abs(R-40) < mean(diff(sort(unique(R))))/2 )
with( subset(paramgrid,usethese), {
            plot( N, sm, col=grey(pmut), pch=20, log='y' )
            con( N, sm, pmut, add=TRUE, col='red', lwd=2 )
            # fcon( N, sm, pmut, levels=seq(0,1,length.out=20), col=grey(seq(0,1,length.out=20)) )
            # points( sm~N, col=grey(pmut), pch=20, log='y' )
        } )
# sm against R
usethese <- with(paramgrid, abs(N-100) < mean(diff(sort(unique(N))))/2 )
with( subset(paramgrid,usethese), {
            plot( R, sm, col=grey(pmut), pch=20, log='y' )
            con( R, sm, pmut, add=TRUE, col='red', lwd=2 )
        } )

