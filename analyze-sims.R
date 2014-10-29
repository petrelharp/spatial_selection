#!/usr/bin/Rscript

source("sim-patchy-selection-fns.R")

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
        points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
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
