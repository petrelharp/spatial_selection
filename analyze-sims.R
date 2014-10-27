#!/usr/bin/Rscript

source("sim-patchy-selection-fns.R")

mutsims <- read.table("mutation-sims-info.tsv",sep='\t',header=TRUE)

names(mutsims)[names(mutsims)=="N"] <- "N.name"
names(mutsims)[names(mutsims)=="N.1"] <- "N"
names(mutsims)[names(mutsims)=="sm"] <- "sm.name"
names(mutsims)[names(mutsims)=="sm.1"] <- "sm"
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

with( subset(mutsims,adapted), plot( muttime, adapttime, log='xy', pch=20 ) )
abline(0,1,untf=TRUE)
abline(coef(muttime.lm),col='blue',lwd=2,untf=TRUE)
for ( smval in 1:nlevels(mutsims$sm.name) ) {
    with( subset(mutsims,adapted & as.numeric(sm.name)==smval), 
        points( 
            sort(unique(muttime)), 
            tapply(adapttime,factor(muttime,levels=sort(unique(muttime))),mean,na.rm=TRUE), 
            cex=2, col=smval )
    )
}
legend("bottomright", legend=levels(mutsims$sm.name), pch=20, col=1:nlevels(mutsims$sm.name) )


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


