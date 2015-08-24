#!/usr/bin/Rscript

source("sim-patchy-selection-fns.R")

maxtime <- 25000  # number of gens these were run for

###
# mutation:
#
# All simulations also used a linear grid of 501 demes with a patch of 99 demes in the center, the migration
# model described in Simulation methods, μ = 10 −5 , and s p = .0023 (calculated as the growth rate as
# described in the text). T adapt is the mean time until 100 B alleles were present in the patch, and p adapted is
# the proportion of the simulations that adapted by the 25,000 generations. The simulation began with no B
# alleles.

## made with collate-sims.R
mutsims <- read.table("mutation-sims-info.tsv", sep='\t', header=TRUE )

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
pdf(file="mutation-times-predicted.pdf", width=4, height=5, pointsize=10)
with( subset(mutsims,adapted), 
        plot( muttime, adapttime, log='xy', xlab="mean time to adaptation by mutation", ylab="time to hit 100 in patch",
                xlim=c(120,120000), col=sm.name, pch=as.numeric(N.name) )
        )
abline(0,1,untf=TRUE)
# abline(coef(muttime.lm),col='blue',lwd=2,untf=TRUE)
for ( smval in 1:nlevels(mutsims$sm.name) ) {
    with( subset(mutsims,adapted & as.numeric(sm.name)==smval),  {
        muttime.fac <- factor(format(muttime,digits=1))
        points( 
            tapply( muttime, muttime.fac, mean, na.rm=TRUE ),
            tapply(adapttime, muttime.fac, mean,na.rm=TRUE), 
            cex=3, pch=20 )
    } )
}
legend("bottomright", legend=c(levels(mutsims$sm.name),levels(mutsims$N.name)), pch=c(rep(1,nlevels(mutsims$sm.name)),1:nlevels(mutsims$N.name)), col=c(1:nlevels(mutsims$sm.name),rep(1,nlevels(mutsims$N.name))) )
dev.off()



# look at residuals
if (interactive()) {
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
}


################
# migration
#
# "All simulations also used a linear grid of 501 demes with two patches of 99 demes each, separated by R demes
# in the center; the migration model described in Simulation methods, μ = 10 −5 , and s p = .0023 (calculated
# as the growth rate as described in the text). T adapt is the mean time until 100 B alleles were present in the
# patch, and p adapted is the proportion of the simulations that adapted by the 25,000 generations. At the start
# of the simulation, one patch was initialized with a frequency of 0.8 B alleles, which were absent elsewhere.

## made with collate-sims.R
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
# migsims$migtime <- with( migsims, 0.5/( 4 * N * gb * abs(sm) * exp(- R * sqrt(2*abs(sm))/sigma ) ) + adjust )
# migsims$migtime <- with( migsims, 1/( 5 * ((1+gb)/(1-exp(-(1+gb)))) * N * pmin(2*gb/(1+gb),abs(sm)) * abs(sm) * exp(- R * sqrt(2*abs(sm))/sigma ) ) + adjust )
migsims$migtime <- with( migsims, 1/( 1 * N * pmin(2*gb/(1+gb),abs(sm)) * abs(sm) * rowSums( exp(- outer(R,0:98,"+")*sqrt(2*abs(sm))/sigma ) ) ) + adjust )
migsims$adapted <- ( 
    ( migsims$dims == "1D" ) &
    ( migsims$hit100.2 < 24000 ) &
    ( migsims$time2 >= migsims$hit100.2 ) )

migsims$valid <- with( migsims, 
        ( R*sqrt(abs(sm))/sigma > 1 ) & 
        ( size2 > (sigma/sqrt(abs(sm)))*atan(sqrt(abs(sm/gb))) ) &
        ( abs(sm) * 2 * sigma * N > 1 )
    )


pdf(file="migration-time-predicted.pdf",width=4, height=8, pointsize=10)
layout((1:2))
with( subset(migsims,adapted), {
        plot( migtime, hit100.2, col=sm.name, pch=as.numeric(N.name), log='xy', 
                xlim=c(10,1e8), ylim=c(5,30000),
                ylab="time to hit 100 in second patch", xlab="mean migration time") 
        # points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, col=adjustcolor("red",.8), cex=2 )
        points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), median, na.rm=TRUE ), pch=20, cex=2 )
        abline(0,1,untf=TRUE)
        abline(v=24000, lty=2)
        legend("bottomright",legend=c(levels(sm.name),levels(N.name)), col=c(1:nlevels(sm.name),rep("black",nlevels(N.name))), pch=c(rep(1,nlevels(sm.name)),1:nlevels(N.name)), bg=adjustcolor("white",.75)  )
    } )
with( subset(migsims,adapted), {
        plot( jitter(R), hit100.2, col=sm.name, pch=as.numeric(N.name), log='xy', ylab="time to hit 100 in second patch", xlab="distance between patches", ylim=c(1,30000))
        points( tapply(R,factor(migtime,levels=sort(unique(migtime))),mean),  tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
        legend("bottomright",legend=c(levels(sm.name),levels(N.name)), col=c(1:nlevels(sm.name),rep("black",nlevels(N.name))), pch=c(rep(1,nlevels(sm.name)),1:nlevels(N.name)), bg=adjustcolor("white",.75) )
    } )
dev.off()


if (interactive()) {
    layout(matrix(1:4,nrow=2))
    with( subset(migsims,adapted), {
            plot( migtime, hit100.2, col=R.name, log='xy', xlim=c(10,1e8)) 
            points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
            abline(0,1,untf=TRUE)
            abline(v=24000)
            legend("bottomright",pch=1,legend=levels(R.name), col=1:nlevels(R.name) )
            plot( migtime, hit100.2, col=N.name, log='xy', xlim=c(10,1e8))
            points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
            abline(0,1,untf=TRUE)
            abline(v=24000)
            legend("bottomright",pch=1,legend=levels(N.name), col=1:nlevels(N.name) )
            plot( migtime, hit100.2, col=sm.name, log='xy', xlim=c(10,1e8)) 
            points( sort(unique(migtime)), tapply( hit100.2, factor(migtime,levels=sort(unique(migtime))), mean, na.rm=TRUE ), pch=20, cex=2 )
            abline(0,1,untf=TRUE)
            abline(v=24000)
            legend("bottomright",pch=1,legend=levels(sm.name), col=1:nlevels(sm.name) )
        } )

    pairs( migsims[c("N", "gb", "gm", "R", "sigma", "migtime")],
        upper.panel=function (x,y,...) { points(x,y,...,col=migsims$sm.name) }, 
        lower.panel=function (x,y,...) { points(x,y,...,col=migsims$N.name) } )

    pairs( migsims[c("final1", "time1", "hit100.1", "final2", "time2", "hit100.2","migtime")],
        upper.panel=function (x,y,...) { points(x,y,...,col=migsims$sm.name) }, 
        lower.panel=function (x,y,...) { points(x,y,...,col=migsims$N.name) } )
}

####
# for paper

require(colorspace)
require(calibrate)

pdf( file="times-predicted-observed.pdf", width=6.5, height=3.3, pointsize=10 )
layout(t(1:2))
par(mar=c(4,4,1,1)+.1)
# MUTATION
    # remove simulation configurations where a substantial portion haven't yet adapted
    usethese <- with( mutsims, adapted & ( muttime < maxtime/10 ) ) 
    with( droplevels(subset(mutsims,usethese)), {
            sm.vals <- (levels(droplevels(sm.name))); sm.col.pal <- sequential_hcl(length(sm.vals))
            sm.pch <- match( levels(sm.name)[tapply(sm.name,paramstring,"[",1)], sm.vals )
            sm.cols <- sm.col.pal[ sm.pch ]
            sm.labs <- parse( text=gsub("sm-","s[m]==",levels(sm.name)) )
            N.vals <- (levels(droplevels(N.name))); N.col.pal <- terrain_hcl(length(N.vals))
            N.cols <- N.col.pal[ match( levels(N.name)[tapply(N.name,paramstring,"[",1)], N.vals ) ]
            N.labs <- parse( text=gsub("N-","rho ==",levels(N.name)) )
            mu.vals <- (levels(droplevels(mu.name))); mu.col.pal <- sequential_hcl(length(mu.vals))
            mu.pch <- match( levels(mu.name)[tapply(mu.name,paramstring,"[",1)], mu.vals )
            mu.cols <- mu.col.pal[ mu.pch ]
            mu.labs <- parse( text=gsub("mu","mu == ",levels(mu.name)) )
            xx <- exp(jitter(log(tapply(muttime,paramstring,mean,na.rm=TRUE))))
            plot( 0, type='n', 
                    log='xy', xlim=c(100,10000),ylim=c(100,10000),
                    xlab=expression(1/lambda['mut']),
                    ylab="time to hit 100 in patch",
            ) 
            segments( x0=xx,
                    y0=tapply(adapttime,paramstring,quantile,.25,na.rm=TRUE),
                    y1=tapply(adapttime,paramstring,quantile,.75,na.rm=TRUE),
                    col=mu.cols
                    )
            points( xx,
                    tapply(adapttime,paramstring,median,na.rm=TRUE), 
                    pch=21+(sm.pch%%5), bg=mu.cols, lwd=2, col=N.cols
            ) 
            abline(0,1,untf=TRUE)
            legend("bottomright", cex=0.5, pt.cex=1,
                legend=c(sm.labs,N.labs,mu.labs),
                pch=c( 21+((1:nlevels(sm.name))%%5), rep(1,nlevels(N.name)), rep(21,nlevels(mu.name)) ), 
                col=c("black",N.col.pal)[c(rep(1,nlevels(sm.name)),1+1:nlevels(N.name),rep(1,nlevels(mu.name)))],
                pt.bg=c(NA,sm.col.pal)[c( rep(1,nlevels(sm.name)), rep(1,nlevels(N.name)), 1+(1:nlevels(mu.name)) ) ] 
            )
        } )

# MIGRATION
    usethese <- with( migsims, adapted  & ( migtime < maxtime ) ) 
    with( droplevels(subset(migsims,usethese)), {
            sm.vals <- (levels(droplevels(sm.name))); sm.col.pal <- sequential_hcl(length(sm.vals))
            sm.pch <- match( levels(sm.name)[tapply(sm.name,paramstring,"[",1)], sm.vals ) 
            sm.cols <- sm.col.pal[ sm.pch ]
            sm.labs <- parse( text=gsub("sm-","s[m]==",levels(sm.name)) )
            R.vals <- (levels(droplevels(R.name))); R.col.pal <- terrain_hcl(length(R.vals))
            R.pch <- match( levels(R.name)[tapply(R.name,paramstring,"[",1)], R.vals )
            R.cols <- R.col.pal[ R.pch ]
            R.labs <- parse( text=gsub("R-","R==",levels(R.name)) )
            N.vals <- (levels(droplevels(N.name))); N.col.pal <- sequential_hcl(length(N.vals))
            N.pch <- match( levels(N.name)[tapply(N.name,paramstring,"[",1)], N.vals )
            N.cols <- N.col.pal[ N.pch ]
            N.labs <- parse( text=gsub("N-","rho ==",levels(N.name)) )
            plot( 0, type='n', log='xy', xlim=c(10,3.5e4), ylim=c(10,3.5e4), 
                    ylab="time to hit 100 in second patch", 
                    xlab=expression(1/lambda['mig']) )
            segments( x0=tapply(migtime,paramstring,mean,na.rm=TRUE),
                    y0=tapply(hit100.2,paramstring,quantile,.25,na.rm=TRUE),
                    y1=tapply(hit100.2,paramstring,quantile,.75,na.rm=TRUE),
                    col=N.cols
                    )
            points( tapply(migtime,paramstring,mean,na.rm=TRUE),
                    tapply(hit100.2,paramstring,median,na.rm=TRUE),
                    pch=21+(sm.pch%%5), col=R.cols, bg=N.cols, lwd=2 )
            abline(0,1,untf=TRUE)
            legend("topleft", cex=0.5, pt.cex=1, 
                legend=c(sm.labs,N.labs,R.labs),
                pch=c( 21+((1:nlevels(sm.name))%%5), rep(21,nlevels(N.name)), rep(21,nlevels(R.name))), 
                col=c("black",R.col.pal)[c(rep(1,nlevels(sm.name)),rep(1,nlevels(N.name)), 1+1:nlevels(R.name))],
                pt.bg=c(NA,sm.col.pal)[c( rep(1,nlevels(sm.name)), 1+(1:nlevels(N.name)), rep(1,nlevels(R.name)))] 
            )
    } )
dev.off()

##
# phase space picture of probability of migration versus mutation

pmut <- function(sm,mu,R,A,sigma,gb,C=10/2) { ( A * mu ) / ( A*mu + 2*C*abs(sm)*pmin(1,(1+gb)/(2*gb))*exp(-R*sqrt(2*abs(sm))/sigma^2) ) }
muval <- 1e-5
Aval <- 99
sigmaval <- median(migsims$sigma)  # i.e. most common one, here.
rval <- 0.3
sbval <- 0.01
paramgrid <- expand.grid( 
    N=seq(25,1600,length.out=100), 
    sm=10^(seq(-3,-1,length.out=100)), 
    R=seq(10,160,length.out=100), 
    mu=muval, A=Aval, sigma=sigmaval, r=rval, sb=sbval )
paramgrid$gb <- getgrowth(paramgrid)$gb
paramgrid$pmut <- with( paramgrid, pmut(sm,mu,R,A,sigma,gb) )
paramgrid$valid <- with(paramgrid, 
        ( R*sqrt(abs(sm))/sigma > 1 ) & 
        ( A > (sigma/sqrt(abs(sm)))*atan(sqrt(abs(sm/gb))) ) &
        ( abs(sm) * 2 * sigma * N > 1 )
    )

simgrid <- expand.grid( 
            N=sort(intersect(unique(mutsims$N),unique(migsims$N))), 
            sm=sort(intersect(unique(mutsims$sm),unique(migsims$sm))),
            mu=muval, 
            R=sort(unique(migsims$R)),
            A=Aval,
            sigma=sigmaval,
            r=rval,
            sb=sbval,
            dims="1D",
            stringsAsFactors=FALSE
        )
simgrid$gb <- getgrowth(simgrid)$gb
simgrid$theory.pmut <- with(simgrid, pmut( sm, mu, R, A, sigma, gb ) )
simgrid <- cbind( simgrid, t( sapply( 1:nrow(simgrid), function (k) { params <- simgrid[k,]
            migtimes <- subset( migsims, (N==params$N) & (sm==params$sm) & (R==params$R) & (dims==params$dims) )$hit100.2
            muttimes <- subset( mutsims, (N==params$N) & (sm==params$sm) & (mu==params$mu) & (dims==params$dims) )$hit100.1
            if ( (length(migtimes)>20) & (length(muttimes)>20) ) {
                pmut <- mean(outer(migtimes,muttimes,">"))
            } else {
                pmut <- NA
            }
            c( migtime=median(migtimes), muttime=median(muttimes), pmut=pmut )
        } ) ) )
simgrid$valid <- with(simgrid, 
        ( R*sqrt(abs(sm))/sigma > 1 ) & 
        ( A > (sigma/sqrt(abs(sm)))*atan(sqrt(abs(sm/gb))) ) &
        ( abs(sm) * 2 * sigma * N > 1 )
    )


# look at results
if (interactive()) {
    layout(t(1:2))
    with( simgrid, { 
                plot(pmut, theory.pmut, col=match(sm,sort(unique(sm))) )
                legend("bottomright", legend=sort(unique(sm)), col=1:length(sort(unique(sm))), pch=1 )
                abline(0,1)
            } )
    with( simgrid, { 
                plot(pmut, with(simgrid,pmut(sm,mu,R,A,sigma,gb,C=5)), col=match(sm,sort(unique(sm))), main="C=5" )
                legend("bottomright", legend=sort(unique(sm)), col=1:length(sort(unique(sm))), pch=1 )
                abline(0,1)
            } )
}

# helper functions: do contour plots from expand.grid-type setup
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

# sm against R
pdf(file="prob-mutation-compared.pdf", width=4, height=3, pointsize=10 )
par(mar=c(4,4,3,1)+.1)
usethese <- with(paramgrid, abs(N-100) < mean(diff(sort(unique(N))))/2 )
with( subset(paramgrid,usethese), {
            plot( R, sm, col=grey(pmut), pch=20, log='y', xlim=c(10,180),
                    xlab=expression(R), ylab=expression(s[m]), main="probability of parallel adaptation" )
            con( R, sm, pmut, add=TRUE, col='black', levels=c(.01,.05,.25,.5,.75,.95,.99) )
        } )
with( subset(simgrid, !is.na(pmut) ), {
            Rlocs <- R
            smlocs <- abs(sm) * ifelse(N==1000, ifelse( sm>-2e-2, 1+.3, 1/(1+.3)), 1 )
            points( Rlocs, abs(smlocs), bg=grey(pmut), pch=ifelse(N==1000,21,24), cex=2 ) 
            text( Rlocs, abs(smlocs), formatC(pmut,digits=2,format="f"), 
                    cex=0.5, col='red', pos=4 )
        } )
legend("topright", pch=c(21,24), legend=c(expression(rho==1000),expression(rho==100)), 
        pt.bg=grey(.75), pt.cex=1.5, bg='white' )
dev.off()

####
# table of sims done

require(xtable)

# mutation
usethese <- with( mutsims, dims=="1D" )
params <- c( "mu", "N", "sb", "sm", "hit100.1", "adapted" )
paramtab <- aggregate( subset(mutsims,usethese)[params], subset(mutsims,usethese)['paramstring'], mean )
nsims <- table(mutsims$paramstring)
paramtab$n <- nsims[match(paramtab$paramstring,names(nsims))]
paramtab <- paramtab[ order(paramtab$N,paramtab$sm), -1 ]
names(paramtab)[match(c("N","sm","hit100.1","adapted","mu","sb"),names(paramtab))] <- c("$\\rho$","$s_m$","$T_\\text{adapt}$","$p_\\text{adapted}$","$\\mu$","$s_p$")

sink(file="mutation-parameters-table.tex")
print.xtable(
    xtable( cbind(paramtab[1:floor(nrow(paramtab)/2),],paramtab[(floor(nrow(paramtab)/2)+1):nrow(paramtab),]),
        align=c("r","|",rep("r",ncol(paramtab)),"|","|",rep("r",ncol(paramtab)),"|"),
        label="stab:mutation_params",
        digits=c(0, rep( c(2,0,2,2,0,2,0), 2 )),
        display=c("s",rep( c("g","f","f","f","f","f","f"), 2 )),
        caption="
            Parameter values used in estimates of mean time to adaptation by mutation of figure~\\ref{fig:sim_times}.
            All simulations also used a linear grid of 501 demes with a patch of 99 demes in the center, 
            the migration model described in \\nameref{ss:simulations},
            $\\mu=10^{-5}$, and $s_p=.0023$ (calculated as the growth rate as described in the text).
            $T_\\text{adapt}$ is the mean time until 100 $B$ alleles were present in the patch,
            and $p_\\text{adapted}$ is the proportion of the simulations that adapted by the 25,000 generations.
            The simulation began with no $B$ alleles.
            "
        ), 
    size="\\tiny",
    include.rownames=FALSE, 
    sanitize.colnames.function=identity,
    sanitize.text.function=identity,
    )
sink(NULL)

# migration
usethese <- with( migsims, dims=="1D" )
params <- c( "R", "N", "sb", "sm", "hit100.2", "adapted" )
paramtab <- aggregate( subset(migsims,usethese)[params], subset(migsims,usethese)['paramstring'], mean )
nsims <- table(migsims$paramstring)
paramtab$n <- nsims[match(paramtab$paramstring,names(nsims))]
paramtab <- paramtab[ order(paramtab$N,paramtab$sm), -1 ]
names(paramtab)[match(c("N","sm","sb","hit100.2","adapted"),names(paramtab))] <- c("$\\rho$","$s_m$","$s_p$","$T_\\text{adapt}$","$p_\\text{adapted}$")


sink(file="migration-parameters-table.tex")
print.xtable(
    xtable( cbind(paramtab[1:floor(nrow(paramtab)/2),],paramtab[(floor(nrow(paramtab)/2)+1):nrow(paramtab),]),
        align=c("r","|",rep("r",ncol(paramtab)),"|","|",rep("r",ncol(paramtab)),"|"),
        display=c("s",rep( c("f","f","f","e","f","f","f"), 2 )),
        digits=c(0, rep( c(0,0,2,0,0,2,0), 2 )),
        label="stab:migration_params",
        caption="
            Parameter values used in estimates of mean time to adaptation by migration of figure~\\ref{fig:sim_times}.
            All simulations also used a linear grid of 501 demes with two patches of 99 demes each,
            separated by $R$ demes in the center;
            the migration model described in \\nameref{ss:simulations},
            $\\mu=10^{-5}$, and $s_p=.0023$ (calculated as the growth rate as described in the text).
            $T_\\text{adapt}$ is the mean time until 100 $B$ alleles were present in the patch,
            and $p_\\text{adapted}$ is the proportion of the simulations that adapted by the 25,000 generations.
            At the start of the simulation, one patch was initialized with a frequency of 0.8 $B$ alleles, which were absent elsewhere.
            "
        ), 
    size="\\tiny",
    include.rownames=FALSE, 
    sanitize.colnames.function=identity,
    sanitize.text.function=identity,
    )
sink(NULL)


####
#  make example trace plots

# mutation
if (!file.exists("example-mutation-sims")) {
    dir.create("example-mutation-sims")
    plotsims <- with( subset(mutsims,dims=="1D"), tapply( seq_along(paramstring), droplevels(paramstring), sample, 1 ) )
    for (k in plotsims) {
        with( subset(mutsims,dims=="1D")[k,],  {
            fcon <- pipe(paste("find ../patchy-sim -name", filename))
            fname <- scan(fcon,what='char')
            close(fcon)
            cat("reading from ", fname, "\n")
            load(fname)
            pdf(file=file.path("example-mutation-sims",gsub("_Rdata",".pdf",gsub(".","_",filename,fixed=TRUE))), width=7,height=3,pointsize=10)
            # png(file=gsub(".Rdata",".png",filename),width=7*144,height=3*144,pointsize=10,res=144)
            layout(t(1:2))
            plotpophist(pophist,plotlegend=FALSE,maxtimes=100,
                main=with(pophist$pop$params,parse(text=sprintf("list(rho == %d, mu == %0.2e)",N,mu))) )
            freqs <- pophist$pophist[,,2,]/pophist$pop$params$N
            xlocs <- floor(seq(1,nrow(freqs),length.out=200))
            matplot(xlocs, freqs[xlocs,floor(seq(1,ncol(freqs),length.out=25))],type='l', 
                col=heat.colors(25), ylim=c(0,1),
                ylab="frequency of B", xlab="space (demes)",
                main=with(pophist$pop,parse(text=sprintf("list(s[m] == %0.2e, s[p] == %0.2e)",abs(params$sm),getgrowth(params)$gb))) )
            ss <- as.vector(pophist$pop$params$s)
            lines( xlocs, ((ss-min(ss))/diff(range(ss)))[xlocs], col='black', lty=2, lwd=2 )
            dev.off()
        } )
    }
}

# migration
if (!file.exists("example-migration-sims")) {
    dir.create("example-migration-sims")
    plotsims <- with( subset(migsims,dims=="1D"), tapply( seq_along(paramstring), droplevels(paramstring), sample, 1 ) )
    for (k in plotsims) {
        with( subset(migsims,dims=="1D")[k,],  {
            fcon <- pipe(paste("find ../patchy-sim -name", filename))
            fname <- scan(fcon,what='char')
            close(fcon)
            cat("reading from ", fname, "\n")
            load(fname)
            pdf(file=file.path("example-migration-sims",gsub("_Rdata",".pdf",gsub(".","_",filename,fixed=TRUE))), width=7,height=3,pointsize=10)
            # png(file=gsub(".Rdata",".png",filename),width=7*144,height=3*144,pointsize=10,res=144)
            layout(t(1:2))
            plotpophist(pophist,plotlegend=FALSE,maxtimes=100,
                main=with(pophist$pop$params,parse(text=sprintf("list(rho == %d, R == %0.2d)",N,R))) )
            freqs <- pophist$pophist[,,2,]/pophist$pop$params$N
            xlocs <- floor(seq(1,nrow(freqs),length.out=200))
            matplot(xlocs, freqs[xlocs,floor(seq(1,ncol(freqs),length.out=25))],type='l', 
                col=heat.colors(25), ylim=c(0,1),
                ylab="frequency of B", xlab="space (demes)",
                main=with(pophist$pop,parse(text=sprintf("list(s[m] == %0.2e, s[p] == %0.2e)",abs(params$sm),getgrowth(params)$gb))) )
            ss <- as.vector(pophist$pop$params$s)
            lines( xlocs, ((ss-min(ss))/diff(range(ss)))[xlocs], col='black', lty=2, lwd=2 )
            dev.off()
        } )
    }
}

