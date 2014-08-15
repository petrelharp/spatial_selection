source("patchy-selection-fns.R")
source("plotting-fns.R")
require(xtable)

#######
# Phase diagram:
#  time to migration and mutations as a function of migration rate, mutation rate, and/or deleterious selection coef

# Plot level sets of time to mut and mig
#  first against characteristic dist and pop density
#  then agains mutation rate and pop density
muttime <- function (env) {with(env,1/(A*mu*2*sb*rho))}
migtime <- function (env) {with(env,1/(const*(sb*4*sd)*rho*(1/sqrt(cl*sqrt(sd)))*exp(-cl*sqrt(sd))))}
rhovals <- logseq(1,1e5,length.out=100) # per km^2
clvals <- seq(1,120,length.out=100)
muvals <- logseq(1e-9,1e-4,length.out=100) 
sdvals <- logseq(1e-6,1e-2,length.out=100)
A <- 100 # km^2
rho <- 10  # per km^2
mu <- 1e-8 
sb <- .1
sd <- .01
const <- .5
cl <- 20  # R/sigma
clplot <- expand.grid(cl=clvals,rho=rhovals)
clplot$mut <- muttime(clplot)
clplot$mig <- migtime(clplot)
muplot <- expand.grid(mu=muvals,rho=rhovals)
muplot$mut <- muttime(muplot)
muplot$mig <- migtime(muplot)
sdplot <- expand.grid(mu=muvals,sd=sdvals)
sdplot$mut <- muttime(sdplot)
sdplot$mig <- migtime(sdplot)
f <- function (x) { matrix(x,nrow=(sqrt(length(x)))) }
pminmax <- function (x,lower,upper) { pmin(upper,pmax(lower,x)) }

pdf(file="phase-diagram-log.pdf",width=6,height=3,pointsize=10)
layout(t(1:2))
par(mar=c(4,3.5,1,1)+.1)
levelsets <- 10^rev(c(0,1,2,3,4,5,6))  # which lines to draw
labelsets <- levelsets
# plot(0,type='n',xlab='',ylab='',xlim=range(clvals),ylim=range(rhovals),log='y')
# contour( clvals, rhovals, f(clplot$mut), col='red', levels=levelsets, labels=labelsets, method="edge", lwd=2, add=TRUE )
logcontour( clvals, rhovals, f(clplot$mut), col='red', levels=levelsets, labels=labelsets, method="edge", log='y' )
logcontour( clvals, rhovals, f(clplot$mig), col='blue', add=TRUE, labels=labelsets, levels=levelsets, method="edge", log='y' )
mtext(side=1,line=2.5,text=expression(R / sigma))
mtext(side=2,line=2.5,expression(rho))
# Tmut == Tmig line:
abline(v=clvals[which.min(abs(migtime(list(cl=clvals))-muttime(.GlobalEnv)))], lwd=2, lty=2)
for (lev in levelsets) { # note rho is log10-scale.
    rholev <- muttime(.GlobalEnv)*rho/lev # mutation time determines rho
    if (rholev>min(rhovals) & rholev<max(rhovals)) {
        polygon( range(clvals)[c(1,1,2,2)], log10(c(max(rhovals),rholev))[c(1,2,2,1)], border=NA, col=adjustcolor("red",.1) )
    }
    miglevs <- 1/(const*(sb*4*sd)*(1/sqrt(clvals*sqrt(sd)))*exp(-clvals*sqrt(sd))*lev)  # migration:
    polygon( c(clvals,rev(clvals)), log10(c( pminmax(miglevs,min(rhovals),max(rhovals)), rep(max(rhovals),length(clvals)) )), border=NA, col=adjustcolor("blue",.1), )
}
legend("topleft",bg="white",legend=c(expression(T[scriptstyle(mut)]),expression(T[scriptstyle(mig)])),fill=c("red","blue"))

# plot(0,type='n',xlab='',ylab='',xlim=range(muvals),ylim=range(rhovals),log='xy')
logcontour( A*muvals, sdvals, f(sdplot$mut), col='red', labels=labelsets, levels=levelsets, method="edge", lwd=2, log='xy' )
logcontour( A*muvals, sdvals, f(sdplot$mig), col='blue', add=TRUE, labels=labelsets, levels=levelsets, lwd=2, log='xy' )
mtext(side=1,line=2.5,expression(mu * A))
mtext(side=2,line=2.5,expression(s[m]))
# Tmut == Tmig line:
lines( log10( muttime(list(A=1,mu=1))/migtime(list(sd=sdvals)) ), log10(sdvals), lwd=2, lty=2 )
for (lev in levelsets) { # note rho and mu are log10-scale.
    # mutation:
    mulev <- muttime(list(mu=1))/lev
    if ( ( mulev>min(muvals) ) & ( mulev < max(muvals) ) ) {
        polygon( log10(A*c(mulev,max(muvals)))[c(1,1,2,2)], log10(range(sdvals))[c(1,2,2,1)], border=NA, col=adjustcolor("red",.1), )
    }
    # migration:
    if (lev > min(sdplot$mig) & lev < max(sdplot$mig)) {
        sdlev <- uniroot(f=function(sd) migtime(list(sd=sd))-lev, lower=min(sdvals), upper=max(sdvals), tol=1e-8)$root
        if (sdlev>min(sdvals) & sdlev<max(sdvals)) {
            polygon( log10(A*range(muvals))[c(1,1,2,2)], log10(c(min(sdvals),sdlev))[c(1,2,2,1)], border=NA, col=adjustcolor("blue",.1), )
        }
    }
}
dev.off()

# ALTERNATIVELY
if (FALSE) {
    logcontour( muvals, rhovals, f(muplot$mut), col='red', labels=labelsets, levels=levelsets, method="edge", lwd=2, log='xy' )
    logcontour( muvals, rhovals, f(muplot$mig), col='blue', add=TRUE, labels=labelsets, levels=levelsets, lwd=2, log='xy' )
    for (lev in levelsets) { # note rho and mu are log10-scale.
        # mutation:
        if (any( (1/(A*muvals*2*sb*lev) < max(rhovals)) & (1/(A*muvals*2*sb*lev) > min(rhovals)) )) {
            polygon( log10(c(muvals,rev(muvals))), log10(c( pminmax(1/(A*muvals*2*sb*lev),min(rhovals),max(rhovals)), rep(max(rhovals),length(muvals)) )), border=NA, col=adjustcolor("red",.1), )
        }
        # migration:
        rholev <- 1/(const*2*sd*(sb+sd)*exp(-cl)*lev)
        if (rholev>min(rhovals) & rholev<max(rhovals)) {
            polygon( log10(range(muvals))[c(1,1,2,2)], log10(c(max(rhovals),rholev))[c(1,2,2,1)], border=NA, col=adjustcolor("blue",.1), )
        }
    }
    mtext(side=1,line=2.5,expression(mu * A))
    mtext(side=2,line=2.5,expression(rho))
    legend("topright",bg="white",legend=c(expression(T[scriptstyle(mut)]),expression(T[scriptstyle(mig)])),fill=c("red","blue"))
}

pdf(file="phase-diagram.pdf",width=6,height=3,pointsize=10)
layout(t(1:2))
par(mar=c(4,3.5,1,1)+.1)
levelsets <- 10^c(0,2,4,6)  # which lines to draw
labelsets <- levelsets
# plot(0,type='n',xlab='',ylab='',xlim=range(clvals),ylim=range(rhovals),log='y')
# contour( clvals, rhovals, f(clplot$mut), col='red', levels=levelsets, labels=labelsets, method="edge", lwd=2, add=TRUE )
logcontour( clvals, rhovals, f(clplot$mut), col='red', levels=levelsets, labels=labelsets, method="edge", lwd=2 )
logcontour( clvals, rhovals, f(clplot$mig), col='blue', add=TRUE, labels=labelsets, levels=levelsets, method="edge", lwd=2 )
mtext(side=1,line=2.5,text=expression(R * sqrt(s[m]) / sigma))
mtext(side=2,line=2.5,expression(rho))
# plot(0,type='n',xlab='',ylab='',xlim=range(muvals),ylim=range(rhovals),log='xy')
logcontour( muvals, rhovals, f(muplot$mut), col='red', labels=labelsets, levels=levelsets, method="edge", lwd=2 )
logcontour( muvals, rhovals, f(muplot$mig), col='blue', add=TRUE, labels=labelsets, levels=levelsets, lwd=2 )
mtext(side=1,line=2.5,expression(mu))
mtext(side=2,line=2.5,expression(rho))
dev.off()


#####
# Prob of establishment as distance from patch

iterroot <- function (genfn,tol=sqrt(.Machine$double.eps),maxiter=1000,...) {
    p0 <- genfn(1/2,...)
    for (k in 1:maxiter) {
        p <- genfn(p0,...)
        if (max(abs(p0-p))<tol) { break } else { p0 <- p }
    }
    attr(p,"niter") <- k
    attr(p,"converged") <- (k<maxiter)
    return(p)
}

dx <- 15            # distance between demes (km)
rangesize <- 500   # half-length of range (km)
xx <- seq(-rangesize,rangesize,by=dx)  # deme locations (km)
sb <- .05
sd <- .05

get.prob.extinct <- function (clinewidth) {
    alpha <- (sb+sd)/clinewidth      # slope of selection with altitude (units of s per km)
    selfn <- function (x) { pmin( sb, pmax( -sd, alpha*x ) ) }
    selx <- selfn(xx)
    # Poisson
    m <- .5
    eqvals <- c( 1, iterroot( function (p) exp((1+sb)*(p-1)) ) )
    genfn <- function (f,selx,migr=m) { exp((1+selx)*(f+migr*diff(c(eqvals[1],f,eqvals[2]),differences=2)-1)) }
    niter <- 1000
    f <- function (x) { 1/2 + sb*tanh(-x/30) }
    ff <- matrix(NA,ncol=niter,nrow=length(xx))
    ff[,1] <- f(xx)
    for (k in 2:niter) { ff[,k] <- genfn(ff[,k-1],selx) }
    return( list( prob.extinct=ff[,niter], selx=selx, eqvals=eqvals, m=m ) )
}

pdf(file="prob-establishment.pdf",width=6,height=3,pointsize=10)
layout(t(1:2))
par(mar=c(4,4,1,1)+.1,mgp=c(2.5,1,0))
with( get.prob.extinct(50), {
        plot(xx,1-prob.extinct,type='n',lwd=2,ylim=c(0,max(1-prob.extinct)),xlab="distance (km)", ylab="prob of establishment")
        polygon(c(xx,rev(xx)),c(1-prob.extinct,rep(0,length(xx))),col=adjustcolor("blue",.5))
        # abline(h=1-eqvals[2], lty=2, lwd=2, col='red')
        # abline(v=c(xx[min(which(selx>-sd))],xx[max(which(selx<sb))]),lty=3,lwd=2)
        polygon(xx[c(min(which(selx>-sd)),max(which(selx<sb)))][c(1,1,2,2)], c(0,par("usr")[4])[c(1,2,2,1)],density=5,col=grey(.5))
        text( xx[min(which(selx>-sd))], max(1-prob.extinct), labels=as.expression(substitute(s[m]==sd,list(sd=-sd))), adj=c(1.2,1) )
        text( xx[max(which(selx<sb))], max(1-prob.extinct), labels=as.expression(substitute(s[b]==sb,list(sb=sb))), adj=c(-0.2,1) )
        arrows( x0=min(xx) + 5*dx, y0=mean(1-prob.extinct), x1=min(xx)+2*dx*m/sqrt(sd) + 5*dx, angle=90, code=3, length=.05 )
        text( min(xx) + 5*dx + (2/2)*dx*m/sqrt(sd), mean(1-prob.extinct), pos=3, labels=expression(2*sigma/sqrt(s)) )
    } )
with( get.prob.extinct(500), {
        plot(xx,1-prob.extinct,type='l',lwd=2,ylim=c(0,max(1-prob.extinct)),xlab="distance (km)", ylab="prob of establishment")
        polygon(c(xx,rev(xx)),c(1-prob.extinct,rep(0,length(xx))),col=adjustcolor("blue",.5))
        # abline(h=1-eqvals[2], lty=2, lwd=2, col='red')
        # abline(v=c(xx[min(which(selx>-sd))],xx[max(which(selx<sb))]),lty=3,lwd=2)
        polygon(xx[c(min(which(selx>-sd)),max(which(selx<sb)))][c(1,1,2,2)], c(0,par("usr")[4])[c(1,2,2,1)],density=5,col=grey(.5))
        text( xx[min(which(selx>-sd))], max(1-prob.extinct), labels=as.expression(substitute(s[m]==sd,list(sd=-sd))), adj=c(.8,1) )
        text( xx[max(which(selx<sb))], max(1-prob.extinct), labels=as.expression(substitute(s[b]==sb,list(sb=sb))), adj=c(0.2,1) )
        arrows( x0=min(xx)+5*dx, y0=mean(1-prob.extinct), x1=min(xx)+2*dx*m/sqrt(sd)+5*dx, angle=90, code=3, length=.05 )
        text( min(xx) + 5*dx + (2/2)*dx*m/sqrt(sd), mean(1-prob.extinct), pos=3, labels=expression(2*sigma/sqrt(s)) )
} )
# legend('topleft',legend=c('prob of establishment','cline location'),lty=c(1,1),col=c('black','green'),lwd=c(2,1),bg='white')
dev.off()



##########

# Helianthus petiolaris/neglectus example:
#  H.p. in southern Colorado; H.n. in Texas-ish

# variance effective population density: total numbers / variance in offspring number
rhovals <- c( 1e4 / 20, 1e5 / 10 )  # per km^2
# dispersal distance: seeds, pollen
sigmavals <- c(.020, .100, 1)  # km / gen
# beneficial selection coefficient within the patches
# --- set so that prob of establishment is 2*sb
sb <- 0.2
# deleterious selection coefficient outside of patches
sm <- .001
# mutation rate (note things will be linear in this
muvals <- c(10^-8,10^-5)
# size of patches is A <- 60 # km^2
# ... but here we want the source area available for new mutations to come from
#    if the dunes are unadapted, i.e. a ring around them of distance, say, 4*sigma
A <- function(sigma) { 4*sigma*2*pi*10 }
# distance between patches
R <- 1000 # km

###
# Table up those values
params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)
exvalues <- t( sapply(1:dim(params)[1], function (k) {
        with( params[k,], everything(mu=mu, rho=rho, sb=sb, sm=sm, sigma=sigma, R=R, A=A(sigma)) )
        } ) )

# write out the table in latex
extable <- xtable(exvalues,digits=0)
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 3
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<.001) | any(abs(x)>1e5) }))] <- -2

filename <- "helianthus-ex-table.tex"
# write("\\documentclass{article} \\usepackage[landscape]{geometry} \\begin{document}", file=filename)
print(extable, file=filename)#, append=TRUE)
# write("\\end{document}", file=filename, append=TRUE)

####
# And look at parallel adaptation within the big species?
source("Spatial_adaptation/standing-variation-fns.R")

# pre-environment-change deleterious selection coefficient
sd <- .01

params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)

exvalues <- t( sapply(1:dim(params)[1], function (k) {
            with(params[k,], everything(mu=mu, rho=rho, sb=sb, sd=sd, sigma=sigma) )
        } ) )

# write out the table in latex
extable <- xtable(exvalues,digits=0)
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 3
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<.001) | any(abs(x)>1e5) }))] <- -2

filename <- "helianthus-standing-ex-table.tex"
# write("\\documentclass{article} \\usepackage[landscape]{geometry} \\begin{document}", file=filename)
print(extable, file=filename)#, append=TRUE)
# write("\\end{document}", file=filename, append=TRUE)


#####
# # Plot against mutation rate  NOT SO USEFUL.
# muvals <- 10^seq.int(-5,-8,length.out=50)
# params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)
# exvalues <- as.data.frame( t( sapply(1:dim(params)[1], function (k) {
#         with( params[k,], everything(mu=mu, rho=rho, sb=sb, sm=sm, sigma=sigma, R=R, A=A) )
#         } ) ) )
# 
# # coplot(ratioMutMig ~ mu | sigma*rho, data=exvalues, overlap=0)
# rhocols <- c("black","red")  # rainbow_hcl(length(rhovals))
# names(rhocols) <- rhovals
# legtext <- as.expression( lapply( rhovals, function(x) substitute( rho==myrho, list(myrho=x) ) ) )
# layout(c(1,length(sigmavals)))
# for (sigmaval in sigmavals) {
#     plot( ratioMutMig ~ mu, data=exvalues[exvalues$sigma==sigmaval,], col=rhocols[paste(rho)], main=as.expression(substitute(sigma==zzz, list(zzz=sigmaval))), xlab=expression(mu), ylab="Mutation/Migration" )
#     if (sigmaval==sigmavals[1]) { legend("topleft",col=rhocols,pch=1, legend=legtext) }
# }
# 


######
# Maize example:
source("patchy-selection-fns.R")

# variance effective population density: total numbers / variance in offspring number
rhovals <- c( 1e4 / 100, 1e5 / 100 )  # per km^2
# dispersal distance: seeds, pollen
sigmavals <- c(.020, .100, 1)  # km / gen
# beneficial selection coefficient within the patches
# --- set so that prob of establishment is 2*sb
sb <- 0.005
# deleterious selection coefficient outside of patches
sm <- .001
# mutation rate (note things will be linear in this
muvals <- c(10^-8,10^-5)
# size of area of influx
A <- 500 # km^2
# distance between patches
R <- 4000 # km
# variance in reproductive number -- sb above is scaled by this.
xisq <- 50

###
# Table up those values
params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)
exvalues <- t( sapply(1:dim(params)[1], function (k) {
        with( params[k,], everything(mu=mu, rho=rho, sb=sb*xisq, sm=sm, sigma=sigma, R=R, A=A) )
        } ) )

# write out the table in latex
extable <- xtable(exvalues,digits=0)
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 3
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<.001) | any(abs(x)>1e5) }))] <- -2

filename <- "maize-ex-table.tex"
# write("\\documentclass{article} \\usepackage[landscape]{geometry} \\begin{document}", file=filename)
print(extable, file=filename)#, append=TRUE)
# write("\\end{document}", file=filename, append=TRUE)

####
# And look at parallel adaptation within the big species?
source("Spatial_adaptation/standing-variation-fns.R")

# pre-environment-change deleterious selection coefficient
sd <- .01

params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)

exvalues <- t( sapply(1:dim(params)[1], function (k) {
            with(params[k,], everything(mu=mu, rho=rho, sb=sb, sd=sd, sigma=sigma) )
        } ) )

# write out the table in latex
extable <- xtable(exvalues,digits=0)
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 3
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<.001) | any(abs(x)>1e5) }))] <- -2

filename <- "maize-standing-ex-table.tex"
# write("\\documentclass{article} \\usepackage[landscape]{geometry} \\begin{document}", file=filename)
print(extable, file=filename)#, append=TRUE)
# write("\\end{document}", file=filename, append=TRUE)

