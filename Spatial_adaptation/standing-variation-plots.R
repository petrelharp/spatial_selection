
source("Spatial_adaptation/standing-variation-fns.R")

# how big to make the figures
figwidth <- 10
figheight <- 8

# playing around
mu <- 10^(-8)
rho <- 10
sb <- .01
sd <- .01
sigmavals <- c(0.1,1,10,20,50)

# plot Standing Proportion Area as a function of sigma
# sPvals <- sapply(sigmavals, function (sigma) { standingProportionArea(mu,rho,sb,sd,sigma)$value } )
# plot(sigmavals, sPvals, type="b", xlab=expression(paste(sigma, " = dispersal SD")), ylab=expression(paste(z[0], " = proportion of patches from standing variation")), ylim=c(0,1))


# sickle cell case from ralphcoop2010
mu <- 10^(-8)
rhovals <- c(2.5,25)
rhocols <- c("black","red")
rholabs <- c(expression(rho==2.5), expression(rho==25))
sb <- .05
sd <- .3
sigmavals <- 1+(1:20)*30
sigmacols <- rainbow(length(sigmavals)+5)[1:length(sigmavals)]
sigmalabs <- c(expression(sigma==0.1), expression(sigma==1.0), expression(sigma==10), expression(sigma==20), expression(sigma==50), expression(sigma==100))
otherlabs <- expression(paste( mu == 10^-8, "  ", s[b] == .05, "  ", s[d] == .3 ))

# plot Proportion Area as a function of sigma
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="proportion-by-sigma.eps", width=figwidth, height=figheight, title=paramString() )
plot(0, 0, type="n", xlab=expression(paste(sigma, " = dispersal SD")), ylab="Proportion from standing variation", xlim=range(sigmavals), ylim=c(0,1), main=otherlabs )
legend("topleft", lty=c(rep(1,length(rhocols)),1,2), col=rhocols, legend=c(rholabs, expression(paste(z[0], " = proportion area")), expression(paste(nu["+"], " = proportion numbers"))) )
for (k in 1:length(rhovals)) {
    rho <- rhovals[k]
    sPvals <- sapply(sigmavals, function (sigma) { standingProportionArea(mu,rho,sb,sd,sigma)$value } )
    lines(sigmavals, sPvals, col=rhocols[k])
    sPvals <- sapply(sigmavals, function (sigma) { standingProportionNumbers(mu,rho,sb,sd,sigma)$value } )
    lines(sigmavals, sPvals, col=rhocols[k], lty=2)
}
dev.off()

# sickle cell case from ralphcoop2010
mu <- 10^(-5)
rho <- 25
sb <- .05
sigmavals <- c(0.1,1,10,20,50,100)
sigmacols <- rainbow(length(sigmavals)+5)[1:length(sigmavals)]
sigmalabs <- c(expression(sigma==0.1), expression(sigma==1.0), expression(sigma==10), expression(sigma==20), expression(sigma==50), expression(sigma==100))
sdvals <- (0:50)/50
otherlabs <- expression(paste( mu == 10^-5, "  ", s[b] == .05, "  ", rho == 25 ))

# plot Proportion Area as a function of s_d
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="proportion-by-sd.eps", width=figwidth, height=figheight, title=paramString() )
par(mar=par("mar")*c(1,1,1,2))
plot(0, 0, type="n", xlab=expression(paste(s[d], " = disadvantage")), ylab="Proportion from standing variation", xlim=range(sdvals), ylim=c(0,1), main=otherlabs, bty="L" )
legend(1,1, lty=c(rep(1,length(sigmacols)),1,2), col=sigmacols, legend=c(sigmalabs, expression(paste(z[0], " = prop area")), expression(paste(nu["+"], " = prop numbers"))), xpd=TRUE, cex=0.5 )
for (k in 1:length(sigmavals)) {
    sigma <- sigmavals[k]
    sPvals <- sapply(sdvals, function (sd) { standingProportionArea(mu,rho,sb,sd,sigma)$value } )
    lines(sdvals, sPvals, col=sigmacols[k])
    sPvals <- sapply(sdvals, function (sd) { standingProportionNumbers(mu,rho,sb,sd,sigma)$value } )
    lines(sdvals, sPvals, col=sigmacols[k], lty=2)
}
dev.off()

# sickle cell case from ralphcoop2010
mu <- 10^(-5)
rho <- 25
sb <- .05
sigmavals <- c(0.1,1,10,20,50,100)
sigmacols <- rainbow(length(sigmavals)+5)[1:length(sigmavals)]
sigmalabs <- c(expression(sigma==0.1), expression(sigma==1.0), expression(sigma==10), expression(sigma==20), expression(sigma==50), expression(sigma==100))
sdvals <- (0:50)/50
otherlabs <- expression(paste( mu == 10^-5, "  ", s[b] == .05, "  ", rho == 25 ))

# plot characteristic length as a function of sd
# with many sigma values
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="charlen-by-sd-sigma.eps", width=figwidth, height=figheight, title=paramString() )
plot(0, 0, type="n", xlab=expression(paste(s[d], " = disadvantage")), ylab=expression(paste(chi, " = characteristic length")), xlim=range(sdvals), ylim=c(0,100), main=otherlabs )
legend("topleft", lty=1, col=sigmacols, legend=sigmalabs)
for (k in 1:length(sigmavals)) {
    sigma <- sigmavals[k]
    CLvals <- sapply(sdvals, function (sd) { charLength(mu,rho,sb,sd,sigma)$value } )
    lines(sdvals, CLvals, col=sigmacols[k])
}
dev.off()

#with one sigma value
mu <- 10^(-5)
rho <- 25
sb <- .05
sigma <- 50
otherlabs <- expression(paste( mu == 10^-5, "  ", s[b] == .05, "  ", rho == 25, "  ", sigma == 50 ))

setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="charlen-by-sd-limit.eps", width=figwidth, height=figheight, title=paramString() )
plot(0, 0, type="n", xlab=expression(paste(s[d], " = disadvantage")), ylab=expression(paste(chi, " = characteristic length")), xlim=range(sdvals), ylim=c(0,100), main=otherlabs )
legend("topleft", lty=1:3, legend=c("characteristic length","only standing variation", "only new mutations"))
CLvals <- sapply(sdvals, function (sd) { charLength(mu,rho,sb,sd,sigma)$value } )
lines(sdvals, CLvals)
lines( sdvals, sapply(sdvals, function(sd) { oldCharLength(mu,rho,sb,sd,sigma)$value } ), lty=2 )
abline(h=newCharLength(mu,rho,sb,sd,sigma), lty=3)
dev.off()


## contour plots
mu <- 10^(-6)
rho <- 25
sb <- .05

# plot characteristic length versus sd and sigma
sigmavals <- 1+(0:30)*4
sdvals <- (0:30)/30
CLvals <- apply( expand.grid(sigmavals,sdvals), 1, function (x) { sigma<-x[1]; sd<-x[2]; charLength(mu,rho,sb,sd,sigma)$value } )
dim(CLvals) <- c(length(sigmavals),length(sdvals))
#
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="charlen-contour.eps", width=figwidth, height=figheight, title=paramString() )
filled.contour(sdvals, sigmavals,CLvals,xlab=expression(paste(s[d], " = disadvantage")), ylab=expression(paste(sigma, " = dispersal SD")), main=expression(paste(chi, " = characteristic length")), sub=paramString(list(mu=mu,rho=rho,sb=sb)))
dev.off()

# plot proportion versus sd and sigma
sigmavals <- 1+(0:30)*4
sdvals <- (0:30)/30
SPvals <- apply( expand.grid(sigmavals,sdvals), 1, function (x) { sigma<-x[1]; sd<-x[2]; standingProportionArea(mu,rho,sb,sd,sigma)$value } )
dim(SPvals) <- c(length(sigmavals),length(sdvals))
#
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="proportion-contour.eps", width=figwidth, height=figheight, title=paramString() )
filled.contour(sdvals, sigmavals,SPvals,xlab=expression(paste(s[d], " = disadvantage")), ylab=expression(paste(sigma, " = dispersal SD")), main=expression(paste(z[0], " = proportion area from standing variation")), sub=paramString(list(mu=mu,rho=rho,sb=sb)) )
dev.off()
