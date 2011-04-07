
source("Spatial_adaptation/standing-variation-fns.R")

# playing around
mu <- 10^(-8)
rho <- 10
sb <- .01
sd <- .01
sigmavals <- c(0.1,1,10,20,50)

# plot Standing Proportion Numbers as a function of sigma
# sPvals <- sapply(sigmavals, function (sigma) { standingProportionNumbers(mu,rho,sb,sd,sigma)$value } )
# plot(sigmavals, sPvals, type="b", xlab=expression(paste(sigma, " = dispersal SD")), ylab=expression(paste(z[0], " = proportion of patches from standing variation")), ylim=c(0,1))


# sickle cell case from ralphcoop2010
mu <- 10^(-8)
rhovals <- c(2.5,25)
rhocols <- c("black","red")
rholabs <- c(expression(rho==2.5), expression(rho==25))
sb <- .05
sd <- .01
sigmavals <- c(0.1,1,10,20,50,100)
otherlabs <- expression(paste( mu == 10^-8, "  ", s[b] == .05, "  ", s[d] == .01 ))

# plot Standing Proportion Numbers as a function of sigma
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript( file="proportion-by-sigma.eps", width=6.5, height=5, title=paramString() )
plot(0, 0, type="n", xlab=expression(paste(sigma, " = dispersal SD")), ylab=expression(paste(z[0], " = proportion of patches from standing variation")), xlim=range(sigmavals), ylim=c(0,1), main=otherlabs )
legend("topleft", lty=1, col=rhocols, legend=rholabs)
for (k in 1:length(rhovals)) {
    rho <- rhovals[k]
    sPvals <- sapply(sigmavals, function (sigma) { standingProportionNumbers(mu,rho,sb,sd,sigma)$value } )
    lines(sigmavals, sPvals, col=rhocols[k])
}

