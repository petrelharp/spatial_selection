#
# Implement theory from standing variation parts

# mu is mutation rate
# rho is population density
# sigma is mean squared dispersal distance
# sb >0 is advantage of mutation for t>0
# sd >0 is disadvantage of mutation for t<0
# 
# ... and
# lambda <- 2 * mu * rho * sb    # is rate of new successful mutations
# lambda0 <- lambda / sd         # is approximate rate of standing successful mutations
# v <- sigma * sqrt(2*sb)        # speed of wave of advance

meanTime <- function (mu, rho, sb, sd, sigma) {
    # mean time til adaptation of a point;
    # E[\tau] is int_0^\infty exp(-lambda0 pi v^2 t^2 - lambda pi v^2 t^3) dt
    f <- function (t) {
        exp( - 4 * mu * rho * sb^2 * sigma^2 * t^2 * (1/sd + t) )
    }
    return( integrate(f, 0, Inf) )
}

meanNewDensity <- function (mu, rho, sb, sd, sigma) {
    # mean density of patches arising from new mutations
    #  is lambda * E[\tau]
    x <- meanTime(mu, rho, sb, sd, sigma)
    x$value <- x$value * 2 * mu * rho * sb 
    return( x )
}

standingProportionArea <- function (mu, rho, sb, sd, sigma) {
    # proportion of space arising from standing variation
    # is 2 * lambda0 * pi * v * t * exp( - 2 * lambda0 * pi * v * t - lambda * pi * v^2 * t^2 )
    f <- function (t) {
        (1/sd) * 4 * sqrt(2) * mu * rho * sb^(3/2) * pi * t * exp( - 4 * mu * rho * sb^(3/2) * sigma * pi * t ( sqrt(2) / sd + sigma * t ) )
    }
    return( integrate(f, 0, Inf) )
}

standingProportionNumbers <- function (mu, rho, sb, sd, sigma) {
    # proportion of patches arising from standing variation
    #    is 1/(1+sd * E[\tau])
    x <- meanTime(mu, rho, sb, sd, sigma)
    x$value <- 1 / (1 + sd * x$value)
    return( x )
}

# playing around
mu <- 10^(-8)
rho <- 10
sb <- .01
sd <- .01
sigmavals <- c(0.1,1,10,20,50)

# plot Standing Proportion Numbers as a function of sigma
sPvals <- sapply(sigmavals, function (sigma) { standingProportionNumbers(mu,rho,sb,sd,sigma)$value } )
plot(sigmavals, sPvals, type="b", xlab=expression(paste(sigma, " = dispersal SD")), ylab=expression(paste(z[0], " = proportion of patches from standing variation")), ylim=c(0,1))


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
plot(0, 0, type="n", xlab=expression(paste(sigma, " = dispersal SD")), ylab=expression(paste(z[0], " = proportion of patches from standing variation")), xlim=range(sigmavals), ylim=c(0,1), main=otherlabs )
legend("topleft", lty=1, col=rhocols, legend=rholabs)
for (k in 1:length(rhovals)) {
    rho <- rhovals[k]
    sPvals <- sapply(sigmavals, function (sigma) { standingProportionNumbers(mu,rho,sb,sd,sigma)$value } )
    lines(sigmavals, sPvals, col=rhocols[k])
}

