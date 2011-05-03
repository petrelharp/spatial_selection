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

charLength <- function (mu, rho, sb, sd, sigma) {
    # the characteristic length, which is
    #   the positive solution to
    #   lambda0 x^2 + lambda x^3 / v = 1/pi
    roots <-  polyroot( c(-1/pi, 0, 2*mu*rho*sb/sd, sqrt(2*sb)*mu*rho / sigma ) )
    # sanity check
    ii <- which.max(Re(roots))
    if (Im(roots[ii])>1000*.Machine$double.eps) { warning("Imaginary root in characteristic length!") }
    return( Re(roots[ii]) )
}

paramString <- function(ps=c("mu","rho","sb","sd","sigma"), pos=-1) {
    # return a string for use in plots of the parameters.
    # takes a vector or list of parameter names
    if ( is.null(names(ps)) ) { names(ps) <- ps }
    if ( mode(ps) != mode(1.0) ) { ps[1:length(ps)] <- 0; mode(ps) <- mode(1.0) }
    for ( x in names(ps) ) { 
        if (exists(x,where=pos,mode=mode(1.0))) {
            ps[x] <- get(x, pos=pos, mode=mode(1.0)) 
        } else {
            ps[x] <- NA
        }
    }
    paste( sapply( 1:length(ps), function (k) { sprintf( "%s=%.3g", names(ps)[k], ps[k] ) } ), collapse=", " )
}

