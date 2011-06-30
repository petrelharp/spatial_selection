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
# lambda0 <- lambda / log(1/(1-sd))         # is approximate rate of standing successful mutations
# v <- sigma * sqrt(2*sb)        # speed of wave of advance

everything <- function (mu, rho, sb, sd, sigma) {
    # bundle up everything
        return( c( 
            c( z0 = standingProportionArea(mu,rho,sb,sd,sigma)$value,
                chi = charLength(mu,rho,sb,sd,sigma)$value,
                Etau = meanTime(mu,rho,sb,sd,sigma)$value),
            c( mu=mu, rho=rho, sb=sb, sd=sd, sigma=sigma )
        ) )
}

meanTime <- function (mu, rho, sb, sd, sigma) {
    # mean time til adaptation of a point;
    # E[\tau] is int_0^\infty exp(-lambda0 pi v^2 t^2 - lambda pi v^2 t^3) dt
    if (sd==0) {
        return( list(value=0) )
    } else if (sd==1) {
        f <- function (t) {
            exp( - 4 * mu * rho * sb^2 * sigma^2 * t^3 / 3) 
        }
    } else { 
        f <- function (t) {
            exp( - 4 * mu * rho * sb^2 * sigma^2 * t^2 * (1/log(1/(1-sd)) + t/3) )
        }
    }
    xx <- 2^(-30:30)
    yy <- sapply(xx, function(x) { f(x) } )
    scale <- xx[ which.min( abs(10-yy*xx) ) ]
    return( integrate(function(x) { f(x*scale)*scale }, 0, Inf) )
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
    # is int_0^\infty 2 * lambda0 * pi * v^2 * t * exp( - lambda0 * pi * v^2 * t^2 - lambda * pi * v^2 * t^3 /3 ) dt
    if (sd==1) {
        return( list(value=0) )
    } else if (sd==0) {
        return( list(value=1) )
    } else { 
        lambda <- 2*mu*rho*sb
        lambdaoh <- lambda/ log(1/(1-sd))
        v <- sigma * sqrt(2*sb)
        f <- function (t) {
            2 * lambdaoh * pi * v^2 * t * exp( - lambdaoh*pi*v^2*t^2 - lambda*pi*v^2*t^3 /3 )
        }
        xx <- 2^(-30:30)
        yy <- sapply(xx, function(x) { f(x) } )
        scale <- xx[ which.min( abs(10-yy*xx) ) ]
        return( integrate(function(x) { f(x*scale)*scale }, 0, Inf) )
    }
}

standingProportionNumbers <- function (mu, rho, sb, sd, sigma) {
    # proportion of patches arising from standing variation
    #    is 1/(1+log(1/(1-sd)) * E[\tau])
    if (sd==1) {
        return( list(value=0) )
    } else { 
        x <- meanTime(mu, rho, sb, sd, sigma)
        x$value <- 1 / (1 + log(1/(1-sd)) * x$value)
        return( x )
    }
}


charLength <- function (mu, rho, sb, sd, sigma) {
    # the characteristic length, which is
    #   the positive solution to
    #   lambda0 x^2 + lambda x^3 / v = 1/pi
    if (sd==1) {
        return( list( value=(sigma/(rho*mu*pi*sqrt(2*sb)))^(1/3), roots=NA) )
    } else if (sd==0) {
        return( list( value=0, roots=NA) )
    } else {
        roots <-  polyroot( c(-1/pi, 0, 2*mu*rho*sb/log(1/(1-sd)), sqrt(2*sb)*mu*rho / sigma ) )
        ii <- which.max(Re(roots[Im(roots)<10*.Machine$double.eps]))
        # sanity check
        if (length(ii)==0) { stop("No real roots for characteristic length?") }
        return( list(value=Re(roots[ii]), roots=roots) )
    }
}

newCharLength <- function (mu, rho, sb, sd, sigma) {
    # the characteristic length *if only looking at new mutations*, which is
    #     x = ( v/(pi lambda) )^(1/3)
    return( list( value=(sigma/(rho*mu*pi*sqrt(2*sb)))^(1/3), roots=NA) )
}

oldCharLength <-  function (mu, rho, sb, sd, sigma) {
    # the characteristic length *if only looking at standing variation*, which is
    #     x = ( 1/(pi lambda0) )^(1/2)
    if (sd==0) {
        return( list( value=0, roots=NA) )
    } else {
        return( list( value=(2*mu*rho*sb/log(1/(1-sd)))^(-1/2), roots=NA) )
    }
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
    paste( sapply( 1:length(ps), function (k) { sprintf( "%s=%.2g", names(ps)[k], ps[k] ) } ), collapse=", " )
}

