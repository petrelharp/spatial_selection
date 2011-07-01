#
# Implement theory from patchy selection parts

# mu is mutation rate
# rho is population density
# sigma is mean squared dispersal distance
# sb >0 is advantage of mutation within patches
# sm <0 is disadvantage of mutation between patches
# R is distance between patches
# A is area of patches
# 
# rearrange so that
# s = max(|sb|,|sm|)
# sb = gb * s
# sm = gm * s

everything <- function(mu, rho, sb, sm, sigma, R, A) {
    # compute all these things
    s <- max(sb,sm)
    gb <- sb/s
    gm <- sm/s
    mutI <- mutationInflux(mu, rho, s, gb, gm, sigma, R, A)
    migI <- migrationInflux(mu, rho, s, gb, gm, sigma, R, A)
    return( c( c( clScale=clineScale(mu, rho, s, gb, gm, sigma, R, A),
        minPatch=minimumPatch(mu, rho, s, gb, gm, sigma, R, A),
        mutIn=mutI,
        migIn=migI,
        ratioMutMig = mutI/migI
        ), c(mu=mu, rho=rho, sb=sb, sm=sm, sigma=sigma, R=R, A=A)
        ) )
}

clineScale <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # sigma / sqrt(2 s) = the scale on which things happen between patches; 
    # see slatkin; barton etc.
    return( sigma / sqrt(2*sm) )
}

minimumPatch <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # the minimal patch area ( = width^2)
    return( ( 2 * atanh( sqrt(gm/gb) ) * sigma / sqrt(2*sm) )^2 )
}

mutationInflux <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # approximation for mutational influx per patch
    # XXX to-do: do numerics?
    return( 2 * sb * A * rho * mu )
}

migrationInflux <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # approximation for migrational influx per patch
    # XXX to-do: do numerics?
    return( 2 * sb * A * rho * exp( - sqrt(abs(s*gm)*R/sigma) ) )
}


paramString <- function(ps=c("mu","rho","s","gb","gm","sigma","R", "A"), pos=-1) {
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

