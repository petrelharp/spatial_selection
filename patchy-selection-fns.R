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
    return( sigma / sqrt(2*s*gm) )
}

minimumPatch <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # the minimal patch area ( = width^2)
    return( ( 2 * atanh( sqrt(gm/gb) ) * sigma / sqrt(2*s*gm) )^2 )
}

mutationInflux <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # approximation for mutational influx per patch
    # XXX to-do: do numerics?
    return( 2 * s*gb * A * rho * mu )
}

migrationInflux <- function(mu, rho, s, gb, gm, sigma, R, A) {
    # approximation for migrational influx per patch
    # XXX to-do: do numerics?
    # return( 2 * s*gb * A * rho * exp( - sqrt(abs(s*gm))*R/sigma ) )
    return( s*gb * sqrt(s*gm/pi) * rho * exp( - sqrt(abs(s*gm))*R/sigma ) )
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

logseq <- function(from,to,...) { 
    # returns a sequence whose logarithms are uniformly spaced,
    #  i.e. a multiplicative sequence
    exp( seq( log(from), log(to), ... ) ) 
}
expseq <- function(from,to,...) { 
    # returns a sequence that when exponentiated is uniformly spaced
    log( seq( exp(from), exp(to), ... ) ) 
}

require(colorspace)
hcolor <- function (z,alpha=.75,cols=adjustcolor(heat_hcl(64,h=c(40,360),l=70,c.=c(70,100)),alpha),nc=length(cols),...) {
    # coloring function for hplot
    cols[as.numeric(cut(pmin(1,pmax(-1,z)),breaks=seq(-1,1,length.out=nc+1),include.lowest=TRUE))]
}
hplot <- function (z,x=as.vector(col(z)),y=as.vector(row(z)),scale=1,max.cex=4,labs,xlabs=labs,ylabs=labs,alpha=.75,...) {
    # Like heatmap but with circles...
    plot( x=x, y=y, pch=20, cex=max.cex*pmin(1,sqrt(abs(z)/scale)), col=hcolor(z/scale,alpha=alpha), xaxt='n', yaxt='n', xlab="", ylab="", ... )
    if (!missing(labs) | !missing(xlabs)) 
        axis(1, at=1:length(xlabs), labels=xlabs, las=2) 
    if (!missing(labs) | !missing(ylabs)) 
        axis(2, at=1:length(ylabs), labels=ylabs, las=2)
}

