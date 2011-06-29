#
# Implement theory from patchy selection parts

# mu is mutation rate
# rho is population density
# sigma is mean squared dispersal distance
# sb >0 is advantage of mutation within patches
# sd >0 is disadvantage of mutation between patches
# R is distance between patches
# 

paramString <- function(ps=c("mu","rho","sb","sd","sigma","R"), pos=-1) {
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

