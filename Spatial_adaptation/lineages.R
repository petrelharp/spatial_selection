
# Given a population history, trace lineages back in time
#  and record statistics about their path

lineages <- function ( pophist, nlin, migrsteps, T=dim(pophist)[4], m=1, linit=NULL ) {

    # simulate the movement of nlin selected lineages back through pophist,
    # beginning at locations in linit 
    # using the dispersal pattern migrsteps
    #  migrsteps is of the form list( (prob, dx, dy) ) with sum(prob)=1.
    #  and m is the probability of migrating
    #  pophist is indexed by (x, y, type, time)
    #    with types either 0 (nonselected) or 1 (selected)
    # returns:
    #  LL = location of lineages through time (2 x T x nlin)

    # NOTE y comes before x
    ny <- dim(pophist)[1]
    nx <- dim(pophist)[2]

    if (is.null(linit)) { # sample initial locations uniformly
        initpop <- pophist[,,2,T]
        dim(initpop) <- dim(pophist)[1:2]
        locs <- sample( nx*ny, nlin, replace=TRUE, prob=initpop )
        linit <- rbind( row(initpop)[locs], col(initpop)[locs] )
    }

    # add a zero step and reweight migrsteps
    mmm <- sum( sapply( migrsteps, function (x) x[1] ) )
    migrsteps <- c( list( c(1-m, 0, 0) ), lapply(migrsteps, function (x) { c( x[1]*m/mmm, x[-1] ) } ) )

    # the Lineage Locations
    #   and the Local Sizes (migration-weighted)
    LL <- array(NA, c(2,T,nlin))
    LS <- array(NA, c(1,T,nlin))
    LL[,1,] <- linit
    LS[,1,] <- pophist[cbind(linit[1,],linit[2,],2,T)]  # ?'[': in x[i], 'i' can be a matrix with as many columns as there are dimensions of ‘x’...
    dimnames(LS)[[2]] <- dimnames(LL)[[2]] <- dimnames(pophist)[[4]]  # label the time entries
    dimnames(LL)[[1]] <- c("y","x")
    for (t in 2:T) {
        # pophist[,,2,T-t+1,drop=FALSE] # is the distribution of selected types in previous generation
        for (k in 1:nlin) {
            # problocs[1,] is vector of weights for where to migrate to
            # problocs[2:3,] is the set of locations
            problocs <- sapply(migrsteps, function (migr) {
                        xy <- LL[,t-1,k] + migr[3:2]   # LL stores (y,x) not (x,y)...
                        return( c(pophist[xy[1],xy[2],2,T-t+1]*m[1],xy) )
                    } )
            LS[,t,k] <- sum(problocs[1,])
            if ( LS[,t,k] <= 0 ) { warning("Ooops!  Nowhere to migrate to.  Aborting."); break }
            LL[,t,k] <- problocs[-1, sample(length(migrsteps),1,prob=problocs[1,])]
        }
    }
    return( list( lineagelocs=LL, localsizes=LS ) )
}


