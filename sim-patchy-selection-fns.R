#!/usr/bin/R
# Forked from spatial-spread.R

###
# Parameters

getsigma <- function (params) {
    # compute sigma (SD of dispersal) from parameters
    probs <- params$m * sapply( params$migrsteps, function (x) x[1] )
    probs <- probs / sum(probs)
    steps <- sapply( params$migrsteps, function (x) sqrt( sum(x[2:3]^2) ) )
    return( sqrt( sum( probs * steps^2 ) ) )
}

getabssigma <- function (params) {
    # compute sigma (SD of dispersal) from parameters
    probs <- params$m * sapply( params$migrsteps, function (x) x[1] )
    probs <- probs / sum(probs)
    steps <- sapply( params$migrsteps, function (x) sqrt( sum(x[2:3]^2) ) )
    return( sum( probs * steps ) )
}

###
# Populations

newpop <- function(params, ntypes=10, nseeds=1 ) {
    # use this to create a new population on a grid of dim'ns range
    # A population is a 3D array, with the first two dimensions indexing spatial location
    #  and the third indexing types
    #  ... so n[i,j,k] is the number of type-k at location (i,j)
    
    n <- array(0, dim=c(params$range,ntypes))
    n[,,2][sample(prod(params$range),nseeds)] <- 1
    n[,,1] <- params$N-n[,,2]
    
    # use this structure to keep track of number of types so far, etc.
    pop <- list(n=n, gen=1, mutations=list(), params=params )
    return(pop)
}


########
# run the simulation

pophistory <- function (pop, nsteps, step=20) {
    # make a 2x1x1-dimensional array of population states,
    # where the fourth dimension is time,
    # and each slice is the population state at time points
    # separated by "step" generations.

    pophist <- array(0,dim=c(dim(pop$n),nsteps))
    pophist[,,,1] <- pop$n
    dimnames(pophist)[[4]]<-step*(1:nsteps)
    
    for (k in 2:nsteps) {
        pop <- generations(pop,step)
        pophist[,,,k] <- pop$n
    }
    return( list(pophist=pophist,pop=pop) )
}


runsim <- function (pop, nsteps, step=10, plotit=TRUE) {
    for (k in 1:nsteps) {
        if (sum(pop$n)==0) { print(k); break }
        pop <- generations(pop, step)
        if (plotit) { plotpop(pop$n) }
    }
    return(pop)
}

generations <- function (pop, ngens) {
        for (k in 1:ngens) { pop <- generation(pop) }
        return(pop)
}

generation <- function (pop) {

        pop$gen <- pop$gen + 1

        # Reproduction
        for (k in 1:dim(pop$n)[3]) {
            # types other than the first type have fitness given by s
            fit <- 1+pop$params$s*(k>1)
            # number of new individuals 
            pop$n[,,k] <- pop$n[,,k] + rrbinom( pop$n[,,k], pop$params$r*fit ) 
        }

        # Mutation
        #   makes some assumptions about there not being too many mutants;
        #   adds new mutants to as-yet-unfilled types, if possible.
        if (pop$params$mu>0) {
            # compute population proportions
            areas <- apply(pop$n,3,sum)/pop$params$N
            # use this is deciding which type to make new mutations
            smallAreas <- ( areas<max(quantile(areas,.1),.01) )
            if ( all( !smallAreas ) ) {  # if none are small put the mutations anywhere
                smallAreas <- rep(TRUE, length(smallAreas))
            }
            for (k in 1:dim(pop$n)[3]) {
                nmuts <- rbinom(1, sum(pop$n[,,k]), pop$params$mu)  # how many mutants are produced from type k
                if (nmuts>maxmuts) { stop("Too many mutations!") }  # can't deal with too many
                if (nmuts>0) for (ell in 1:nmuts) {
                    newtype <- sample((1:dim(pop$n)[3])[smallAreas],1)  # what type is it
                    newloc <- which(rmultinom(1,1,pop$n[,,k]) > 0)      # where is it
                    pop$n[,,newtype][newloc] <- pop$n[,,newtype][newloc] + 1  # add it in
                    pop$mutations <- c(pop$mutations, list(c(pop$gen,k,newtype,newloc))) # record this event
                }
            }
        }

        # Migration
        #  -- also an approximation: this allows some to migrate
        #   more than once, but with probability < m^2.
        for (k in 1:dim(pop$n)[3]) {
            for (migr in pop$params$migrsteps) {
                pop$n[,,k] <- migrate(pop$n[,,k,drop=FALSE], pop$params$m, migr)
            }
        }

        # Resample down to N
        #   (approximate hypergeometric by multinomial)
        pop$n <- aperm( apply(pop$n, c(1,2), function (nn) { rmultinom(1, pop$params$N, nn) } ), c(2,3,1))

        return(pop)
}

migrate <- function ( n, m, migr ) {
    # does a single migration step on a matrix n.
    # with total probability per individual of migration m.
    # expects migr to be a triplet of the format (probability, dx, dy)

    prob <- migr[1]*m
    dx <- migr[2]
    if (length(migr)==2) { dy <- 0 } else { dy <- migr[3] }
    if ( dim(n)[3] != 1 ) { 
        warn( "Passing multiple types to migrate!" ) 
    } else { dim(n) <- dim(n)[1:2] }

    migrants <- rrbinom( n, prob )
    outmigrants <- migrants  # where they came from (need to subtract this off)

    nc <- dim(n)[2] # x-extent of populations
    nr <- dim(n)[1] # y-extent of populations

    if (dx>0) {
    	migrants <- cbind( matrix(0,nrow=dim(migrants)[1],ncol=dx), migrants[,-((nc-dx+1):nc),drop=FALSE] )
    } else if (dx<0) {
    	migrants <- cbind( migrants[,-1:dx,drop=FALSE], matrix(0,nrow=dim(migrants)[1],ncol=-dx) )
    }
    if (dy>0) {
    	migrants <- rbind( matrix(0,ncol=dim(migrants)[2],nrow=dy), migrants[-((nr-dy+1):nr),,drop=FALSE] )
    } else if (dy<0) {
    	migrants <- rbind( migrants[-1:dy,,drop=FALSE], matrix(0,ncol=dim(migrants)[2],nrow=-dy) )
    }

    #while (dx > 0) { migrants <- cbind(0,migrants[,-nc,drop=FALSE]); dx <- dx-1 } # shift right
    #while (dx < 0) { migrants <- cbind(migrants[,-1,drop=FALSE],0); dx <- dx+1 } # shift left

    #while (dy > 0) { migrants <- rbind(0,migrants[-nr,,drop=FALSE]); dy <- dy-1 } # shift down
    #while (dy < 0) { migrants <- rbind(migrants[-1,,drop=FALSE],0); dy <- dy+1 } # shift up

    return( n + migrants - outmigrants )
}


####
# Plotting

colorRampTrans <- function (basecol, n=12) {
    # return a vector of n colors
    # ranging from tranparent-white (at the start) to the color (at the end)
    # like colorRamp does.

    ramp <- colorRamp(c("#FFFFFFFF",basecol), space="Lab")
    x <- ramp(seq.int(0, 1, length.out = n))
    x <- cbind(x, seq.int(0,255,length.out = n) )
    rgb(x[, 1], x[, 2], x[, 3], x[,4], maxColorValue = 255)
}


plotpophist <- function(pophist,timeslice=TRUE,...) {
    # Plot one-dimensional population history with time on the x-axis
    if (!any(dim(pophist$pophist)[1:2]==1)) { stop("Population is not one-dimensional.") } 
    zimg <- plotpop( aperm(pophist$pophist[,,,timeslice,drop=TRUE],c(3,1,2)), params=pophist$pop$params, xlab="time", ylab="space", ... )
    axis(1); axis(2)
    return( invisible(zimg) )
}


animpophist <- function (pophist,delay=0,nsteps=120,...) {
    ptimes <- floor(seq(1,dim(pophist$pophist)[4],length.out=min(nsteps,dim(pophist)[4])))
    # plot the population history as an animation
    if ( dim(pophist$pophist)[1]==1 ) {
        lapply( ptimes, function (t) { plotlinearpop(pophist$pophist[,,,t],params=pophist$pop$params,...); Sys.sleep(delay) } )
    } else {
        # apply(pophist$pophist,4,function(x) { plotpop(x,params=pophist$pop$params,plotlegend=FALSE,...); Sys.sleep(delay) } )
        lapply( ptimes, function (t) { plotpop(pophist$pophist[,,,t],params=pophist$pop$params,plotlegend=FALSE,...); Sys.sleep(delay) } )
    }
    return(invisible(dim(pophist$pophist)[4]))
}


plotlinearpop <- function(n, thresh=.05, cols, ...) {
    # plot (colored) curves representing proportions.
    if ("n" %in% names(n)) { n <- n$n }
    if (length(dim(n))==4) { n <- n[,,,,drop=TRUE] }
    else if (length(dim(n))==3) { n <- n[,,,drop=TRUE] }
    if ( length(dim(n))!=2 ) { stop("plotlinearpop: incorrect dimensions") }
    rangelen <- dim(n)[1]
    if (missing(cols)) { cols <- c(list(c("#FFFFFF")), lapply(rainbow(dim(n)[2]), colorRampTrans, n=12) ) }
    
    # renormalize to proportions
    n <- n/apply(n,1,sum)

    plot( 1, 1, type="n", xlim=c(1,rangelen), ylim=c(0,1) )
    for (k in 2:dim(n)[2]) {
        pcolor <- cols[[k]][length(cols[[k]])]
        lines( n[,k], col=pcolor )
    }
}

plotpop <- function (n, params, maxtimes=200, thresh=.05, cols, coltrans, ...) {
    # Layer the colors on top of eachother (multiplicatively) 
    # for the different types

    if (missing(params)) { params <- n$params }
    if ("n" %in% names(n)) { n <- n$n }
    if (length(dim(n))==4 & dim(n)[4]==1) { dim(n) <- dim(n)[1:3] }

    # subsample if too many columns
    if ((dim(n)[1]>maxtimes)) { n <- n[(1:(dim(n)[1]))%%floor(dim(n)[1]/maxtimes)==0,,] }
    if ((dim(n)[2]>maxtimes)) { n <- n[,(1:(dim(n)[2]))%%floor(dim(n)[2]/maxtimes)==0,] }
    # color picking
    if (missing(cols)) { cols <- c(list(c("#FFFFFF")), lapply(rainbow(dim(n)[3]), colorRampTrans, n=12) ) }
    if (missing(coltrans)) { coltrans <- function (x,colors,ncols=length(colors)) { colors[ceiling((ncols-1)*x/params$N+1)] } }

    # Bulk plotting
    img <- array(1,dim=c(3,dim(n)[1:2]))  # has r,g,b layers above population
    for (k in 2:dim(n)[3]) {
        # image(1-log(1+n[,,k])/log(N),zlim=c(0,1), col=cols[[k]], new=k>1)
        # image(n[,,k]/N,zlim=c(0,1), col=cols[[k]], new=k>1)
        newcolors <- sapply(coltrans(n[,,k],cols[[k]]),col2rgb)/255
        dim(newcolors) <- dim(img)
        img <- img * newcolors
    }
    # img <- apply(img,c(2,3),function(x)rgb(x[1],x[2],x[3]))
    img <- rgb(img[1,,],img[2,,],img[3,,])
    dim(img) <- dim(n)[1:2]
    plotpopcircles(img,dim(n),colors=sapply(cols,function(x)x[1]), ...)

    # Add points for transient populations
    #for (k in 2:dim(n)[3]) {
    #    if ( any(n[,,k] > 0) & all(n[,,k] < N*thresh) ) {
    #        mutlocs <- which(n[,,k] > 0)
    #        points( row(n[,,k])[mutlocs], col(n[,,k])[mutlocs], col=cols[[k]][floor(length(cols[[k]])/2)], pch="." )
    #    }
    #}

    return(invisible(img))
}

plotpopcircles <- function(img,dims=dim(img),colors,plotlegend=TRUE,...) {
    # Do the plotting.
    cex <- 3.5*min( par("pin") / (par("cin")*dims[1:2]) )
    plot(row(img),col(img),col=img,pch=20,cex=cex,xaxt="n",yaxt="n",xlim=c(1,(dims[1]+1)), ...)
    if (plotlegend) {
        legend("topright",legend=1:length(colors),pch=20,col=sapply(colors,tail,n=1))
    }
}



#####
# Helper functions

rrbinom <- function (sizes, prob) {
    # a helper function to vectorize rbinom
    # in the way we use it here most often:
    #
    # for each element of sizes, 
    # draws a single binomial with probability prob
    # respecting the properties of sizes.
    
    dims <- dim(sizes)
    dimn <- dimnames(sizes)
    
    sizes <- rbinom(length(sizes), sizes, prob)
    dim(sizes) <- dims
    dimnames(sizes) <- dimn
    
    return(sizes)
}
