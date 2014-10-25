#!/usr/bin/R
# Forked from spatial-spread.R

###
# Parameters

getsigma <- function (params) {
    # compute sigma (SD of dispersal) from parameters
    probs <- sapply( params$migrsteps, function (x) x[1] )
    probs <- params$m * probs / sum(probs)
    stepsq <- sapply( params$migrsteps, function (x) ( sum(x[2:3]^2) ) )
    return( sqrt( sum( probs * stepsq ) ) )
}

getgrowth <- function (params) {
    # growth rate, i.e. numbers behave as exp(t*g):
    return( list( gb=log((1+params$r*(1+params$sb))/(1+params$r)), gm=log((1+params$r*(1+params$sm))/(1+params$r)) ) )
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

pophistory <- function (pop, nsteps, start=1, step=1, progress=0, burnin=0) {
    # make a 2+1+1-dimensional array of population states,
    # where the fourth dimension is time,
    # and each slice is the population state at time points
    # separated by "step" generations.
    if (pop$gen != start) { 
        cat("Setting generation to", start, ".\n")
        pop$gen <- start
    }
    pophist <- array(0,dim=c(dim(pop$n),nsteps))
    pophist[,,,1] <- pop$n
    attr(pophist,"gens") <- start + (step*(1:nsteps)-1)
    occupation <- array(0,dim=dim(pop$n))
    
    for (k in 2:(nsteps*step)) {
        pop <- generation(pop)
        if (k>burnin) { occupation <- occupation + pop$n }
        if (k %% step == 0) { pophist[,,,k%/%step] <- pop$n }
        if (progress>0 & k %% progress == 0) { cat(paste("..",k,"\n")) }
    }
    return( list(pophist=pophist,pop=pop,occupation=occupation,burnin=burnin) )
}


runsim <- function (pop, nsteps, step=10, plotit=TRUE) {
    for (k in 1:nsteps) {
        if (sum(pop$n)==0) { print(k); break }
        pop <- generations(pop, step)
        if (plotit) { plotpop(pop$n) }
    }
    return(pop)
}

generations <- function (pop, ngens, occupation) {
        for (k in 1:ngens) { 
            pop <- generation(pop) 
        }
        return(pop)
}

generation <- function (pop) {

        pop$gen <- pop$gen + 1
        if (!is.null(pop$params$fitfn)) {
            pop$params$s <- pop$params$fitfn(pop$gen)
        }

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
                if (nmuts>100) { stop("Too many mutations!") }  # can't deal with too many
                if (nmuts>0) for (ell in 1:nmuts) {
                    newtype <- sample((1:dim(pop$n)[3])[smallAreas],1)  # what type is it
                    newloc <- which(rmultinom(1,1,pop$n[,,k]) > 0)      # where is it
                    pop$n[,,newtype][newloc] <- pop$n[,,newtype][newloc] + 1  # add it in
                    pop$mutations <- c(pop$mutations, list(c(pop$gen,k,newtype,newloc))) # record this event
                }
            }
        }

        # Migration
        #   doesn't allow multiple migrations by accumulating migrants in 'inmigrants'
        #   and adding these back in at the end
        mprobs <- sapply( pop$params$migrsteps, "[[", 1 )  # prob of each migration step
        mprobs <- pop$params$m * mprobs/sum(mprobs)
        for (k in seq_along(mprobs)[-1]) { mprobs[k] <- mprobs[k] / prod(1-mprobs[1:(k-1)]) }  # convert to sequence of binomials
        for (k in 1:dim(pop$n)[3]) {
            inmigrants <- array(0,dim=dim(pop$n)[1:2])
            for (j in seq_along(mprobs)) {
                migrants <- migrate(pop$n[,,k,drop=FALSE], prob=mprobs[j], migr=pop$params$migrsteps[[j]])
                pop$n[,,k] <- pop$n[,,k] - migrants$outmigrants 
                inmigrants <- inmigrants + migrants$inmigrants
            }
            pop$n[,,k] <- pop$n[,,k] + inmigrants
        }

        # Resample down to N
        #   (approximate hypergeometric by multinomial)
        pop$n <- aperm( apply(pop$n, c(1,2), function (nn) { rmultinom(1, pop$params$N, nn) } ), c(2,3,1))

        return(pop)
}

migrate <- function ( n, m, prob, migr ) {
    # does a single migration step on a matrix n.
    # with total probability per individual of migration m.
    # expects migr to be a triplet of the format (probability, dx, dy)

    if (missing(prob)) { prob <- migr[1]*m }
    dx <- migr[2]
    if (length(migr)==2) { dy <- 0 } else { dy <- migr[3] }
    if ( dim(n)[3] != 1 ) { 
        warn( "Passing multiple types to migrate!" ) 
    } else { dim(n) <- dim(n)[1:2] }

    outmigrants <- rrbinom( n, prob )  # where they came from (need to subtract this off)
    migrants <- matrix( 0, nrow=nrow(outmigrants), ncol=ncol(outmigrants) )

    nc <- dim(n)[2] # x-extent of populations
    nr <- dim(n)[1] # y-extent of populations
    fromcols <-  1:(nc-abs(dx)) + if (dx<0) { abs(dx) } else { 0 }
    tocols <-  1:(nc-abs(dx)) + if (dx<0) { 0 } else { abs(dx) }
    fromrows <-  1:(nr-abs(dy)) + if (dy<0) { abs(dy) } else { 0 }
    torows <-  1:(nr-abs(dy)) + if (dy<0) { 0 } else { abs(dy) }
    migrants[torows,tocols] <- outmigrants[fromrows,fromcols]

    return( list( outmigrants=outmigrants, inmigrants=migrants ) )
}


#####
# Find equilibrium freqs

meanmigrate <- function ( n, m, migr ) {
    # does as migrate, but in expectation.
    prob <- migr[1]*m
    dx <- migr[2]
    if (length(migr)==2) { dy <- 0 } else { dy <- migr[3] }
    if ( dim(n)[3] != 1 ) { 
        warn( "Passing multiple types to migrate!" ) 
    } else { dim(n) <- dim(n)[1:2] }
    migrants <- n*prob
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
    return( n + migrants - outmigrants )
}

getequilibrium <- function (params,nreps=1000,keep.convergence=FALSE, init=0.5 ) {
    # find equilibrium solution
    # matches (absorbing) boundary conditions of the simulation above
    if ('pop'%in%names(params)) { params <- params$pop$params }
    if ('params'%in%names(params)) { params <- params$params }
    s <- params$s
    migrsteps <- params$migrsteps
    m <- params$m
    n <- array( c(1-init,init), dim=c(dim(s),2) )
    f <- function (n) {
        n[,,1] <- n[,,1] * (1 + params$r)
        n[,,2] <- n[,,2] * (1 + (1+s) * params$r)  #  see above
        for (k in 1:2) {
            for (migr in migrsteps) {
                n[,,k] <- meanmigrate(n[,,k,drop=FALSE], m, migr)
            }
        }
        n <- sweep( n, c(1,2), rowSums(n,dim=2), "/" )
        return(n)
    }
    if (keep.convergence) {
        xx <- array( 0, dim=c(dim(s),nreps) )
        xx[,,1] <- n[,,2]
        for (k in 2:dim(xx)[3]) {
            n <- f(n)
            xx[,,k] <- n[,,2]
        }
        # matplot(xx[1,,floor(seq(1,nreps,length.out=100))],type='l',col=rainbow(1.2*100))
    } else {
        for (k in 2:nreps) {
            n <- f(n)
        }
        xx <- matrix( init, nrow=nrow(s), ncol=ncol(s) )
        xx[] <- n[,,2]
    }
    return (xx)
}


####
# Finding adaptation times in sims

adapttime.1D <- function (pophist) {
    # Find the patches
    s <- as.vector(pophist$pop$params$s)
    patchnum <- (s>0)*c(0,cumsum(diff(s>0)>0))  # works in 1D
    patchsize <- table(patchnum)
    # mean frequency in each patch in leat 80% of runs
    finalfreqs <- tapply( rowMeans( pophist$pophist[,,2,seq(.8*dim(pophist$pophist)[4],dim(pophist$pophist)[4])] ), patchnum, mean )/pophist$pop$params$N
    # trajectory of mean frequencies in each patch
    meanfreqs <- apply( pophist$pophist[,,2,] , 2, tapply, patchnum, mean )/pophist$pop$params$N
    # these patches seem to have adapted
    adapted <- ( finalfreqs > 2*finalfreqs[1] )
    # when is "adaptation"? last time passes max for unadapted patch
    thresh <- max(meanfreqs[1,])
    adapttimes <- apply( meanfreqs, 1, function (x) { length(x)-which.min(rev(x>thresh)) } )
    return( data.frame( patch=seq_along(finalfreqs)-1, final=finalfreqs, adapted=adapted, time=adapttimes, size=as.vector(patchsize) ) )
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


plotpophist <- function(pophist,timeslice=1:dim(pophist$pophist)[4],transposed=FALSE,...) {
    # Plot one-dimensional population history with time on the x-axis
    if (!any(dim(pophist$pophist)[1:2]==1)) { stop("Population is not one-dimensional.") } 
    if (!transposed) {
        popn <- aperm(pophist$pophist[,,,,drop=TRUE],c(3,1,2))
        zimg <- plotpop( popn, params=pophist$pop$params, xlab="time (generations)", ylab="space (demes)", userows=timeslice, ... )
    } else {
        popn <- aperm(pophist$pophist[,,,,drop=TRUE],c(1,3,2))
        zimg <- plotpop( popn, params=pophist$pop$params, xlab="space (demes)", ylab="time (generations)", usecols=timeslice, ... )
    }
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

plotpop <- function (n, params, maxtimes=200, thresh=.05, cols, coltrans, userows=1:dim(n)[1], usecols=1:dim(n)[2], ...) {
    # Layer the colors on top of eachother (multiplicatively) 
    # for the different types

    if (missing(params)) { params <- n$params }
    if ("n" %in% names(n)) { n <- n$n }
    if (length(dim(n))==4 & dim(n)[4]==1) { dim(n) <- dim(n)[1:3] }
    if (is.logical(userows)) { userows <- which(userows) }
    if (is.logical(usecols)) { usecols <- which(usecols) }

    # subsample if too many columns
    if (length(userows)>maxtimes) { 
        userows <- userows[unique(seq(1,length(userows),length.out=maxtimes))]
    }
    if (length(usecols)>maxtimes) { 
        usecols <- usecols[unique(seq(1,length(usecols),length.out=maxtimes))]
    }
    n <- n[userows,usecols,] 
    # color picking
    if (missing(cols)) { cols <- c(list(c("#FFFFFF")), lapply(rainbow(dim(n)[3]), colorRampTrans, n=12) ) }
    if (missing(coltrans)) { coltrans <- function (x,colors,ncols=length(colors)) { colors[ceiling((ncols-1)*x/params$N+1)] } }
    # merge colors
    img <- array(1,dim=c(3,dim(n)[1:2]))  # has r,g,b layers above population
    for (k in 2:dim(n)[3]) {
        # image(1-log(1+n[,,k])/log(N),zlim=c(0,1), col=cols[[k]], new=k>1)
        # image(n[,,k]/N,zlim=c(0,1), col=cols[[k]], new=k>1)
        newcolors <- sapply(coltrans(n[,,k],cols[[k]]),col2rgb)/255
        dim(newcolors) <- dim(img)
        img <- img * newcolors
    }
    img <- rgb(img[1,,],img[2,,],img[3,,])
    dim(img) <- dim(n)[1:2]
    plotpopcircles(img,dim(n),colors=sapply(cols,function(x)x[1]), rowtimes=userows, coltimes=usecols, ...)

    # Add points for transient populations
    #for (k in 2:dim(n)[3]) {
    #    if ( any(n[,,k] > 0) & all(n[,,k] < N*thresh) ) {
    #        mutlocs <- which(n[,,k] > 0)
    #        points( row(n[,,k])[mutlocs], col(n[,,k])[mutlocs], col=cols[[k]][floor(length(cols[[k]])/2)], pch="." )
    #    }
    #}

    return(invisible(img))
}

plotpopcircles <- function(img,dims=dim(img),colors,plotlegend=TRUE,rowtimes=1:nrow(img),coltimes=1:ncol(img),...) {
    # Do the plotting.
    cex <- 3.5*min( par("pin") / (par("cin")*dims[1:2]) )
    plot(rowtimes[row(img)],coltimes[col(img)],col=img,pch=20,cex=cex,xaxt="n",yaxt="n",...)
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
