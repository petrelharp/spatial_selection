require(colorspace)

# Given a population history, trace lineages back in time
#  and record statistics about their path

lineages <- function ( pophist, nlin, migrsteps, T=dim(pophist)[4], m, linit=NULL, coal=TRUE ) {
    # simulate the movement of nlin selected lineages back through pophist,
    # beginning at locations in linit 
    #   which is a nlin x 3 matrix recording (x,y,t)
    # using the dispersal pattern migrsteps
    #  migrsteps is of the form list( (prob, dx, dy) ) with sum(prob)=1.
    #  and m is the probability of migrating
    #  pophist is indexed by (x, y, type, time)
    #    with types either 0 (nonselected) or 1 (selected)
    # returns:
    #    LL = location of lineages through time (2 x T x nlin)
    if (missing(migrsteps)) { migrsteps <- pophist$pop$params$migrsteps }
    if (missing(m)) { m <- pophist$pop$params$m }
    if ("pophist" %in% names(pophist)) { pophist <- pophist$pophist }

    # NOTE y comes before x
    ny <- dim(pophist)[1]
    nx <- dim(pophist)[2]

    if (is.null(linit)) { # sample initial locations uniformly
        initpop <- pophist[,,2,T]
        dim(initpop) <- dim(pophist)[1:2]
        locs <- sample( nx*ny, nlin, replace=TRUE, prob=initpop )
        linit <- cbind( col(initpop)[locs], row(initpop)[locs], T )
    } else if (ncol(linit)==2) {
        linit <- cbind(linit,T)
    }

    # add a zero step and reweight migrsteps
    mmm <- sum( sapply( migrsteps, "[", 1 ) )
    msteps <- do.call( rbind, c( list( c(1-m, 0, 0) ), lapply(migrsteps, function (x) { c( x[1]*m/mmm, x[-1] ) } ) ) )

    # the Lineage Locations
    #   and the Local Sizes (migration-weighted)
    #   and coalescent events
    T <- linit[,3]
    LL <- array(NA, c(2,max(T),nrow(linit)))  # dimensions are: (x,y) , (time ago) , (lineage number)
    LS <- array(NA, c(1,max(T),nrow(linit)))  # dimensions are: (size) , (time ago) , (lineage number)
    LL[cbind(1:2,rep(T,each=2),rep(1:nrow(linit),each=2))] <- t(linit[,2:1])
    LS[cbind(1,T,1:nrow(linit))] <- pophist[cbind(linit[,2:1],2,T)]  # ?'[': in x[i], 'i' can be a matrix with as many columns as there are dimensions of ‘x’...
    dimnames(LL) <- list( c("y","x"), dimnames(pophist)[[4]][1:max(T)], NULL )  
    dimnames(LS) <- list( NULL, dimnames(pophist)[[4]][1:max(T)], NULL )  
    coalevents <- if (coal) { data.frame( orig=1:nrow(linit), coalto=NA, x=NA, y=NA, t=NA ) } else { NULL }
    for (t in (max(T)-1):1) {
        # pophist[,,2,T-t+1,drop=FALSE] # is the distribution of selected types in previous generation
        active <- which(!is.na(LS[1,t+1,]))
        for (k in seq_along(active)) {
            # problocs[,1] is vector of weights for where to migrate to
            # problocs[,2:3] is the set of locations
            problocs <- msteps[,c(1,3,2)]
            problocs[,-1] <- rep(LL[,t+1,active[k]],each=nrow(problocs)) - problocs[,-1] # minus, since its backwards movement.
            problocs <- problocs[ (pmin(problocs[,2],problocs[,3])>=1) & (problocs[,2]<=ny) & (problocs[,3]<=nx), , drop=FALSE ]
            problocs[,1] <- problocs[,1] * pophist[cbind(problocs[,-1,drop=FALSE],2,t)]
            LS[,t,active[k]] <- sum(problocs[,1])
            if ( LS[,t,active[k]] <= 0 ) { stop("Ooops!  Nowhere to migrate to at t=", t, "; loc=", paste(LL[,t+1,active[k]],collapse=",")) }
            moveto <- problocs[sample(1:nrow(problocs),1,prob=problocs[,1]),-1]
            LL[,t,active[k]] <- moveto
            if (k>1 & coal) {  # coalesce (in lookdown fashion)?
                cancoal <- which( (LL[1,t,active[1:(k-1)]] == moveto[1]) & (LL[2,t,active[1:(k-1)]] == moveto[2] ) )
                coalwith <- sample(c(cancoal,0), 1, prob=c( rep(1,length(cancoal)), pophist[moveto[1],moveto[2],2,t]-length(cancoal) ) )
                # if (pophist[moveto[1],moveto[2],2,t]==1) { print( c(t=t,k=k,ak=active[k],moveto=moveto,cancoal=cancoal,coalwith=coalwith) ) }
                if (coalwith>0) {
                    coalevents[active[k],-1] <- c(active[coalwith],moveto[2],moveto[1],t)
                    LL[,t,active[k]] <- LS[,t,active[k]] <- NA
                }
            }
        }
    }
    return( list( linit=linit, T=T, lineagelocs=LL, localsizes=LS, coalevents=coalevents ) )
}

getfams <- function (eps=1, lins, params, LL=lins$lineagelocs, coalevents=lins$coalevents) {
    # pull out all lineages in lins coalescing before hitting the patch plus eps * sigma
    coal <- subset(coalevents,params$patchdist[as.matrix(coalevents[,c("y","x")])]>params$patchsize+eps*params$sigma)
    fams <- lapply( sort(setdiff(coal$coalto,coal$orig)), function (k) {
            thisLL <- t(LL[,,k]) 
            exittime <- max( c( min(which(!is.na(thisLL[,1]))), which( params$patchdist[thisLL]<=params$patchsize ) ) ) # last exit from patch
            coaltime <- min( coal$t[coal$coalto==k] )
            usetimes <- min(coaltime,exittime):max(which(!is.na(thisLL[,1])))
            return( data.frame( k=k, t=usetimes, thisLL[ usetimes, ] ) )
        } )
    names(fams) <- sort(setdiff(coal$coalto,coal$orig))
    for (j in 1:nrow(coal)) {
        fams[[paste(coal$coalto[j])]] <- rbind( fams[[paste(coal$coalto[j])]], 
                data.frame( k=coal$orig[j], t=which(!is.na(LL[1,,coal$orig[j]])), na.omit(t(LL[,,coal$orig[j]])) )
            )
    }
    return(fams)
}

plotlins <- function( lineages, coalevents, linit, cols, jit=TRUE, add=FALSE, subset=TRUE, dims=sum((apply(lins$lineagelocs,1,max,na.rm=TRUE))>1), do.coalevents=FALSE, do.linit=FALSE, plotit=TRUE, ... ) {
    # plot lineages, on top of either an image or plotpophistory, depending on dimension
    LL <- if ("lineagelocs" %in% names(lineages)) { lineages$lineagelocs } else { lineages }
    if (missing(coalevents) & do.coalevents & "coalevents" %in% names(lineages)) { coalevents <- lineages$coalevents }
    if (missing(linit) & do.linit & "linit" %in% names(lineages)) { linit <- lineages$linit }
    if (missing(cols) & !do.coalevents) { cols <- adjustcolor( rainbow_hcl(dim(LL)[3]), min(1,max(.05,20/dim(LL)[3])) ) }
    if (missing(cols) & do.coalevents) {
        cols <- rep(NA, dim(LL)[3])
        cols[ is.na(coalevents$coalto) ] <- rainbow_hcl(sum(is.na(coalevents$coalto)))
        while( any(is.na(cols)) ) {
            cols[ is.na(cols) ] <- cols[ coalevents$coalto[ is.na(cols) ] ]
        }
    }
    if (jit) {
        LL <- jitter(LL)
        if (do.coalevents) { coalevents[c("x","y")] <- lapply( coalevents[c("x","y")], jitter ) }
    }
    dothese <- (1:dim(LL)[3])[subset]
    cols <- rep(cols,length.out=sum(dothese))
    if (plotit) {
        if (dims==2) {
            if (!add) { plot( 0, 0, type='n', xlim=c(1,max(LL[2,,],na.rm=TRUE)), ylim=c(1,max(LL[1,,],na.rm=TRUE)), xlab='', ylab='', ... ) }
            lapply( dothese, function (k) { lines( LL[2,,k], LL[1,,k], col=cols[k] ) } )
            if (do.linit) { points( linit[dothese,2], linit[dothese,1], pch=20, cex=.5, col=adjustcolor("black",.5) ) }
            if (do.coalevents) { points( coalevents$x[dothese], coalevents$y[dothese], pch="*", cex=2, col=adjustcolor("black",.5) ) }
        } else {
            lapply( dothese, function (k) { lines( LL[2,,k], col=cols[k] ) } )
            if (do.linit) { points( linit[dothese,3], linit[dothese,1], pch=20, cex=.5, col=adjustcolor("black",.5) ) }
            if (do.coalevents) { points( coalevents$t[dothese], coalevents$x[dothese], pch="*", cex=2, col=adjustcolor("black",.5) ) }
        }
    }
    return( invisible( LL[,,dothese] ) )
}

