#!/usr/bin/Rscript

simfiles <- list.files(pattern="*-pophistory-run.Rdata")
names(simfiles) <- gsub( "([^-]*)-.*","\\1", simfiles)
filetimes <- file.info(simfiles)$mtime
simfiles <- simfiles[ order(filetimes) ]
filetimes[ order(filetimes) ]

simparams <- lapply( simfiles, function (x) {
        load(x)
        c( pophist$pop$params, list( has.occupation=( !is.null(pophist$occupation) | max(pophist$occupation)==0 ) ) )
    } )


simparam.df <- data.frame( t( sapply( lapply( simparams, "[", c("mu","r","m","N","range","patchsize","sb","sm","sigma",'has.occupation') ), unlist ) ) )
simparam.df$nsteps <- sapply( sapply( simparams, "[", "nsteps" ), function (x) { if (is.null(x)) { NA } else { x[[1]] } } )
simparam.df$stepsize <- sapply( sapply( simparams, "[", "stepsize" ), function (x) { if (is.null(x)) { NA } else { x[[1]] } } )
simparam.df$ngens <- with(simparam.df, nsteps*stepsize)

write.table(simparam.df, file='siminfo.tsv', sep='\t', quote=FALSE)
