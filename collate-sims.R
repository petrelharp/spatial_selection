#!/usr/bin/Rscript

if (!interactive()) { basedir <- commandArgs(TRUE)[1] }

source("sim-patchy-selection-fns.R")

raw.simfiles <- list.files( basedir, recursive=TRUE, pattern=".*Rdata" )
list.simfiles <- strsplit(raw.simfiles,"/")
nvars <- max( sapply( list.simfiles, length ) )

if (FALSE) { # for cleanup 
    movethese <- raw.simfiles[sapply(list.simfiles,length)==4]
    for (x in movethese) { file.rename( file.path( "migration", x), file.path( "mutation", x ) ) }
}

# df.simfiles <- tapply( list.simfiles, sapply( list.simfiles, length ), function (x) do.call( rbind, x ) )
# simfiles <- data.frame( rbind( df.simfiles[[2]], cbind( df.simfiles[[1]][,1], NA, df.simfiles[[1]][,-1] ) ) )
simfiles <- as.data.frame( do.call( rbind, list.simfiles ) )
if (ncol(simfiles)==5) {
    names(simfiles) <- c("sm","R","N","dims","filename")
} else {
    names(simfiles) <- c("sm","N","dims","filename")
}

sim.params <- lapply( file.path(basedir,raw.simfiles), function (x) {
            load(x)
            atime.vals <- rep(NA, 4*3)
            names(atime.vals) <- c( paste("final",0:2,sep=''), paste("time",0:2,sep=''), paste("size",0:2,sep=''), paste("hit100.",0:2,sep='') )
            if (min(pophist$pop$params$range)==1) {
                atimes <- adapttime.1D(pophist)
                atime.vals[1:nrow(atimes)] <- atimes$final
                atime.vals[3+(1:nrow(atimes))] <- atimes$time
                atime.vals[2*3+(1:nrow(atimes))] <- atimes$size
                atime.vals[3*3+(1:nrow(atimes))] <- atimes$hit100
            }
            c( pophist$pop$params[c("mu","r","m","N","range","ntypes","sb","sm","nsteps","stepsize","sigma")], atime.vals )
        } )

df.params <- do.call(rbind, lapply(sim.params,unlist))

all.params <- cbind(simfiles,df.params)

write.table(all.params, file=paste(basedir,"-sims-info.tsv",sep=''), sep='\t', row.names=FALSE, quote=FALSE)


