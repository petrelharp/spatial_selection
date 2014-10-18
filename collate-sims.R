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
simfiles <- do.call( rbind, list.simfiles )

names(simfiles) <- c("sm","R","N","dims","filename")

sim.params <- lapply( file.path(basedir,raw.simfiles), function (x) {
            load(x)
            atimes <- adapttime.1D(pophist)
            atime.vals <- rep(NA, 6)
            names(atime.vals) <- c( paste("final",0:2,sep=''), paste("time",0:2,sep='') )
            atime.vals[1:nrow(atimes)] <- atimes$final
            atime.vals[3+(1:nrow(atimes))] <- atimes$time
            c( pophist$pop$params[c("mu","r","m","N","range","ntypes","sb","sm","nsteps","stepsize","sigma")], atime.vals )
        } )

df.params <- do.call(rbind, lapply(sim.params,unlist))

all.params <- cbind(simfiles,df.params)

write.table(all.params, file=paste(basedir,"-sims-info.tsv",sep=''), sep='\t', row.names=FALSE, quote=FALSE)


