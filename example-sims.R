#!/usr/bin/Rscript
#
# make a pdf with example simulations from various parameter values

if (!interactive()) { basedir <- commandArgs(TRUE)[1] }

source("sim-patchy-selection-fns.R")

oname <- paste( gsub("/$","",basedir), "-examples", sep='' )

# number of each to choose
nsims <- 2
ntimes <- 300

recurse.dir <- function (dirname,pos=list()) {
    newpos <- c(pos,dirname)
    newpath <- do.call(file.path,newpos)
    cat(newpath,"\n")
    subdirs <- list.dirs(path=newpath,recursive=FALSE,full.names=FALSE)
    subdirs <- subdirs[order( as.numeric(gsub("[^0-9]","",subdirs)) )]
    for (subdir in subdirs) { recurse.dir(subdir,newpos) }
    raw.simfiles <- list.files( newpath, pattern=".*Rdata" )
    for ( fname in raw.simfiles[sample.int(length(raw.simfiles),min(length(raw.simfiles),nsims))] ) {
        load(file.path(newpath,fname))
        if (min(pophist$pop$params$range)==1) {
            cat("  ", fname,"\n")
            timeslice <- seq(1,dim(pophist$pophist)[4],length.out=ntimes)
            plotpophist(pophist, timeslice=timeslice, main=fname)
            freqs <- pophist$pophist[,,2,timeslice]/pophist$pop$params$N
            matplot(freqs[,floor(seq(1,ncol(freqs),length.out=25))],type='l', col=heat.colors(25), ylim=c(0,1), main=paste(pos,collapse=', '))
            ss <- as.vector(pophist$pop$params$s)
            lines( (ss-min(ss))/diff(range(ss)), col='black', lty=2, lwd=2 )
        }
    }
}

# png(file=paste(oname,".png",sep=''), width=8*144, height=10.5*144, res=144, pointsize=10 )
pdf(file=paste(oname,".pdf",sep=''), width=8, height=10.5, pointsize=10 )

layout(matrix(1:8,ncol=2,byrow=TRUE))
recurse.dir(basedir)

dev.off()
