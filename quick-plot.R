#!/usr/bin/Rscript
#
# make images of simulation runs
#  (works for 1D ones atm)

fnames <- commandArgs(TRUE)
onames <- gsub("[.][^.]*$",".png",fnames)

source("sim-patchy-selection-fns.R")

for (k in seq_along(fnames)) {
    fname <- fnames[k]
    oname <- onames[k]
    load(fname)
    png(file=oname,width=12*144,height=6*144,res=144,pointsize=10)
    layout(t(1:2))
    plotpophist(pophist)
    freqs <- pophist$pophist[,,2,]/pophist$pop$params$N
    matplot(freqs[,floor(seq(1,ncol(freqs),length.out=25))],type='l', col=heat.colors(25), ylim=c(0,1))
    ss <- as.vector(pophist$pop$params$s)
    lines( (ss-min(ss))/diff(range(ss)), col='black', lty=2, lwd=2 )
    dev.off()
}
