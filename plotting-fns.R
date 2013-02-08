# misc plotting fns

logcontour <- function(x,y,z,log="",add=FALSE,...) {
    # contour plot with log scale axes
    logx <- grepl("x",log)
    logy <- grepl("y",log)
    xlims <- range( if (logx) { log10(x) } else { x } )
    ylims <- range( if (logy) { log10(y) } else { y } )
    if (!add) { 
        plot(0,type='n',xlab='',ylab='',xlim=xlims, ylim=ylims, xaxt='n', yaxt='n') 
        if (logx) {
            xat <- unique(floor(pretty(c(ceiling(xlims[1]-diff(ylims)/1e5),floor(xlims[2]+diff(ylims)/1e5)))))
            # xat <- axisTicks(range(x),log=logx)
            axis(1, at=xat, labels=10^xat)
        } else {
            axis(1)
        }
        if (logy) {
            yat <- unique(floor(pretty(c(ceiling(ylims[1]-diff(ylims)/1e5),floor(ylims[2]+diff(ylims)/1e5)))))
            # yat <- axisTicks(ylims,log=logy)
            axis(2, at=yat, labels=10^yat)
        } else {
            axis(2)
        }
    }
    contour( if(logx){log10(x)}else{x}, if(logy){log10(y)}else{y}, z, add=TRUE, ... )
}


