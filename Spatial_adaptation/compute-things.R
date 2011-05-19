source("Spatial_adaptation/standing-variation-fns.R")
require(xtable)

examples <- list( 
    "HbS0" = list(
        mu = 10^(-8),
        rho = 2.5,
        xi = 2,
        # sb = .05,  # this is what we had before?
        sb = .15,  # but it seems to actually be this (currat et al 2002, hedrick review)
        sd = .1,  # ???
        sigma = 10,  # what we had before
        R = 4000
    ),
    "HbS weak" = list(
        mu = 10^(-8),
        rho = 2.5,
        xi = 2,
        sb = .15,
        sd = .1,  # ???
        sigma = 50,  # closer to estimated value of 100?  gets right char length
        R = 4000
    ),
    "HbS strong" = list(
        mu = 10^(-8),
        rho = 2.5,
        xi = 2,
        sb = .15,  
        sd = .5,  # doesn't matter much
        sigma = 50,
        R = 4000
    ),
    "HbC" = list(
        mu = 10^(-8),
        rho = 2.5,
        xi = 2,
        sb = .07,  # Wood et al 2005 (.04-.09)
        sd = .1,  # ???
        sigma = 50,
        R = 4000
    ),
    "G6PD strong" = list(
        mu = 150*10^(-8),  # 150 bp
        rho = 2.5*0.75,  # on the X
        xi = 2.5,
        sb = .25,  # Slatkin et al 2008
        sd = .1,
        sigma = 50,
        R = 4000
    ),
    "G6PD weak" = list(
        mu = 150*10^(-8),  # 150 bp
        rho = 2.5*0.75,  # on the X
        xi = 2,
        sb = .04,  # Tiskkoff et al 2001
        sd = .1,
        sigma = 50,
        R = 4000
    )
)

examples <- lapply(examples, function (x) {
    with(x, {
            x$z0 <- standingProportionArea(mu,rho,sb,sd,sigma)$value
            x$chi <- charLength(mu,rho,sb,sd,sigma)$value
            x$Etau <- meanTime(mu,rho,sb,sd,sigma)$value
            x$mixing <- (R/sigma)^2  # time until pattern mixes
            x
        } ) } )

dx <- length(examples[[1]])
examples.mat <- unlist(examples)*1.0  # make sure they're all in the same numeric mode for xtable
dim(examples.mat) <- c(dx,length(examples.mat)/dx)
dimnames(examples.mat) <- list(names(examples[[1]]),names(examples))
examples.mat <- t(examples.mat)
extable <- xtable(examples.mat,digits=0)
digits(extable)[1+which(apply(examples.mat, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 2
digits(extable)[1+which(apply(examples.mat, 2, function(x) { any(abs(x)<.01) }))] <- -2

filename <- "examples.tex"
write("\\documentclass{article} \\begin{document}", file=filename)
print(extable, file=filename, append=TRUE)
write("\\end{document}", file=filename, append=TRUE)

