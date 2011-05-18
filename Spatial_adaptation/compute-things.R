source("Spatial_adaptation/standing-variation-fns.R")
require(xtable)

examples <- list( 
    "HbS" = list(
        mu = 10^(-8),
        rho = 2,
        xi = 2,
        sb = .05,
        sd = .1,
        sigma = 10,
        R = 4000
    ),
    "G6PD" = list(
        mu = 10^(-8),
        rho = 2,
        xi = 2,
        sb = .05,
        sd = .1,
        sigma = 10,
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
digits(extable)[1+which(apply(examples.mat, 2, function(x) { any(abs(x)<1) }))] <- 2

filename <- "examples.tex"
write("\\documentclass{article} \\begin{document}", file=filename)
print(extable, file=filename, append=TRUE)
write("\\end{document}", file=filename, append=TRUE)

