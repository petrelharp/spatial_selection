source("patchy-selection-fns.R")
require(xtable)

# Helianthus petiolaris/neglectus example:
#  H.p. in southern Colorado; H.n. in Texas-ish

# variance effective population density: total numbers / variance in offspring number
rhovals <- c( 1e4 / 20, 1e5 / 10 )  # per km^2
# dispersal distance: seeds, pollen
sigmavals <- c(.020, .100, 1)  # km / gen
# beneficial selection coefficient within the patches
# --- set so that prob of establishment is 2*sb
sb <- 0.2
# deleterious selection coefficient outside of patches
sm <- .001
# mutation rate (note things will be linear in this
muvals <- c(10^-8,10^-5)
# size of patches is A <- 60 # km^2
# ... but here we want the source area available for new mutations to come from
#    if the dunes are unadapted, i.e. a ring around them of distance, say, 4*sigma
A <- function(sigma) { 4*sigma*2*pi*10 }
# distance between patches
R <- 1000 # km

###
# Table up those values
params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)
exvalues <- t( sapply(1:dim(params)[1], function (k) {
        with( params[k,], everything(mu=mu, rho=rho, sb=sb, sm=sm, sigma=sigma, R=R, A=A(sigma)) )
        } ) )

# write out the table in latex
extable <- xtable(exvalues,digits=0)
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 3
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<.001) | any(abs(x)>1e5) }))] <- -2

filename <- "Spatial_adaptation/helianthus-ex-table.tex"
# write("\\documentclass{article} \\usepackage[landscape]{geometry} \\begin{document}", file=filename)
print(extable, file=filename)#, append=TRUE)
# write("\\end{document}", file=filename, append=TRUE)

####
# And look at parallel adaptation within the big species?
source("Spatial_adaptation/standing-variation-fns.R")

# pre-environment-change deleterious selection coefficient
sd <- .01

params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)

exvalues <- t( sapply(1:dim(params)[1], function (k) {
            with(params[k,], everything(mu=mu, rho=rho, sb=sb, sd=sd, sigma=sigma) )
        } ) )

# write out the table in latex
extable <- xtable(exvalues,digits=0)
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<1 | abs((x-floor(x))/x)>.01) }))] <- 3
digits(extable)[1+which(apply(exvalues, 2, function(x) { any(abs(x)<.001) | any(abs(x)>1e5) }))] <- -2

filename <- "Spatial_adaptation/helianthus-standing-ex-table.tex"
# write("\\documentclass{article} \\usepackage[landscape]{geometry} \\begin{document}", file=filename)
print(extable, file=filename)#, append=TRUE)
# write("\\end{document}", file=filename, append=TRUE)


#####
# # Plot against mutation rate  NOT SO USEFUL.
# muvals <- 10^seq.int(-5,-8,length.out=50)
# params <- expand.grid(rho=rhovals, sigma=sigmavals, mu=muvals)
# exvalues <- as.data.frame( t( sapply(1:dim(params)[1], function (k) {
#         with( params[k,], everything(mu=mu, rho=rho, sb=sb, sm=sm, sigma=sigma, R=R, A=A) )
#         } ) ) )
# 
# # coplot(ratioMutMig ~ mu | sigma*rho, data=exvalues, overlap=0)
# rhocols <- c("black","red")  # rainbow_hcl(length(rhovals))
# names(rhocols) <- rhovals
# legtext <- as.expression( lapply( rhovals, function(x) substitute( rho==myrho, list(myrho=x) ) ) )
# layout(c(1,length(sigmavals)))
# for (sigmaval in sigmavals) {
#     plot( ratioMutMig ~ mu, data=exvalues[exvalues$sigma==sigmaval,], col=rhocols[paste(rho)], main=as.expression(substitute(sigma==zzz, list(zzz=sigmaval))), xlab=expression(mu), ylab="Mutation/Migration" )
#     if (sigmaval==sigmavals[1]) { legend("topleft",col=rhocols,pch=1, legend=legtext) }
# }
# 
