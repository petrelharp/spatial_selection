#######
# Find prob of survival with nonconstant selection

iterroot <- function (genfn,tol=sqrt(.Machine$double.eps),maxiter=1000,...) {
    p0 <- genfn(1/2,...)
    for (k in 1:maxiter) {
        p <- genfn(p0,...)
        if (max(abs(p0-p))<tol) { break } else { p0 <- p }
    }
    attr(p,"niter") <- k
    attr(p,"converged") <- (k<maxiter)
    return(p)
}

dx <- 15            # distance between demes (km)
rangesize <- 1000   # half-length of range (km)
xx <- seq(-rangesize,rangesize,by=dx)  # deme locations (km)
sb <- .1
sd <- .00001
elevends <- c(500,2500)     # endpoints of elevation (m)
clinewidth <- 500           # a much wider cline
clinewidth <- 26            # altitudinal cline width (km): S.Am
clinewidth <- 62            # altitudinal cline width (km): Mexico
alpha <- (sb+sd)/clinewidth      # slope of selection with altitude (units of s per km)
selfn <- function (x) { pmin( sb, pmax( -sd, alpha*x ) ) }
selx <- selfn(xx)

# binary branching
m <- .5
eqvals <- c( 1, iterroot( function (p) (1/2)*(1-sb) + (1/2)*(1+sb)*p^2 ) )
genfn <- function (f,selx,migr=m) { (1/2)*(1-selx) + (1/2)*(1+selx)*f*( f + migr*diff(c(eqvals[1],f,eqvals[2]),differences=2) ) }
# Poisson
m <- .5
eqvals <- c( 1, iterroot( function (p) exp((1+sb)*(p-1)) ) )
genfn <- function (f,selx,migr=m) { exp((1+selx)*(f+migr*diff(c(eqvals[1],f,eqvals[2]),differences=2)-1)) }

# maize
N <- 103359  # plants per field
Nf <- 561    # ears to make next generation
pe <- 0.068  # prob of total seed stock replacement
mg <- .0083  # pollen migration
psig <- 0    # pollen migration distance
ssig <- 50   # seed stock "migration" distance
pm <- .02    # prob of partial stock addition
rm <- 0.2     # proportion of stock replaced
poisgenfn <- function (f,lambda) { 
    ans <- exp(lambda*(f-1)) 
    if (any(!is.finite(ans))) { stop("oops! nonfinite??") } else { return(ans) }
}
migrave <- function (f,pmig) {
    # average of f[x+X] over X chosen from [-sig,sig] with pmig giving the probs of moving
    #  (zero distance is what's left over)
    pzero <- 1-sum(pmig)
    sig <- length(pmig)
    ix <- outer( seq_along(f), (-sig:sig), "+" ) - .5
    ix <- .5 + ( ( (1-2*((ix%/%length(f))%%2))*ix ) %% length(f) )
    fx <- f[ix]
    dim(fx) <- dim(ix)
    return( rowSums(fx*c(rev(pmig/2),pzero,pmig/2)[col(fx)]) )
}
genfn <- function (f,selx,migr=1) {
    s <- (1+2*selx)/2
    polm <- 0
    seedm <- (1-1/(1+(ssig/dx)))^(0:(6*ssig/dx)); seedm <- migr*(seedm/sum(seedm))[-1]  # geometric with mean ssig
    return( poisgenfn( migrave(f,polm), mg*s ) * # pollen migration
        ( pe +  # local extinction
        ( (1-pe)*(1-pm) * poisgenfn(f,(1-mg)*s) * # local pollen
            poisgenfn(poisgenfn(f,(N/Nf)*s),Nf/N) + # local seeds
          (1-pe)*pm * poisgenfn(f,(1-rm)*(1-mg)*s) *       # local pollen, partial new stock
            poisgenfn(poisgenfn(f,(N/Nf)*s),(1-rm)*Nf/N)  # local seeds, partial new stock
        ) *
        ( (pe/(1-pe)) * poisgenfn(migrave(f,seedm),(1-mg)*s) * # replant other field, pollen
            poisgenfn(poisgenfn(migrave(f,seedm),(Nf/N)*s),Nf/N) + # replant other field, seeds
          pm * poisgenfn(migrave(f,seedm),rm*(1-mg)*s) *  # partial replant other field, pollen
            poisgenfn(poisgenfn(migrave(f,seedm),(N/Nf)*s),rm*Nf/N) + # partial replant other field, seeds
          (1-pe/(1-pe)-pm) ) # no replanting elsewhere
        ) )
}

eqvals <- c( 1, iterroot( function (p) genfn(p,selx=sb,migr=0) ) )
diffx <- 1-c(.00001,0)
tmpsl <- diff(genfn(diffx,selx=sb,migr=0))/diff(diffx)
meanoff <- diff(genfn(1-c(.00001,0),selx=sb,migr=0))/diff(1-c(.00001,0))
varoff <- (1/2)*diff(genfn(1-c(.00002,.00001,0),selx=sb,migr=0),differences=2)/prod(diff(1-c(.00002,.00001,0))) - meanoff^2

# calculate mean of n log n


niter <- 1000

f <- function (x) { 1/2 + sb*tanh(-x/30) }
ff <- matrix(NA,ncol=niter,nrow=length(xx))
ff[,1] <- f(xx)
for (k in 2:niter) { ff[,k] <- genfn(ff[,k-1],selx) }

layout(t(1:2))
matplot(xx,ff[,floor(seq(2,ncol(ff),length.out=100))],type='l',ylim=1-2*rev(1-eqvals),col=rainbow(120)[1:100],xlab="distance (km)", ylab="prob of exinction")
lines(xx,ff[,1],col='red')
abline( h=eqvals[2], lty=2, lwd=2 )
text( min(xx), eqvals[2], labels="equilibrium value", adj=c(0,0) )
abline(v=0,col='green')
# par(new=TRUE)
# plot(xx,selfn(xx),col='green',type='l',yaxt='n',ylab='')
# axis(4)
plot(xx,1-ff[,niter],type='l',lwd=2,ylim=c(0,max(1-ff[,niter],2*sb/varoff)),xlab="distance (km)", ylab="prob of survival")
abline(h=1-eqvals[2], lty=2, lwd=2, col='red')
text( min(xx), 1-eqvals[2], labels="equilibrium value", adj=c(0,0) )
abline(v=c(0,xx[min(which(selx>-sd))],xx[max(which(selx<sb))]),col='green')
lines(xx,1-iterroot( genfn, selx=selx ), col='red', lty=2 )
lines(xx,2*selx/varoff)


pdf(paste("prob-estab-",clinewidth,".pdf",sep=''), width=3, height=3, pointsize=10)
par(mar=c(5,4,1,1)+.1)
plot(xx,1-ff[,niter],type='l',lwd=2,ylim=c(0,max(1-ff[,niter],2*sb/varoff)),xlab="distance (km)", ylab="prob of survival")
abline(h=1-eqvals[2], lty=2, lwd=2, col='red')
abline(v=c(0,xx[min(which(selx>-sd))],xx[max(which(selx<sb))]),col='green')
lines(xx,1-iterroot( genfn, selx=selx ), col='red', lty=2 )
lines(xx,2*selx/varoff)
legend('topleft',legend=c('truth','2*s/var','cline location'),lty=c(1,1,1),col=c('black','black','green'),lwd=c(2,1,1),bg='white')
dev.off()

## Look at the fixed-s generating functions
layout(t(1:2))
tmpx <- seq(0,1,length.out=1000)
diffx <- 1-c(.00001,0)
plot( tmpx, genfn(tmpx,selx=0,migr=0), ylim=c(0,1), type='l', col='blue' ) # deleterious
tmpsl <- diff(genfn(diffx,selx=0,migr=0))/diff(diffx); abline( 1-tmpsl, tmpsl, col='blue', lty=3 )
lines( tmpx, genfn(tmpx,selx=-sd,migr=0), col='red' ) # neutral
tmpsl <- diff(genfn(diffx,selx=-sd,migr=0))/diff(diffx); abline( 1-tmpsl, tmpsl, col='red', lty=3 )
lines( tmpx, genfn(tmpx,selx=sb,migr=0), col='green' )  # beneficial
tmpsl <- diff(genfn(diffx,selx=sb,migr=0))/diff(diffx); abline( 1-tmpsl, tmpsl, col='green', lty=3 )
tmpvar <- (1/2)*diff(genfn(1-c(.00002,.00001,0),selx=sb,migr=0),differences=2)/prod(diff(1-c(.00002,.00001,0))) - tmpsl^2
abline(0,1,lty=2)
legend("bottomright",lty=1,col=c('green','blue','red'),legend=paste("s=",c(sb,0,-sd)))
tmpx <- seq(1-2*(1-eqvals[2]),1,length.out=1000)
plot( tmpx, genfn(tmpx,selx=sb,migr=0)-tmpx, type='l', col='green' )
lines( tmpx, genfn(tmpx,selx=0,migr=0)-tmpx, col='blue' )
lines( tmpx, genfn(tmpx,selx=-sd,migr=0)-tmpx, col='red' )
points( c(eqvals[2],1-2*(tmpsl-1)/tmpvar), c(0,0), pch="*", cex=2, col=c("black","red") )
legend("topright",pch="*",col=c("black","red"),legend=c("true","approx"))
abline(h=0,lty=2)

########
## old stuff
if (FALSE) {


# From Jeff:
N <- 103359
Nf <- 561
m <- 0.2
e <- 0.068
mg <- .0083
pm <- .02
n <- 325

# from their paper
N <- 100000
Nf <- 500
m <- 0.3
e <- 0.1
mg <- .018
pm <- .02
# number of demes = ???
n <- 325  


# calculated
Ne <- 4 * Nf * N / (3*Nf + N)
Nfm <- Nf*m

a1 <- (1/4) * pm * (1-e) * (1-m)^2
a2 <- (1/4) * pm * (1-e) * m^2
a3 <- (1/4) * (1- (1-e) * pm)
a4 <- ( 3/4 - mg*(1-mg/4) ) * (1 - 2 * pm * m * (1-e) * (1-m)) + (mg * (2-mg) * (1-e) * pm * m * (1-m) + (1/4) * mg^2 )/ (n-1)
a <- c(a1,a2,a3,a4)

b4 <- ( 3/4 - mg*(1-mg/4) ) * (1 - (1-e)^2 * (1-pm*m)^2 ) / ( n *(1-e) - 1 ) + mg*(1-mg/4) / (n-1)
b5 <- (1/4)* ( 1 - (1-e)^2 * (1-pm*m)^2 ) / ( n *(1-e) - 1 )
b <- c(b4,b5)
Q <- b5/(b4+b5)

P1 <- 1/(2*(Nf-Nfm))
P2 <- 1/(2*Nfm)
P3 <- 1/(2*Nf)
P4 <- 1/(2*N)
P <- c(P1,P2,P3,P4)

T0 <- ( (1 - sum(a))/sum(b) + 1 ) / ( P4 * Q * (1-sum(a)) + sum(a*P) )
T1 <- T0 * (1-Q*P4) + 1/sum(b)

Fst <- (1-1/n) / ( (1-1/n) + sum(b) * T0 )

## Results:
# Jeff's parameters:
# T0 = 134640
# T1 = 136724
# Fst = 0.01519433

# paper's parameters:
# T0 = 140435
# T1 = 141809
# Fst = 0.009655918


###
# Branching process
# probability $$(1/2)$$: mean (1+s) (pollen)
# probability $$(1/2) (1-N_f/N) = .4975$$: 0 (not seed corn)
# probability $$(1/2) N_f/N = .0025$$: $$(1+s)N/Nf=200(1+s)$$ (seed corn)
#
# note mean is s
#  and variance is (1/2)*(1+s)( 1 + (1+s)) + (1/2)*N_f/N* (1+s) (N/Nf) (1 + (1+s) N/Nf ) - s^2
#                 = (1/2)(1+s)(2+s) - s^2 + (1/2) (1+s) (1+(1+s)(N/Nf)) , i.e. approximately N/2Nf = 100.

s <- .05
# mixture of Poissons:
genfn <- function (x) { (1/2) * ( exp((1+s)*(x-1)) + (1-Nf/N) + Nf/N * exp((x-1)*(1+s)*N/Nf) ) }
# Find the root of genfn(x)-x=0 :
# solution is approx 0.9999 , i.e. 1-2*s/100.


}
