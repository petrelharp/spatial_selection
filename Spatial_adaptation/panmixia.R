panmictic <- function (t,s,N) {
    # logistic growth beginning from 1/N
    1/(1 + (N-1)*exp(-s*t))
}

panmicticinv <- function (x,s,N) {
    # inverse of logistic function above:
    # panmictic( panmicticinv(x,s,N), s, N) = x
    log( N*x*(1-1/N)/(1-x) ) / s
}

quadratic <- function (t,s,N,sigma,rho) {
    # as panmictic until rho*sigma^2/N;
    # as pi*(v*t)^2 after that
    # with v = sigma * sqrt(s)
    t0 <- panmicticinv(rho*sigma^2/N,s,N)
    # needs to be vectorized to be able to integrate
    ans <- panmictic(t,s,N)
    ans[t>t0] <- pmin( 1, (pi*rho*sigma^2*s/N) * (t-t0+1/sqrt(s))^2 )[t>t0]
    return(ans)
}

N <- 1e8
s <- .05
sigma <- 20
rho <- 3
ttt <- 4*(1:100)*log(N)/(100*s)

plot(ttt, sapply(ttt, function(t) panmictic(t,s,N)), type="l", ylim=c(0,1))
lines(ttt, sapply(ttt, function(t) quadratic(t,s,N,sigma,rho)), col="red")

integrate(function (t) {panmictic(t,s,N)-quadratic(t,s,N,sigma,rho)}, 0, Inf)

ss <- sqrt(2)^(-(20:1)); 
rr <- (1:1000)
ssrr <- expand.grid(ss,rr)
ddd <- apply(ssrr,1,function(x) { log( log(x[2])/(x[2]*sqrt(x[1]*2)) ) } )
dim(ddd)<-c(length(ss),length(rr))
# contour(x=ss,y=rr,z=ddd,xlab="s",ylab="r")

# png(file="min-s-versus-r.png",width=600,height=400)
pdf(file="min-s-versus-r.pdf",width=6,height=4)
rr <- 10:1000
plot(rr,8*((log(rr)/rr)^2), type="l", xlab=expression(r==sqrt(A/sigma^2)), ylab=expression(s== 8 * ( log(r)/r )^2), log="xy")
dev.off()

NN <- exp(8:20)
plot(NN, 1/(4*NN*log(NN)), type="l")

# png(file="f-versus-gamma.png",width=600,height=400)
pdf(file="f-versus-gamma.pdf",width=6,height=4)
xx <- (1:999)/1000
plot((1-xx^2)/xx^3, xx, type="l", log="xy", xlab=expression(gamma==(1-x)/x^(3/2)), ylab=expression(x==(chi/chi[0])^2) )
dev.off()

ss <- sqrt(2)^(-(20:6)); 
mu <- 1e-8
# mu <- 1e-6
sigsig <- (10:100)
gg <- apply( expand.grid(ss,sigsig), 1, function (ssig) log( sqrt(ssig[1]) / ( ssig[2] * sqrt(4*pi*mu) * log(ssig[2]^2)^(3/2) ) ) )
dim(gg) <- c(length(ss), length(sigsig))
# png(file="gamma-contour-1e-8.png", width=600, height=400)
pdf(file="gamma-contour-1e-8.pdf", width=6, height=4)
contour(x=ss,y=sigsig,z=gg, xlab="s", ylab=expression(sigma), main=expression(paste(log(gamma),", ", mu==10^-8)))
dev.off()

