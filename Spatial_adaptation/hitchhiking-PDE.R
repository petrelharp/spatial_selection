require(deSolve)

stationaryPhi <- function (t,y,parms,...) {
    # solving the ODE for the stationary wave front
    #   -v f' = f'' + f(1-f)
    #   so here y = c(f,f')
    yp <- y[2]
    ypp <- - parms$v * y[2] - y[1]*(1-y[1])
    return( list( c(yp,ypp) ) )
}

# intital value is approximately (1+t*exp(t),(1+t)*exp(t)) for large negative t
#  yinit <- c(1+min(tt)*exp(min(tt)),100*(1+min(tt))*exp(min(tt)))
# or phi' \approx (1-phi)/v
tt <- (-2000:2000)/100
eps <- .0001
veps <- 0
yinit <- c(1-eps,-eps/sqrt(2))
yytt <- lsoda(y=yinit, times=tt, func=stationaryPhi, parms=list(v=sqrt(2)+veps) )   # dips below zero.  speed not great enough??
plot(yytt[,1], yytt[,2],type="l") 
lines(yytt[,1],yytt[,3]+1)
abline(h=c(0,1))
vepsvals <- (1:40)/40
vepscols <- rainbow(length(vepsvals))
miny <- rep(NA,length(vepsvals))
for (k in 1:length(vepsvals)) {
    yytt <- lsoda(y=yinit, times=tt, func=stationaryPhi, parms=list(v=sqrt(2)+vepsvals[k]) )
    miny[k] <- min(yytt[,2])
    lines(yytt[,1], yytt[,2],col=vepscols[k])
    lines(yytt[,1],yytt[,3]+1,col=vepscols[k])
}
veps <- min(vepsvals[miny>0]) # like 0.5 !?!

vEst <- sqrt(2)+veps
phiEst <- lsoda(y=yinit, times=tt, func=stationaryPhi, parms=list(v=vEst))

####
# try dynamic version

dynamicPhi <- function(t,y,parms,xx,dx,eps=.0001) {
    # solving the PDE
    # and waiting for stationarity
    #   d/dt y = v y' + y'' + y (1-y)
    #  here y is the vector phi(xx)
    #  and yp is the vector phi'(xx)
    n <- length(xx)
    # dx <- diff(xx)  # length n-1
    dy <- diff(y)/dx  # length n-1
    # f'(x) approx ( (b/a+b) (f(x+a)-f(x))/a + (a/a+b) (f(x)-f(x-b))/b )
    yp <- c(dy[1], (dx[-1]*dy[-(n-1)]+dx[-(n-1)]*dy[-1])/(dx[-1]+dx[-(n-1)]) , dy[n-1])  # length n
    # keep boundary conditions fixed to avoid running off?
    y[1] <- 1-eps
    y[n] <- eps
    yp[1] <- yp[n] <- -eps
    # and f''(x) approx 2 * ( (f(x+a)-f(x))/a - (f(x)-f(x-b))/b ) / (a+b)
    dyp <- diff(yp)/dx
    ypp <- c(dyp[1], (dx[-1]*dyp[-(n-1)]+dx[-(n-1)]*dyp[-1])/(dx[-1]+dx[-(n-1)]) , dyp[n-1])
    z <- parms$v * yp + ypp + y*(1-y)
    # boundaries again?
    z[1] <- z[n] <- 0
    # if (any(abs(z)>10)) { browser() }
    return( list( z ) )
}

xx <- (-1000:1000)/50
dx <- diff(xx)
yinit <- (1-tanh(xx/2))/2
eps <- .0001
yinit <- c( rep(1-eps,floor(length(xx)/3)), sort(runif(floor(length(xx)/3)),decr=TRUE), rep(0+eps,length(xx)-2*floor(length(xx)/3)) )
v <- sqrt(2)+2
yytt <- ode.1D(yinit, (1:1200)/20, dynamicPhi,parms=list(v=v),nspec=1,dimens=length(yinit),xx=xx, dx=dx, maxsteps=50000, atol=0 )

speed <- 0
plot( xx-yytt[1,1]*speed, yytt[1,-1], type="l" )
for (k in 2:(dim(yytt)[1]-1)) { lines(xx-yytt[k,1]*speed,yytt[k,-1],col=rainbow(dim(yytt)[1]+5)[k]) }


### useful for debugging
ggg <- function () {
   with( parent.frame(), {
     ylims <- range(c(y,yp,ypp)); plot(xx,y,ylim=ylims)
     points(xx,yp,col="red")
     points(xx,ypp,col="green")
   } )
}
trace(dynamicPhi,exit=ggg)


### PROBLEMS ABOVE. with boundary conditions?
# Transform with logit.

logitPhi <- function(t,y,parms,xx,dx) {
    # solving the PDE
    # and waiting for stationarity
    #   d/dt y = v y' + y'' + y (1-y)
    #  here y is the vector phi(xx)
    #  and yp is the vector phi'(xx)
    ## Try putting through logit function:
    # if df/dt = v f' + f'' + f(1-f)
    # and y = log(f/(1-f))
    # then dy/dt = v y' + (y')^2 + (y'')/(1+e^{-y}) + 1
    n <- length(xx)
    # dx <- diff(xx)  # length n-1
    dy <- diff(y)/dx  # length n-1
    # f'(x) approx ( (b/a+b) (f(x+a)-f(x))/a + (a/a+b) (f(x)-f(x-b))/b )
    yp <- c(dy[1], (dx[-1]*dy[-(n-1)]+dx[-(n-1)]*dy[-1])/(dx[-1]+dx[-(n-1)]) , dy[n-1])  # length n
    # and f''(x) approx 2 * ( (f(x+a)-f(x))/a - (f(x)-f(x-b))/b ) / (a+b)
    dyp <- diff(yp)/dx
    ypp <- c(dyp[1], (dx[-1]*dyp[-(n-1)]+dx[-(n-1)]*dyp[-1])/(dx[-1]+dx[-(n-1)]) , dyp[n-1])
    z <- parms$v * yp + yp^2 + ypp/(1+exp(-y)) + 1
    # if (any(abs(z)>10)) { browser() }
    return( list( z ) )
}

logit <- function (x) { log(x/(1-x)) }
logist <- function (x) { 1 / (1 + exp(-x)) }

xx <- (-1000:1000)/100
dx <- diff(xx)
finit <- (1-tanh(xx/2))/2
# OR
finit <- sort(runif(length(xx)),decr=TRUE)
yinit <- logit(finit)
v <- 4
yytt <- ode.1D(yinit,(1:100)/20,dynamicPhi,parms=list(v=v),nspec=1,dimens=length(yinit),xx=xx, dx=dx, maxsteps=50000 )

speed <- 0
plot( xx-yytt[1,1]*speed, logist(yytt[1,-1]), type="l" )
for (k in 2:(dim(yytt)[1]-1)) { lines(xx-yytt[k,1]*speed, logist(yytt[k,-1]), col=rainbow(dim(yytt)[1]+5)[k]) }
