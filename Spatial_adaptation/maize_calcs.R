
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


