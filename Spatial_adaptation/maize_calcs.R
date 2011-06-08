
# estimated
N <- 100000
Nf <- 500
m <- 0.3
e <- 0.1
mg <- .018
pm <- .01
# number of demes = ???
# n <- 10000  # gets Fst = .011
n <- 28  # gets Fst = .008


# calculated
Ne <- 4 * Nf * N / (3*Nf + N)
Nfm <- Nf*m

a1 <- (1/4) * pm * (1-e) * (1-m)^2
a2 <- (1/4) * pm * (1-e) * m^2
a3 <- (1/4) * (1- (1-e) * pm)
a4 <- ( 3/4 - mg*(1-mg/4) ) * (1 - 2 * pm * m * (1-e) * (1-m)) + mg * (2-mg) * (1-e) * pm * m * (1-m) + (1/4) * mg^2 / (n-1)
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
# with n = 10000
# T0 = 3690781
# T1 = 3734148
# Fst = 0.01161358

# with n = 28
# T0 = 13858
# T1 = 13975
# Fst = 0.008052662

