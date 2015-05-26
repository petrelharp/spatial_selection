ninds <- 1e3
ngens <- 6e2

# reviewer's assertion: families of size K live for K generations and so contribute K^2 to the total
# my assertion:
#  E[ sum_{t \ge 0} Z_t | Z_0=K ] = K sum_{t \ge 0} (1-s)^t = K / s 
# oh: but time TO HIT K starting from 1 is proportional to K, given that you do
#      and total contribution in this case is indeed K^2

f <- function (ninit,s) {
    gens <- matrix( 0, nrow=ninds, ncol=ngens )
    gens[,1] <- ninit
    for (k in 2:ncol(gens)) {
        gens[,k] <- rpois( nrow(gens), gens[,k-1]*(1-s) )
    }

    # this does total time
    if (FALSE) {
        t.zero <- 1+rowSums(gens!=0)
        hist(t.zero,main=paste("Z0=", ninit, "s=",s),xlim=range(t.zero,ninit))
        abline(v=ninit, col='red',lwd=2)
        abline(v=mean(t.zero),lwd=2,col='green')
    }

    # layout(t(1:2))
    # plot(colMeans(ttt,na.rm=TRUE))
    # matplot(t(gens[apply(gens,1,max)>10,1:100]),type='l')

    # and this does total contribution
    total.contrib <- rowSums(gens)
    hist(total.contrib,main=paste("Z0=", ninit, "s=",s),xlim=range(total.contrib,ninit^2))
    abline(v=ninit^2, col='red',lwd=2)
    abline(v=mean(total.contrib),lwd=2,col='green')
    abline(v=ninit/s,lwd=2,col='purple',lty=2)
}


layout(matrix(1:6,nrow=2))

f(20,.02)
f(20,.2)
f(40,.02)
f(40,.2)
f(80,.02)
f(80,.2)
