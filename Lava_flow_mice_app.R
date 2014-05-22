char.length= c(3,30)


A=100

#pdf(file="~/Dropbox/Talks_posters/coopibarra_patchy_2013/figures/prob_parallel_hoekstra_NEW.pdf")
R = seq(5,600,length=10000) ##distance in Km
#layout(t(1:2))
mu = 1e-8
prob.parallel <- mu*A/(mu*A + (1/char.length[1])* exp(-R/char.length[1]))  ##assuming sigma=1
plot(R, prob.parallel,type="l",lwd=3, xlab="R, km",ylab="Prob. parallel mutation.",cex.lab=1.5,cex.axis=1.5,log="x") 
mu=1e-5
prob.parallel <- mu*A/(mu*A + (1/char.length[1])* exp(-R/char.length[1]))  ##assuming sigma=1

lines(R, prob.parallel,lwd=3,col="red")


R = seq(5,600,length=1000) ##distance in Km
mu = 1e-8
prob.parallel <- mu*A/(mu*A + (1/char.length[2])* exp(-R/char.length[2]))  ##assuming sigma=1
lines(R, prob.parallel,lty=2,lwd=3)

mu=1e-5
prob.parallel <- mu*A/(mu*A + (1/char.length[2])* exp(-R/char.length[2]))  ##assuming sigma=1

lines(R, prob.parallel,lwd=3,col="red",lty=2)

legend("topleft", legend=expression(mu == 10^-8, mu == 10^-5, sigma/sqrt(s[m]) == 3,sigma/sqrt(s[m]) == 30),col=c("black","red","grey","grey"),lty=c(1,1,1,2),lwd=2)
#dev.off()
