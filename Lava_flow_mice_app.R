char.length= c(3,30)  ## sigma /sqrt(s)


A=100

#pdf(file="~/Dropbox/postdocs/Peter/Spatial_adaptation/Lava_flow_mice_prob_parallel.pdf")
layout((1:2))
par(mar=c(3,3,0.2,1))
R = seq(5,600,length=10000) ##distance in Km

#layout(t(1:2))
mu = 1e-8
prob.parallel <- mu*A/(mu*A + (1/char.length[1])* exp(-R/char.length[1]))  ##assuming sigma=1
plot(R, prob.parallel,type="l",lwd=3, xlab="R, km",ylab="",cex.lab=1.3,cex.axis=1.3,log="x") 
mtext(side=2,line=2,"Prob. parallel mutation",cex=1.3)
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


####haplotype lengths
##assuming sigma=1
length.hap.1<-sqrt(2)*1/(char.length[1])/R +1/R^2
length.hap.2<-sqrt(2)*1/(char.length[2])/R +1/R^2
plot(R,length.hap.1*100,type="l",log="x",ylab="", xlab="",lty=1,lwd=2,cex.lab=1.3,cex.axis=1.3)
mtext(side=2,line=2,"Shared Haplotype Length, cM",cex=1.3)
mtext(side=1,line=2,"Distance between patches, R km",cex=1.3)
lines(R,length.hap.2*100,lty=2,lwd=2)

legend("topright", legend=expression( sigma/sqrt(s[m]) == 3,sigma/sqrt(s[m]) == 30),col=c("grey","grey"),lty=c(1,2),lwd=2)

dev.copy2pdf(file="~/Dropbox/postdocs/Peter/Spatial_adaptation/Lava_flow_mice_prob_parallel.pdf")