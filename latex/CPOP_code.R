devtools::install_github("grosed/cpop")
library(cpop)

####Section 3.2
library(changepoint)
set.seed(1)
x=1:400
y=simulate(x,changepoints=0:7*50,change.slope=c(0.5,(-1)^(1:7))/5)
mu=simulate(x,changepoints=0:7*50,change.slope=c(0.5,(-1)^(1:7))/5,sigma=0)

res=cpop(y)
##difference data
dy=diff(y)
out=cpt.mean(dy,method="PELT",penalty="CROPS",pen.value=c(1*log(399),2*2*log(399)))

##with default parameters -- no changes detected -- not sure if this code needs 
##to do in the paper -- just the plot
plot(1:399,dy)
abline(v=1:7*50,lwd=2,col=2) ##true changes
abline(v=cpts.full(out)[3,],col=4) ##cpt.mean with (closest to the correct) number of changes

##trend filtering plot
library(genlasso)
out=trendfilter(y,x,ord=1)
cv=cv.trendfilter(out)
##NOT SURE THE FOLLOWING PLOTS NEED TO BE IN PAPER.
fit1=out$fit[,out$lambda==cv$lambda.min]
sum(abs(diff(diff(fit1)))>0.001) ##17 changes
plot(x,y,pch=".")
abline(v=(1:399)[abs(diff(diff(fit1)))>0.001],col="grey",lty=2)
lines(x,mu,col=4,lwd=2)
lines(x,fit1,col=2,lwd=2)

cps=rep(NA,dim(out$fit)[2])
for(i in 1:dim(out$fit)[2]){
  cps[i]=sum(abs(diff(diff(out$fit[,i])))>0.001)
}
  
plot(x,y,pch=".")
lines(x,mu,col=4,lwd=2)
j=38
fit2=out$fit[,j]
abline(v=(1:399)[abs(diff(diff(fit2)))>0.001],col="grey")
lines(x,fit2,col=2,lwd=2)

###SECTION 4.1
####Example for non-uniform grid locations
#### [FOR HELP FILE ON CPOP: (i) Perhaps use simulate to simulate the data; (ii) also I would explicitly simulate new data for the non-uniform grid]

set.seed(1)
x=(1:50/5)^2
y=simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sd=1)
mu=simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sd=0)

res=cpop(y,x,beta=2*log(length(y)),sd=1)
plot(res)
## lines(x,mu,col=4) ##can we add a line of the truth?

#####SECTION 4.2
###Example with heterogeneous data
set.seed(1)
sd=(1:50)/25
y=simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sd=sd)

###analysis assume constant noise standard deviation
res=cpop(y,x,beta=2*log(length(y)),sd=sqrt(mean(sd^2)))
plot(res)
###analysis with true noise standard deviation
res.true=cpop(y,x,beta=2*log(length(y)),sd=sd)
plot(res.true) 
##ADD LINE OF TRUTH?

####RETURN TO THIS IN CROPS SECTION

##SECTION 4.3
##Choice of grid

###SIMULATION STUDY -- CODE SNIPPETS NOT FOR PAPER
if(F){
  set.seed(1)
n.st=200*2^(0:5)
##scenario 1 -- one changepoint
K=10 ##number of replications
time1=matrix(NA,nrow=K,ncol=length(n.st))
time1g=matrix(NA,nrow=K,ncol=length(n.st))
for(i in 1:length(n.st)){
  n=n.st[i]
  for(k in 1:K){
    y=simulate(1:n,n/2,0.5,1)
    time1[k,i]=(system.time(cpop(y)))[1]
    time1g[k,i]=(system.time(cpop(y,grid=(1:n.st[1])*(n/n.st[1]))))[1]
    }
}

n.st=200*2^(0:5)
##scenario 2 -- linear increasing changepoint
K=10 ##number of replications
time2=matrix(NA,nrow=K,ncol=length(n.st))
time2g=matrix(NA,nrow=K,ncol=length(n.st))
for(i in 1:length(n.st)){
  n=n.st[i]
  for(k in 1:K){
    m=2*n/n.st[1]
    y=simulate(1:n,0:(m-1)*n/m,c(0.05,0.1*(-1)^(1:(m-1))),1)
    time2[k,i]=(system.time(cpop(y)))[1]
    time2g[k,i]=(system.time(cpop(y,grid=(1:n.st[1])*(n/n.st[1]))))[1]
  }
}
save(time1,time2,time1g,time2g,file="/Users/paulfearnhead/Dropbox/Apps/Overleaf/JSS: CPOP/time.Rdata")
}else{
load("/Users/paulfearnhead/Dropbox/Apps/Overleaf/JSS: CPOP/time.Rdata")
}

###EXAMPLE FIGURE FOR PAPER-- Log-log plot of CPU vs n; with error bars
t1=apply(time1,2,mean);v1=apply(time1,2,var)
t2=apply(time2,2,mean);v2=apply(time2,2,var)
t1g=apply(time1g,2,mean);v1g=apply(time1g,2,var)
t2g=apply(time2g,2,mean);v2g=apply(time2g,2,var)
min1=apply(time1,2,min);max1=apply(time1,2,max)
min1g=apply(time1g,2,min);max1g=apply(time1g,2,max)
min2=apply(time2,2,min);max2=apply(time2,2,max)
min2g=apply(time2g,2,min);max2g=apply(time2g,2,max)

plot(n.st,t1,type="l",lwd=2,log="xy",ylab="CPU time (sec)",xlab="n",ylim=c(min(t2g),max(t1)))
lines(n.st,t1g,pch="x",col=2,lwd=2)
lines(n.st,t2,pch="x",col=3,lwd=2)
lines(n.st,t2g,pch="x",col=4,lwd=2)
for(i in 1:length(n.st)){
  lines(c(n.st[i],n.st[i]),c(min1[i],max1[i]),col=1)
  lines(c(n.st[i],n.st[i]),c(min1g[i],max1g[i]),col=2)
  lines(c(n.st[i],n.st[i]),c(min2[i],max2[i]),col=3)
  lines(c(n.st[i],n.st[i]),c(min2g[i],max2g[i]),col=4)
}
lm(log(t1)~log(n.st))
lm(log(t2)~log(n.st))
lm(log(t1g)~log(n.st))
lm(log(t2g)~log(n.st))

#######CODE FOR PAPER
####SHOW ACCURACY OF GRID
set.seed(1)
n=1000;m=9
cps=sample(1:99,size=m,replace=F)*10+5
y=simulate(1:n,c(0,cps),c(0.05,0.1*(-1)^(1:m)),1)

res=cpop(y,beta=2*log(length(y))) 
plot(res)
grid=1:100*10
res.grid=cpop(y,grid=grid,beta=2*log(length(grid)))
plot(res.grid) ##Ideally overlay estimates on the same plot as res.

cps.est=changepoints(res.grid)$location
adaptive.grid=NULL
for(i in 1:length(cps.est)) adaptive.grid=c(adaptive.grid,cps.est[i]+(-4):5)
res.adaptive.grid=cpop(y,grid=adaptive.grid,beta=2*log(length(y)))

##compare changepoint estimates with the full
changepoints(res)$location
changepoints(res.adaptive.grid)$location

##Can we overlay plots -- e.g. estimate of res/res.grid on one plot;
##Here the results for res.adaptive.grid and res are almost identical except that the 
##formere has a change at 685 and the latter at 686.

##2nd advantage is that you can have a more even grid
index=c(1:100,100+18*1:50)
x=index
y.sub=y[index]

res=cpop(y.sub,x,beta=2*log(length(y.sub))) 
res.grid=cpop(y.sub,x,grid=10*1:100,2*log(100))
plot(res)
plot(res.grid)
sort(cps)
changepoints(res)$location
changepoints(res.grid)$location

cost(res)
cost(res.grid)

#####SECTION 4.4: MINIMUM SEGMENT LENGTH
x=2*1:500
mu=simulate(x,c(0,cps),c(0.05,0.1*(-1)^(1:m)),0)
y=mu+rnorm(length(mu))*0.9
y[100*1:5]=y[100*1:5]+5

res=cpop(y,x,beta=2*log(length(y)))
plot(res)
changepoints(res)$location
res.min=cpop(y,x,beta=2*(log(length(y))),minseglen=10)
changepoints(res.min)$location
plot(res.min) ##again overlay on earlier results would be good.
system.time(cpop(y,x))[1]
system.time(cpop(y,x,beta=2*(log(length(y))),minseglen=10))[1]

###SECTION 4.5:CROPS

##NEW VERSION WITH UPDATE TO CPOP.CROPS
##example 1 -- unknown variance -- return to example from Section 4.1

set.seed(1)
x=(1:50/5)^2
y=2*simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sd=1)
mu=simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sd=0)

res=cpop(y,x,beta=2*log(length(y)))
plot(res) 

res.crops=cpop.crops(y,x,beta_min=1.5*log(length(y)),beta_max=100*log(length(y)))
res.crops=unique(res.crops) ###Delete when this is default in CROPS
plot(res.crops) 


##calculate the BIC under a model of unknown variance
models=cpop.crops.models(res.crops)

M=length(models)
BIC=rep(NA,M)
ncps=rep(NA,M)
n=length(y)
for(i in 1:M){
  ##is there a better way to code this (e.g. to extract number of changepoints)
  ncps[i]=length(changepoints(models[[i]])[,1])
  BIC[i]=n*log(mean((residuals(models[[i]]))^2))+2*ncps[i]*log(n)
}
plot(ncps,BIC,type="l") ##not sure if this is needed?
i=which.min(BIC)
plot(models[[i]])

###similarly for heterogeneous example
set.seed(1)
sd=(1:50)/25
y=simulate(x,changepoints=c(10,50),change.slope=c(0.25,-0.25),sd=sd)

res.crops=cpop.crops(y,x,sd=(1:50)/50,beta_min=log(length(y)),beta_max=10*log(length(y)))
res.crops=unique(res.crops) ##CUT WHEN DEFAULT FOR CROPS
plot(res.crops)

####correlated noise
set.seed(1)
n=500
x=1:n;m=10
mu=simulate(x,changepoints=(n/(m+1))*0:m,change.slope=c(0.1,0.2*(-1)^(1:m)),sd=0)
epsilon=rnorm(n+2)
y=mu+(epsilon[1:n]+epsilon[2:(n+1)]+epsilon[3:(n+2)])/sqrt(3)

res=cpop(y,x,beta=2*log(length(y)))
plot(res)
res.crops=cpop.crops(y,x,beta_min=0.5*log(length(y)),beta_max=40*log(length(y)))
res.crops=unique(res.crops)
plot(res.crops)

segs=segmentations(res.crops)
models=cpop.crops.models(res.crops)

#"elbow" plot
plot(segs[,"m"],segs[,"Qm"],lwd=2,type="l",xlab="No. of changepoints",ylab="Unpenalised Cost")
abline(v=m,col=2)

##plot model given by "elbow"
plot(models[[12]])

###Section 5 -- real data 
data("wavenumber_spectra")

par(mfrow=c(1,1))
plot(wavenumber_spectra[,1],wavenumber_spectra[,2],type="l",xlab="wavenumber (cyc/m)",ylim=range(wavenumber_spectra[,-1]),log="xy",ylab="spectra")
for(j in 2:4) lines(wavenumber_spectra[,1],wavenumber_spectra[,j+1],col=j)

x=log(wavenumber_spectra[-(1:3),1])
y=log(wavenumber_spectra[-(1:3),4])

##naive estimator of variance
sig2=mean(diff(diff(y))^2)/6 ##simple estimate of variance 

res=cpop(y,x,grid=seq(from=min(x),to=max(x),length=200),minseglen=0.5,sd=sqrt(sig2),beta=2*log(200))
plot(res) ##can we change plot so that we can label axes such as this?

auto.corr=acf(residuals(res))$acf
lambda=2*sum(auto.corr[1:21])-1
r2=residuals(res)^2
res=cpop(y,x,grid=seq(from=min(x),to=max(x),length=200),sd=sqrt(mean(r2)),beta=2*lambda*log(200))
plot(res)
###estimate of sigma as b*exp(ax)
auto.corr=acf(residuals(res))$acf
lambda=2*sum(auto.corr[1:21])-1
r2=residuals(res)^2

##estimate variance as exp(a+bx)
loglik=function(par){ ##functino for minus twice log-likelihood
  return(length(r2)*(par[1])+par[2]*sum(x)+sum(r2/(exp(par[1]+par[2]*x))) )
}

est.hat=optim(c(0,0),loglik)
sig2=exp(est.hat$par[1]+est.hat$par[2]*x)

res=cpop(y,x,grid=seq(from=min(x),to=max(x),length=200),sd=sqrt(sig2),beta=2*lambda*log(200))

plot(res)

##Perhaps put three plots of res on the same figure 3x1 plot about 1/3 to 1/2 page in size; ideally with 
##axes labelled.

res.crops <- cpop.crops(y,x,grid=seq(from=min(x),to=max(x),length=200),sd=sqrt(sig2),beta_min=2*log(200),beta_max=2*lambda*log(200))

###LIDAR
library(SemiPar)
data("lidar")
x <- lidar[,1]
y <- lidar[,2]

sig2=mean(diff(diff(y))^2)/6 ##simple estimate of variance 
res=cpop(y,x,grid=seq(from=min(x),to=max(x),by=1),minseglen=0.5,sd=sqrt(sig2),beta=2*log(331))
plot(res) ##can we change plot so that we can label axes such as this?
r2=residuals(res)^2

##estimate variance as exp(a+bx)
loglik=function(par){ ##functino for minus twice log-likelihood
  return(length(r2)*(par[1])+par[2]*sum(x)+sum(r2/(exp(par[1]+par[2]*x))) )
}

est.hat=optim(c(0,0),loglik)
sig2=exp(est.hat$par[1]+est.hat$par[2]*x)

res=cpop(y,x,grid=seq(from=min(x),to=max(x),by=1),minseglen=0.5,sd=sqrt(sig2),beta=2*log(331))

res.crops <- cpop.crops(y,x,grid=seq(from=min(x),to=max(x),length=200),sd=sqrt(sig2),beta_min=log(200),beta_max=200*log(200))
res.crops <- cpop.crops(y,x,grid=seq(from=min(x),to=max(x),length=200),sd=sqrt(x)/mean(x),beta_min=log(200),beta_max=200*log(200))

res1=cpop(y,x,grid=c(500,550,600,650),minseglen=0.5,sd=sqrt(sig2),beta=0)


