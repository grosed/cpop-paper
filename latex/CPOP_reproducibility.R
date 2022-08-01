## cpop: Detecting changes in piecewise-linear signals
#
## Paul Fearnhead and Dan Grose

##Install package
library(cpop)
library(ggplot2)

#################SECTION 1

##FIGURE 1
##additional packages needed for Figure 1
library(changepoint) 
library(genlasso) 

#simulate data
set.seed(1)
x <- 1:400
y <- simulate(x,changepoints=0:7*50,change.slope=c(0.5,(-1)^(1:7))/5)
mu <- simulate(x,changepoints=0:7*50,change.slope=c(0.5,(-1)^(1:7))/5,sd=0)

###implement two methods -- difference and treat as change-in-mean, and trend-filtering

##(1) difference data and analyse as change-in-mean
dy <- diff(y)
##with default parameters -- no changes detected, so run with CROPS and choose "best" segmentation
out <- cpt.mean(dy,method="PELT",penalty="CROPS",pen.value=c(1*log(399),2*2*log(399)))


##(2) trend filtering
out.tf <- trendfilter(y,x,ord=1)
##show two segmentations -- one based on minimising CV
cv <- cv.trendfilter(out.tf)
fit1 <- out.tf$fit[,out.tf$lambda==cv$lambda.min]
sum(abs(diff(diff(fit1)))>0.001) ##17 changes
##one based on the correct number of changes

#to find this we calculate the number of changes as we vary the penalty in trend-filtering
cps <- rep(NA,dim(out.tf$fit)[2])
for(i in 1:dim(out.tf$fit)[2]){
  cps[i] <- sum(abs(diff(diff(out.tf$fit[,i])))>0.001)
}
##inspection gives that the 38th penalty value gives 7 changes.
fit2 <- out.tf$fit[,38]

par(mfrow=c(2,2))

##PLOT
##1a -- plot data 
plot(x,y,pch=".",cex=2,ylab="Data",xlab="Time")
lines(x,mu,lwd=2,col=4)


##1b--differenced data
plot(1:399,dy,ylab="First differences",xlab="Time",pch=".",cex=2)
abline(v=1:7*50-1,lwd=2,col=4,lty=2) ##true changes
abline(v=cpts.full(out)[3,],col=2,lty=2) ##cpt.mean with (closest to the correct) number of changes

##1c-Trend-filtering with CV
plot(x,y,pch=".",cex=2,ylab="Data",xlab="Time")
abline(v=(1:399)[abs(diff(diff(fit1)))>0.001]+1,col=2,lty=2,lwd=1)
lines(x,mu,col=4,lwd=2)
lines(x,fit1,col=2,lwd=2)

##1d -Trend-filtering with correct number of changes
plot(x,y,pch=".",cex=2,ylab="Data",xlab="Time")
lines(x,mu,col=4,lwd=2)
abline(v=(1:399)[abs(diff(diff(fit2)))>0.001]+1,col=2,lty=2,lwd=1)
lines(x,fit2,col=2,lwd=2)

############SECTION 3

##FIGURE 2
##specify mean function by change-points and change in slope
set.seed(1)
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200
sd <- 0.8
y <- simulate(x, changepoints, change.slope, sd) #simulate data

df <- data.frame("x" = x, "y" = y)
p <- ggplot(data = df, aes(x = x, y = y))
p <- p + geom_point(size=1)
p <- p + geom_vline(xintercept = changepoints, color = "red",  linetype = "dashed")

#for the plot: obtain the true mean function by setting sd=0.
mu <- simulate(x, changepoints, change.slope, sd = 0)
p <- p + geom_line(aes(y = mu), color = "blue")
p <- p + theme_bw()
print(p)

##analysis by cpop
res <- cpop(y, x, sd = 0.8)
summary(res)

##Figure 3
p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

print(p)

##Example code for post-processing

changepoints(res)
estimate(res, x = c(0.1,2.7,51.6))
fitted(res)
head(residuals(res))

##############SECTION 4

##FIGURE 4
##Code as for Figure 3 except the values of x are changed
set.seed(1)
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- (1:200)^(2)/(200) ###x-values now unevenly spaced
sd <- 0.8
y <- simulate(x, changepoints, change.slope, sd) #simulate data

#for the plot: obtain the true mean function by setting sd=0.
mu <- simulate(x, changepoints, change.slope, sd = 0)

##analysis by cpop
res <- cpop(y, x, sd = 0.8)
#summary(res)

p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

print(p)

### FIGURE 5
set.seed(1)
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200 
sd <- x/100
y <- simulate(x, changepoints, change.slope, sd) #simulate data
mu <- simulate(x, changepoints, change.slope, sd = 0)

##analysis by cpop
res <- cpop(y, x, sd = sqrt(mean(sd^2)) )
#summary(res)

res.true <- cpop(y, x, sd=sd)
p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

print(p)

p.true <- plot(res.true)
p.true <- p.true + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.true <- p.true + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

print(p.true)

#####SECTION 4.3

##short set of simulations to calculate time of running CPOP
## as n increases: 200, 400, 800, 1600, 3200, 6400
## two scenarios (i) 1 changepoint; (2) segments of length 100
## compare: grid = x; and grid is evenly spaced set of 200 points.


##This takes about 30 minutes.
  set.seed(1)
  n.st <- 200*2^(0:5)
  ##scenario 1 -- one changepoint
  K <- 10 ##number of replications
  time1 <- matrix(NA,nrow=K,ncol=length(n.st))
  time1g <- matrix(NA,nrow=K,ncol=length(n.st))
  for(i in 1:length(n.st)){
    n <- n.st[i]
    for(k in 1:K){
      y <- simulate(1:n,n/2,0.5,1)
      time1[k,i] <- (system.time(cpop(y)))[1]
      time1g[k,i] <- (system.time(cpop(y,grid=(1:n.st[1])*(n/n.st[1]))))[1]
      cat(".")
    }
    cat("\n")
  }
  
  n.st <- 200*2^(0:5)
  ##scenario 2 -- linear increasing changepoint
  K <- 10 ##number of replications
  time2 <- matrix(NA,nrow=K,ncol=length(n.st))
  time2g <- matrix(NA,nrow=K,ncol=length(n.st))
  for(i in 1:length(n.st)){
    n <- n.st[i]
    for(k in 1:K){
      m <- 2*n/n.st[1]
      y <- simulate(1:n,0:(m-1)*n/m,c(0.05,0.1*(-1)^(1:(m-1))),1)
      time2[k,i] <- (system.time(cpop(y)))[1]
      time2g[k,i] <- (system.time(cpop(y,grid=(1:n.st[1])*(n/n.st[1]))))[1]
      cat(".")
    }
    cat("\n")
  }

###average the times
t1 <- apply(time1,2,mean)
t2 <- apply(time2,2,mean)
t1g <- apply(time1g,2,mean)
t2g <- apply(time2g,2,mean)

##FIGURE 6
plot(c(200,6400),c(0.024,170),type="n",xlab="n",ylab="CPU time (sec)",log="xy",yaxt="n")
axis(2,at=c(0.1,1,10),labels=c(0.1,1,10))
lines(n.st,t1,lwd=2)
lines(n.st,t2,lwd=2,col=2)
lines(n.st,t1g,lwd=2,lty=2)
lines(n.st,t2g,lwd=2,col=2,lty=2)
lines(n.st,t2[1]*n.st^1.7/n.st[1]^1.7,col=2,lty=3)
lines(n.st,t1[1]*n.st^2.5/n.st[1]^2.5,col=1,lty=3)

###Example of running cpop on coarse grid and refining
set.seed(1)
x <- 1:6400
y <- simulate(x, changepoints = 0:31*200, change.slope = c(0.05,0.1*(-1)^(1:(31))), sd = 1)

res.coarse <- cpop(y, x, grid=1:399*16, beta = 2*log(400))
cps <- unlist( changepoints(res.coarse) )

grid <- NULL
for(i in 1:length(cps)){
  grid <- c(grid, cps[i] + (-7):8 )
}

res.fine <- cpop(y, x, grid, beta = 2*log(length(x)))

##now compare to true analysis
res.true <- cpop(y, x, beta = 2*log(length(x)))

cp1=unlist(changepoints(res.fine))
cp2=unlist(changepoints(res.true))
cp2-cp1 ##see error in changepoint locations

cost(res.fine) ##look at cost of fitted changepoints -- two stage small grid
cost(res.true) ##look at cost of fitted changepoints -- use full grid

#####SECTION 4.4
set.seed(1)
##simulate data as for Section 3
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200
mu <- simulate(x, changepoints, change.slope, sd = 0) #mean

##simulate data with t noise
y <- mu + rt(length(x), df = 4)


res <- cpop(y,x, beta=2*log(length(y)), sd = sqrt(4/2)) ##variance of t_d is d/(d-2)

##Now run with a minimum segment length
res.min=cpop(y,x,beta=2*(log(length(y))),minseglen=10, sd = sqrt(4/2))
##compare CPU cost
system.time(cpop(y,x, beta=2*log(length(y)), sd = sqrt(4/2) ) )[1]
system.time(cpop(y,x,beta=2*(log(length(y))),minseglen=10, sd = sqrt(4/2) ) )[1]

##run with too large a minimum segment length
res.min40=cpop(y,x,beta=2*(log(length(y))),minseglen=40, sd = sqrt(4/2))
res.min30=cpop(y,x,beta=2*(log(length(y))),minseglen=30, sd = sqrt(4/2))

###FIGURE 7
p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)

p <- plot(res.min)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)

p <- plot(res.min30)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)

p <- plot(res.min40)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)

####SECTION 4.5
set.seed(1)
##simulate data as for Section 3
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200
mu <- simulate(x, changepoints, change.slope,  0) #mean
y <- simulate(x, changepoints, change.slope, 1.5) #data with sd=1.5

###run crops
res.crops <- cpop.crops(y , x, beta_min= 5 ,beta_max= 50, sd =1)

##calculate the BIC under a model of unknown variance
models <- cpop.crops.models(res.crops)

M <- length(models)
BIC <- rep(NA, M)
ncps <- segmentations( res.crops )[,4]
n <- length(y)
for(i in 1:M){
  BIC[i] <- n*log( mean( ( residuals(models[[i]]) )^2) )+ 2 * ncps[i]  * log(n)
}

###FIGURE 8

plot(res.crops)

i <- which.min(BIC)
p <- plot(models[[i]])
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)



####correlated noise
set.seed(1)
n <- 500
x <- 1:n
mu <- simulate(x,changepoints = 45*0:10,change.slope = c(0.15,0.3*(-1)^(1:10)) ,sd = 0)
epsilon <- rnorm(n+2)
y <- mu + (epsilon[1:n] + epsilon[2:(n+1)] + epsilon[3:(n+2)]) /sqrt(3)


res.crops <- cpop.crops(y,x,beta_min=8,beta_max=200)
segs <- segmentations(res.crops)

plot(res.crops)
#"elbow" plot
plot(segs[,"m"],segs[,"Qm"],lwd=2,type="l",xlab="No. of changepoints",ylab="Unpenalised Cost")
abline(v = 10,col = 2)

##plot model given by "elbow"
models <- cpop.crops.models(res.crops)
p <- plot(models[[5]])
p <- p + geom_vline(xintercept = 45*1:10, color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)

####SECTION 5

data("wavenumber_spectra") ##Load data

##example for one curve -- fitting the decay.
x <- log(wavenumber_spectra[-(1:3),1],base=10)
y <- log(wavenumber_spectra[-(1:3),4],base=10)

grid <- seq(from=min(x),to=max(x),length=200)
##naive estimator of variance
sig2 <- mean( diff( diff(y) )^2 )/6 ##simple estimate of variance 

res <- cpop(y, x, grid, sd=sqrt(sig2), minseg = 0.09, beta=2*log(200))

r2 <- residuals(res)^2
##estimate variance as exp(a+bx)
##function for minus twice log-likelihood
loglik <- function(par){ 
  return(length(r2) * par[1] + par[2] * sum(x) + sum( r2/ (exp(par[1]+par[2]*x) )) )
}

est.hat <- optim( c(0,0) , loglik)
sig2 <- exp(est.hat$par[1] + est.hat$par[2]*x)

res2 <- cpop(y, x, grid , sd=sqrt(sig2), minseg= 0.09, beta=2 * log(200))

######FUNCTION TO AUTOMATE THIS ANALYSIS
wavenumber_est <- function(x,y){
  grid <- seq(from=min(x),to=max(x),length=200)
  ##naive estimator of variance
  sig2 <- mean( diff( diff(y) )^2 )/6 ##simple estimate of variance 
  
  res <- cpop(y, x, grid, sd=sqrt(sig2), minseg = 0.09, beta=2*log(200))
  
  r2 <- residuals(res)^2
  est.hat <- optim( c(0,0) , loglik)
  sig2 <- exp(est.hat$par[1] + est.hat$par[2]*x)
  
  res2 <- cpop(y, x, grid , sd=sqrt(sig2), minseg = 0.09, beta=2 * log(200))
  
  return(res2)
}


##FIGURE 9
##plot of raw data

par(mfrow=c(1,1))
plot(wavenumber_spectra[,1],wavenumber_spectra[,2],type="l",xlab="wavenumber (cyc/m)",ylim=range(wavenumber_spectra[,-1]),log="xy",ylab="spectra")
for(j in 2:4) lines(wavenumber_spectra[,1],wavenumber_spectra[,j+1],col=j)

##initial estimate
plot(res,xlab="log_10 (wavenumber)",ylab="log_10 (spectra)") ##can we change plot so that we can label axes such as this?

## second estimate
plot(res2,xlab="log_10 (wavenumber)",ylab="log_10 (spectra)")

##plot all functions
y.est <- matrix(NA,nrow=4,ncol=length(x))
for(j in 1:4){
  est <- wavenumber_est(x, log(wavenumber_spectra[-(1:3),1+j],base=10))
  y.est[j, ] <- estimate(est,x)[,2]
}

par(mfrow=c(1,1))
plot(wavenumber_spectra[-(1:3),1],wavenumber_spectra[-(1:3),2],type="l",xlab="wavenumber (cyc/m)",ylim=range(wavenumber_spectra[,-1]),log="xy",ylab="spectra")
for(j in 2:4) lines(wavenumber_spectra[-(1:3),1],wavenumber_spectra[-(1:3),j+1],col=j)
for(j in 1:4) lines(10^(x),10^(y.est[j,]),lty=2,lwd=2,col=j)
