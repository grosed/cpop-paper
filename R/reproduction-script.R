# cpop: Detecting Changes in Piecewise-Linear Signals
## Paul Fearnhead and Daniel Grose
## Section 1
### Figure 1
#### simulate data

library(cpop)
set.seed(1)
x <- 1:400
changepoints <- 0:7*50
y <- simchangeslope(x,changepoints=changepoints,change.slope=c(0.5,(-1)^(1:7))/5)
mu <- simchangeslope(x,changepoints=changepoints,change.slope=c(0.5,(-1)^(1:7))/5,sd=0)

#### implement three different methods
#### 1. cpop

res <- cpop(y,sd=1) 

#### 2. difference data and analyse as change in mean

dy <- diff(y)

library(changepoint)
out <- cpt.mean(dy,method="PELT",penalty="CROPS",pen.value=c(1*log(399),2*2*log(399)))

#### 3. trend filtering

library(genlasso)
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

library(ggplot2)

df <- data.frame("x" = x, "y" = y,"dy"=  c(NA,dy))

p.1 <- ggplot(data = df, aes(x = x, y = y))
p.1 <- p.1 + geom_point(alpha=0.4)
p.1 <- p.1 + geom_line(aes(y = mu), color = "blue")
p.1 <- p.1 + theme_bw()
p.1 <- p.1 + xlab("Time") + ylab("Data")

p.2 <- ggplot(data = df, aes(x = x, y = dy))
p.2 <- p.2 + geom_point(alpha=0.4)
p.2 <- p.2 + geom_vline(xintercept = changepoints[-1],color = "blue",linetype = "dashed")
p.2 <- p.2 + geom_vline(xintercept = cpts.full(out)[3,],color = "red",linetype = "dashed")
p.2 <- p.2 + theme_bw()
p.2 <- p.2 + xlab("Time") + ylab("First Differences")

p.3 <- ggplot(data = df, aes(x = x, y = y))
p.3 <- p.3 + geom_point(alpha=0.4)
p.3 <- p.3 + geom_line(aes(y = mu), color = "blue")
p.3 <- p.3 + geom_line(aes(y = fit1), color = "red")
p.3 <- p.3 + geom_vline(xintercept = changepoints[-1],color = "blue",linetype = "dashed")
p.3 <- p.3 + geom_vline(xintercept = (1:399)[abs(diff(diff(fit1)))>0.001]+1,color = "red",linetype = "dashed")
p.3 <- p.3 + theme_bw()
p.3 <- p.3 + xlab("Time") + ylab("Data")

p.4 <- ggplot(data = df, aes(x = x, y = y))
p.4 <- p.4 + geom_point(alpha=0.4)
p.4 <- p.4 + geom_line(aes(y = mu), color = "blue")
p.4 <- p.4 + geom_line(aes(y = fit2), color = "red")
p.4 <- p.4 + geom_vline(xintercept = changepoints[-1],color = "blue",linetype = "dashed")
p.4 <- p.4 + geom_vline(xintercept = (1:399)[abs(diff(diff(fit2)))>0.001]+1,color = "red",linetype = "dashed")
p.4 <- p.4 + theme_bw()
p.4 <- p.4 + xlab("Time") + ylab("Data")

library(gridExtra)

g <- grid.arrange(p.1,p.2,p.3,p.4 ,nrow = 2,ncol=2)

ggsave(file="change_in_slope_examples_ggplot.pdf", g)

## Section 3.1

set.seed(1)
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200
sd <- 0.8
y <- simchangeslope(x, changepoints, change.slope, sd) #simulate data
df <- data.frame("x" = x, "y" = y)
p <- ggplot(data = df, aes(x = x, y = y))
p <- p + geom_point(alpha=0.4)
p <- p + geom_vline(xintercept = changepoints, color = "red",  linetype = "dashed")

#for the plot: obtain the true mean function by setting sd=0.
mu <- simchangeslope(x, changepoints, change.slope, sd = 0)
p <- p + geom_line(aes(y = mu), color = "blue")
p <- p + theme_bw()
print(p)

#### Figure 2

ggsave(file="simulate_example_ggplot.pdf",p)

## Section 3.2

##analysis by cpop
res <- cpop(y, x, beta = 2*log(length(y)), sd = 0.8)
summary(res)

p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)

#### Figure 3

ggsave(file="cpop_example1_ggplot.pdf",p)

## Section 3.3

changepoints(res)

estimate(res, x = c(0.1,2.7,51.6))

fitted(res)

head(residuals(res))

## Section 4.1

set.seed(1)
x <- (1:200)^(2)/(200) ###x-values now unevenly spaced
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
sd <- 0.8
y <- simchangeslope(x, changepoints, change.slope, sd) #simulate data
##analysis by cpop
res <- cpop(y, x, beta = 2*log(length(y)), sd = 0.8)

#### Figure 4

#for the plot: obtain the true mean function by setting sd=0.
mu <- simchangeslope(x, changepoints, change.slope, sd = 0)
p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
print(p)
ggsave(file="cpop_example_uneven_ggplot.pdf",p)


## Section 4.2

set.seed(1)
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200 
sd <- x/100
y <- simchangeslope(x, changepoints, change.slope, sd) #simulate data
##analysis by cpop
res <- cpop(y, x, beta = 2*log(length(y)), sd = sqrt(mean(sd^2)) )
#summary(res)
res.true <- cpop(y, x, beta = 2*log(length(y)), sd=sd)

#### Figure 5

mu <- simchangeslope(x, changepoints, change.slope, sd = 0)
p <- plot(res)
p <- p + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

p.true <- plot(res.true)
p.true <- p.true + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.true <- p.true + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
p <- p + theme(aspect.ratio=1/1)
p.true <- p.true + theme(aspect.ratio=1/1)
g <- grid.arrange(p,p.true,nrow=1,ncol=2)
ggsave(file="cpop_uneven_examples_ggplot.pdf",g)

##  Section 4.3


set.seed(1)
x <- 1:6400
y <- simchangeslope(x, changepoints = 0:31*200, change.slope = c(0.05,0.1*(-1)^(1:(31))), sd = 1)

res.coarse <- cpop(y, x, grid=1:399*16, beta = 2*log(400), sd = 1)
cps <- unlist( changepoints(res.coarse) )

grid <- NULL
for(i in 1:length(cps)){
  grid <- c(grid, cps[i] + (-7):8 )
}

res.fine <- cpop(y, x, grid, beta = 2*log(length(x)), sd = 1)

res<- cpop(y, x, beta = 2*log(length(y)), sd = 1)

##  compare output with default run of cpop 
summary(res.fine) 
summary(res)

#### Figure 6

 set.seed(1)
  n.st <- 200*2^(0:5)
  ##scenario 1 -- one changepoint
  K <- 10 ##number of replications
  time1 <- matrix(NA,nrow=K,ncol=length(n.st))
  time1g <- matrix(NA,nrow=K,ncol=length(n.st))
  for(i in 1:length(n.st)){
    n <- n.st[i]
    for(k in 1:K){
      y <- simchangeslope(1:n,n/2,0.5,1)
      time1[k,i] <- (system.time(cpop(y,beta=2*log(length(y)),sd=1)))[1]
      time1g[k,i] <- (system.time(cpop(y,beta=2*log(length(y)),sd=1,grid=(1:n.st[1])*(n/n.st[1]))))[1]
    }
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
      y <- simchangeslope(1:n,0:(m-1)*n/m,c(0.05,0.1*(-1)^(1:(m-1))),1)
      time2[k,i] <- (system.time(cpop(y,beta=2*log(length(y)),sd=1)))[1]
      time2g[k,i] <- (system.time(cpop(y,beta=2*log(length(y)),sd=1,grid=(1:n.st[1])*(n/n.st[1]))))[1]
    }
  }
#  save(time1,time2,time1g,time2g,file="/Users/paulfearnhead/Dropbox/Apps/Overleaf/JSS: CPOP/time.Rdata")
#    load("/Users/paulfearnhead/Dropbox/Apps/Overleaf/JSS: CPOP/time.Rdata")

###average the times
t1 <- apply(time1,2,mean)
t2 <- apply(time2,2,mean)
t1g <- apply(time1g,2,mean)
t2g <- apply(time2g,2,mean)

df <- data.frame("x" = n.st,
                "t1" = t1,
                "t2" = t2,
                "t1g" = t1g,
                "t2g" = t2g,
                "d1" = t2[1]*n.st^1.7/n.st[1]^1.7,
                "d2" = t1[1]*n.st^2.5/n.st[1]^2.5)
p <- ggplot(data = df, aes(x = x))
p <- p + geom_line(aes(y = t1))
p <- p + geom_line(aes(y = t2),color="red")
p <- p + geom_line(aes(y = t1g),linetype = "dashed")
p <- p + geom_line(aes(y = t2g),linetype = "dashed",color="red")
p <- p + geom_line(aes(y = d1),linetype = "5515",color="red")
p <- p + geom_line(aes(y = d2),linetype = "5515")
p <- p + scale_x_continuous(trans = 'log10')
p <- p + scale_y_continuous(trans = 'log10')
p <- p + xlab("n") + ylab("CPU time (sec)")
p <- p + theme_bw()
print(p)
ggsave(file="cpop_CPU_ggplot.pdf", p)

## Section 4.4

set.seed(1)
##simulate data as for Section 3
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200
mu <- simchangeslope(x, changepoints, change.slope, sd = 0) #mean

##simulate data with t noise
y <- mu + rt(length(x), df = 4)
res.min <- cpop(y,x, beta=2*log(length(y)), minseglen = 10, sd = sqrt(2)) ##variance of t_d is d/(d-2)

#### Figure 7

res <- cpop(y,x, beta=2*log(length(y)), sd = sqrt(2)) ##variance of t_d is d/(d-2)

##Now run with a minimum segment length
#res.min=cpop(y,x,beta=2*(log(length(y))),minseglen=10, sd = sqrt(4/2))
##compare CPU cost
system.time(cpop(y,x, beta=2*log(length(y)), sd = sqrt(2) ) )[1]
system.time(cpop(y,x,beta=2*(log(length(y))),minseglen=10, sd = sqrt(2) ) )[1]

##run with too large a minimum segment length
res.min40=cpop(y,x,beta=2*(log(length(y))),minseglen=40, sd = sqrt(2))
res.min30=cpop(y,x,beta=2*(log(length(y))),minseglen=30, sd = sqrt(2))



p.1 <- plot(res)
p.1 <- p.1 + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.1 <- p.1 + geom_line(aes(y = mu), color = "blue", linetype = "dashed")


p.2 <- plot(res.min)
p.2 <- p.2 + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.2 <- p.2 + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

p.3 <- plot(res.min30)
p.3 <- p.3 + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.3 <- p.3 + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

p.4 <- plot(res.min40)
p.4 <- p.4 + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.4 <- p.4 + geom_line(aes(y = mu), color = "blue", linetype = "dashed")


g <- grid.arrange(p.1,p.2,p.3,p.4 ,nrow = 2,ncol=2)

ggsave(file="cpop_minseg_ggplot.pdf", g)

## Section 4.5

set.seed(1)
##simulate data as for Section 3
changepoints <- c(0, 25, 50, 100)
change.slope <- c(0.2, -0.3, 0.2, -0.1)
x <- 1:200
mu <- simchangeslope(x, changepoints, change.slope,  0) #mean
y <- simchangeslope(x, changepoints, change.slope, 1.5) #data with sd=1.5
###run crops
res.crops <- cpop.crops(y , x, beta_min= 5 ,beta_max= 50, sd = 1)
plot(res.crops)

##calculate the BIC under a model of unknown variance
models <- cpop.crops.models(res.crops)
M <- length(models)
BIC <- rep(NA, M)
ncps <- segmentations( res.crops )[,4]
n <- length(y)
for(i in 1:M){
  BIC[i] <- n*log( mean( ( residuals(models[[i]]) )^2) )+ 2 * ncps[i]  * log(n)
}

#### Figure 8

p.1 <- plot(res.crops)
i <- which.min(BIC)
p.2 <- plot(models[[i]])
p.2 <- p.2 + geom_vline(xintercept = changepoints[-1], color = "blue", linetype = "dashed")
p.2 <- p.2 + geom_line(aes(y = mu), color = "blue", linetype = "dashed")
p.1 <- p.1 + theme(aspect.ratio=1/1)
p.2 <- p.2 + theme(aspect.ratio=1/1)
library(cowplot)
g <- plot_grid(p.1, p.2, align = "v", nrow = 1)
ggsave(file="cpop_crops_example_ggplot.pdf",g)
plot(g)

set.seed(1)
n <- 500
x <- 1:n
mu <- simchangeslope(x,changepoints = 45*0:10,change.slope = c(0.15,0.3*(-1)^(1:10)) ,sd = 0)
epsilon <- rnorm(n+2)
y <- mu + (epsilon[1:n] + epsilon[2:(n+1)] + epsilon[3:(n+2)]) /sqrt(3)

res.crops <- cpop.crops(y,x,beta_min=8,beta_max=200,sd=1)
segs <- segmentations(res.crops)

p <- ggplot(data = segs, aes(x = m))
p <- p + geom_line(aes(y = Qm))
p <- p + geom_vline(xintercept = 10,color = "red")
p <- p + xlab("No. of changepoints") + ylab("unpenalised cost")
plot(p)

### Figure 9

set.seed(1)
n <- 500
x <- 1:n
mu <- simchangeslope(x,changepoints = 45*0:10,change.slope = c(0.15,0.3*(-1)^(1:10)) ,sd = 0)
epsilon <- rnorm(n+2)
y <- mu + (epsilon[1:n] + epsilon[2:(n+1)] + epsilon[3:(n+2)]) /sqrt(3)
p.1 <- plot(res.crops)
p.2 <- p + theme_bw()
models <- cpop.crops.models(res.crops)
p <- plot(models[[5]])
p <- p + geom_vline(xintercept = 45*1:10, color = "blue", linetype = "dashed")
p.3 <- p + geom_line(aes(y = mu), color = "blue", linetype = "dashed")

g <- grid.arrange(p.1,p.2,p.3,layout_matrix=rbind(c(1,2),c(3,3)))
ggsave(file="cpop_crops_ggplot.pdf",g)

## Section 5 

data("wavenumber_spectra") ##Load data

##example for one curve -- fitting the decay.
x <- log(wavenumber_spectra[-(1:3),1],base=10)
y <- log(wavenumber_spectra[-(1:3),4],base=10)

grid <- seq(from=min(x),to=max(x),length=200)
##naive estimator of variance
sig2 <- mean( diff( diff(y) )^2 )/6 ##simple estimate of variance 

res <- cpop(y, x, grid, sd=sqrt(sig2), minseglen = 0.09, beta=2*log(200))

r2 <- residuals(res)^2
##estimate variance as exp(a+bx)
##function for minus twice log-likelihood
loglik <- function(par){ 
  return(length(r2) * par[1] + par[2] * sum(x) + sum( r2/ (exp(par[1]+par[2]*x) )) )
}

est.hat <- optim( c(0,0) , loglik)
sig2 <- exp(est.hat$par[1] + est.hat$par[2]*x)

res2 <- cpop(y, x, grid , sd=sqrt(sig2), minseglen = 0.09, beta=2 * log(200))

#### Figure 10

library(scales)
library(latex2exp)

######FUNCTION TO AUTOMATE THIS ANALYSIS
wavenumber_est <- function(x,y){
  grid <- seq(from=min(x),to=max(x),length=200)
  ##naive estimator of variance
  sig2 <- mean( diff( diff(y) )^2 )/6 ##simple estimate of variance 
  
  res <- cpop(y, x, grid, sd=sqrt(sig2), minseglen = 0.09, beta=2*log(200))
 
  r2 <- residuals(res)^2
  est.hat <- optim( c(0,0) , loglik)
  sig2 <- exp(est.hat$par[1] + est.hat$par[2]*x)
  
  res2 <- cpop(y, x, grid , sd=sqrt(sig2), minseglen = 0.09, beta=2 * log(200))
  
  return(res2)
}

y.est <- matrix(NA,nrow=4,ncol=length(x))
for(j in 1:4){
  est <- wavenumber_est(x, log(wavenumber_spectra[-(1:3),1+j]))
  y.est[j, ] <- estimate(est,x)[,2]
}

p <- ggplot(data = wavenumber_spectra, aes(x = wavenumber))
p <- p + geom_line(aes(y = power_spectra_Feb2000),color="black")
p <- p + geom_line(aes(y = power_spectra_Aug2000),color="red")
p <- p + geom_line(aes(y = power_spectra_Feb2100),color="green")
p <- p + geom_line(aes(y = power_spectra_Aug2100),color="blue")
# p <- p + scale_x_continuous(trans = 'log10')
p <- p + scale_x_continuous(trans = 'log10',breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))
p <- p + scale_y_continuous(trans = 'log10',breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))
# p <- p + xlab("wavenumber (cyc/m)") + ylab("log_10 spectra")
# p <- p + xlab("wavenumber (cyc/m)") +    ylab(TeX("$\\log_{10}$ spectra"))
p <- p + xlab("wavenumber (cyc/m)") +    ylab("spectra")                            
p.1 <- p + theme_bw()

p.2 <- plot(res)
p.2 <- p.2 + xlab(TeX("$\\log_{10}(wavenumber)$")) +    ylab(TeX("$\\log_{10}(spectra)$"))                            

p.3 <- plot(res2)
p.3 <- p.3 + xlab(TeX("$\\log_{10}(wavenumber)$")) +    ylab(TeX("$\\log_{10}(spectra)$")) 
#p.3 <- p.3 + scale_y_continuous(labels=function(.) sprintf("%.1f", .))

df <- wavenumber_spectra[-(1:3),]
p <- ggplot(data = df, aes(x = wavenumber))
p <- p + geom_line(aes(y = power_spectra_Feb2000),color="black")
p <- p + geom_line(aes(y = power_spectra_Aug2000),color="red")
p <- p + geom_line(aes(y = power_spectra_Feb2100),color="green")
p <- p + geom_line(aes(y = power_spectra_Aug2100),color="blue")
p <- p + geom_line(aes(y = exp(y.est[1,])),color="black",linetype="dashed")
p <- p + geom_line(aes(y = exp(y.est[2,])),color="red",linetype="dashed")
p <- p + geom_line(aes(y = exp(y.est[3,])),color="green",linetype="dashed")
p <- p + geom_line(aes(y = exp(y.est[4,])),color="blue",linetype="dashed")
#p <- p + scale_x_continuous(trans = 'log10')
p <- p + scale_x_continuous(trans = 'log10',breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))
p <- p + scale_y_continuous(trans = 'log10',breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))
p <- p + xlab("wavenumber (cyc/m)") + ylab("spectra")
p.4 <- p + theme_bw()

g <- grid.arrange(p.1,p.2,p.3,p.4,nrow=2,ncol=2)

print(g)
ggsave(file="cpop_real_data_ggplot.pdf",g)
