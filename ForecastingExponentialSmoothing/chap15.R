require(expsmooth)
require(Mcomp)

source("chap15_functions.R")

## AMN
set.seed(8)
alpha = 0.1
beta = 0.05
e <- rnorm(100)
y <- numeric(100)
l <- b <- numeric(101)
l[1] <- 0.1
b[1] <- 1
for(i in 1:100)
{
    y[i] <- l[i]*b[i]+ e[i]
    l[i+1] <- l[i]*b[i] + alpha*e[i]
    b[i+1] <- b[i] + beta*e[i]/l[i]
    if(i==5)
        l[i+1] <- 3e-4
}


par(mfrow=c(3,1))
plot(ts(c(NA,y),s=0),xlim=c(0,9.6),ylim=c(-500,500),ylab="y")
plot(ts(l[1:11],s=0),xlim=c(0,9.6),ylim=c(-1,1),ylab="ell");abline(0,0,lty=3)
plot(ts(b,s=0),ylab="b",xlim=c(0,9.6))



# MNN Model, whose sample paths converge to zero with normal and non-normal errors

mod <- ets(rnorm(3)+10,model="MNN",alpha=0.3)
mod$sigma2 <- .3^2

par(mfrow=c(3,1))
set.seed(4);y <- simulate.ets(mod, 2000, initstate=10, bootstrap=FALSE); plot(y,)
y <- simulate.ets(mod, 2000, initstate=10, bootstrap=FALSE); plot(y)
y <- simulate.ets(mod, 2000, initstate=10, bootstrap=FALSE); plot(y)


# LOG NORMAL
set.seed(100)
w <- 0.05
alpha <- 0.3
case <- 1
cc <- alpha*w/4
n <- 2000
y1 <- simulate.mnn(n,alpha,1,"lognormal",c(cc,w),case)$y
y2 <- simulate.mnn(n,alpha,1,"lognormal",c(0,w),case)$y
y3 <- simulate.mnn(n,alpha,1,"lognormal",c(-cc,w),case)$y
y4 <- simulate.mnn(n,alpha,1,"lognormal",c(-cc*1.5,w),case)$y
y5 <- simulate.mnn(n,alpha,1,"lognormal",c(-2*cc,w),case)$y
y6 <- simulate.mnn(n,alpha,1,"lognormal",c(-3*w/5,w),case)$y


par(mfrow=c(3,2))
plot(y1,ylab="(a)",xlab="t")
plot(y2,ylab="(b)",xlab="t")
plot(y3,ylab="(c)",xlab="t")
plot(y4,ylab="(d)",xlab="t")
plot(y5,ylab="(e)",xlab="t")
plot(y6,ylab="(f)",xlab="t")


#Freight

plot(ts(c(M3$N0193$x,M3$N0193$xx),s=1947),xlab="Year",ylab="New freight cars")


## JEWELRY DATA
yy=jewelry
n=nrow(yy)
NS <- ncol(yy)
h=25
#opt1=opt2=opt3=ss=list(NS)
#ss1=ss2=ss3=list(NS)
sigma1=sigma2=sigma3=y=numeric(NS)
mape1=mape2=mape3=y=numeric(NS)
mase1=mase2=mase3=y=numeric(NS)

for(j in 1:NS)
{
    y=ts(yy[1:(n-h),j])
    yhout=ts(yy[(n-h+1):n,j])

    # Find starting point
    fit <- ses(y)
    alpha <- fit$model$par[[1]]
    init.lev <- fit$model$par[[2]]
    startparam <- c(alpha,init.lev)

    opt1 <- nlm(lik,startparam,y=y,model=1)
    ss1 <- fcast.mnn(yhout,opt1$estimate[1],opt1$estimate[2],model=1)
    mase1[j] = ss1$mase
    sigma1[j] <- ss1$mse
    mape1[j]=ss1$mape
    opt2 <- nlm(lik,startparam,y=y,model=2)
    ss2 <- fcast.mnn(yhout,opt2$estimate[1],opt2$estimate[2],model=2)
    mase2[j] = ss2$mase
    sigma2[j] <- ss2$mse
    mape2[j]=ss2$mape
    opt3 <- nlm(lik,startparam,y=y,model=3)
    ss3 <- fcast.mnn(yhout,opt3$estimate[1],opt3$estimate[2],model=3)
    mase3[j] = ss3$mase
    sigma3[j] <- ss3$mse
    mape3[j]=ss3$mape
}


## Plotting MASE

ul=max(mase1,mase2,mase3)
ll=min(mase1,mase2,mase3)
j <- order(apply(cbind(mase1,mase2,mase3),1,median))
plot(mase1[j],xlab="Series",ylab="MASE",main="",ylim=c(ll,ul))
points(mase2[j],pch=2)
points(mase3[j],pch=3)
abline(h=mean(mase1),lty=1)
abline(h=mean(mase2),lty=2)
abline(h=mean(mase3),lty=3)
legend("topleft",legend=c("Model 1","Model 2","Model 3"),lty=1:3,pch=1:3,merge=FALSE)


par(mfrow=c(3,1),mar=c(4,4,2,1))
x1=y1=seq(ll-0.1,ul+.1,0.01)
plot(mase1,mase3,xlim=c(ll,ul),ylim=c(ll,ul),xlab="MASE of model 1",ylab="MASE of model 3",main="")
lines(x1,y1,col="gray")
#text(3,3,"MASE 1 = MASE 3")
plot(mase1,mase2,xlim=c(ll,ul),ylim=c(ll,ul),xlab="MASE of model 1",ylab="MASE of model 2",main="")
lines(x1,y1,col="gray")
#text(3,3,"MASE 1 = MASE 2")
plot(mase2,mase3,xlim=c(ll,ul),ylim=c(ll,ul),xlab="MASE of model 2",ylab="MASE of model 3",main="")
lines(x1,y1,col="gray")
#text(3,3,"MASE 2 = MASE 3")

