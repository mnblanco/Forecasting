## Likelihood calculation
lik <- function(param,y,model)
{
    alpha=param[1]
    init.lev=param[2]
    
    if(alpha<=0 | alpha>1)
    return(1e11) #Return huge number

    # Compute error series
    ssm=est.mnn(y,alpha,init.lev,model)

    # Return -likelihood
    return(ssm$llik)
}

## estimation part
est.mnn <- function(y,alpha,init.lev,model)
{
    n <- length(y)
    lev <- yhat <- eps <- numeric(n)
    lev[1] <- init.lev

        for(i in 1:n)
        {
            if(model==1)
            {
                yhat[i] <- lev[i]
                eps[i] <- (y[i]-yhat[i])/yhat[i]
                lev[i+1] <- lev[i]*(1+alpha*eps[i])
            }
            else if(model==2)
            {
                yhat[i] <- lev[i]
                eps[i] <- y[i]/yhat[i]
                lev[i+1] <- lev[i]*(1-alpha+eps[i]*alpha)
            }
            else if(model==3)
            {
                yhat[i] <- lev[i]
                eps[i] <- y[i]/yhat[i]
                lev[i+1] <- lev[i]*eps[i]^alpha
            }
            else stop("Unknown Model")
        }
    if(model==1)
    {
        mse=mean(eps^2)
        llik=n*log(sum(eps^2))+2*sum(log(abs(lev)))
    }
    else
    {
        mse=mean((eps-1)^2)
        llik=n*log(sum((log(eps)^2)))+2*sum(log(eps))+2*sum(log(abs(lev)))
    }
    return(list(mse=mse,llik=llik))
}

## forecast part
fcast.mnn <- function(y,alpha,init.lev,model)
{
    n <- length(y)
    lev <- yhat <- eps <- numeric(n)
    lev[1] <- init.lev

        for(i in 1:n)
        {
            if(model==1)
            {
                yhat[i] <- lev[i]
                eps[i] <- (y[i]-yhat[i])/yhat[i]
                lev[i+1] <- lev[i]*(1+alpha*eps[i])
            }
            else if(model==2)
            {
                yhat[i] <- lev[i]
                eps[i] <- y[i]/yhat[i]
                lev[i+1] <- lev[i]*(1-alpha+eps[i]*alpha)
            }
            else if(model==3)
            {
                yhat[i] <- lev[i]
                eps[i] <- y[i]/yhat[i]
                lev[i+1] <- lev[i]*eps[i]^alpha
            }
            else stop("Unknown Model")
        }

        eps1=y-yhat
        mae=mean(abs(diff(y)))
        mase=mean(abs(eps1/mae))
        mape=mean(abs(eps1/y*100))
        mse=mean(eps1^2)

    return(list(yhat=yhat,eps=eps,lev=lev,mse=mse,mape=mape,mase=mase))
}

## MNN

simulate.mnn = function(n,alpha,ilev,error,param,case=NA)
{
    if(error=="normal")
        eps <- rnorm(n,param[1],param[2])
    else if(error=="lognormal")
        eps <- rlnorm(n,param[1],param[2])
    else if(error=="gamma")
        eps <- rgamma(n,shape=param[1],scale=param[2])
    else
        stop("Unknown distribution")
    l <- y <- numeric(n+1)
    y[1] <- l[1] <- ilev

    for(i in 1:n)
    {
        if(error=="normal")
        {
            l[i+1] <- l[i]*(1+alpha*eps[i])
            y[i+1] <- l[i]*(1+eps[i])
        }
        else
        {
            if(case==1)
            {
                l[i+1] <- l[i]*eps[i]^alpha
                y[i+1] <- l[i]*eps[i]
            }
            else
            {
                l[i+1] <- l[i]*(1-alpha+eps[i]*alpha)
                y[i+1] <- l[i]*eps[i]
            }
        }
    }
    y <- ts(y[-1])
    return(list(y=y,eps=eps,l=l))
}

