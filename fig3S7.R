

library(rstan)
load("bin_2prop.rda") # b2p

nChains=10; nItr=50000; thin=50

## computing log ratio, 95% CI and the corresponding p-value of two proportions
## from m.dat
m.2prop.Fn <- function(m.dat, nChains=4, nItr=5000, thin=5)
{
    b2p.fit <- sampling(b2p, data=m.dat, iter=nItr, thin=thin, chains=nChains)
    m <- summary(b2p.fit, pars=c("logRat"), probs=c(0.025,0.975))$summary
    logRat.ci <- m[,c(1,4,5)]

    logRat <- as.vector(extract(b2p.fit, pars=c("logRat"), permuted = TRUE)$logRat)
    ## str(logRat)
    qqnorm(logRat)
    qqline(logRat, col=2)
    myHist(logRat)

    logRat.bar <- median(logRat)

    logRat.pFn <- function(p)
    {
        ret <- 0
        p2 <- p/2
        mu.q <- quantile(logRat, probs=c(p2,1-p2))
        if ( logRat.bar > 0 )
        {
            ret <- mu.q[1]
        }
        else
        {
            ret <- mu.q[2]
        }
        ret
    }

    pFn <- function(p)
    {
        ret <- 0
        ##p2 <- p/2
        p2 <- p
        mu.q <- quantile(logRat, probs=c(p2,1-p2))
        if ( logRat.bar > 0 )
        {
            ret <- mu.q[1]
        }
        else
        {
            ret <- mu.q[2]
        }
        ret
    }

    fn.pval <- NA
    if ( sign(pFn(0)*pFn(1))==-1 )
    {
        fn.pval <- uniroot(pFn, c(0, 1))$root
    }

    logRat.pval <- NA
    if (  logRat.bar < 0 )
    {
        logRat.pval <- 1 - pnorm(0, logRat.bar, mad(logRat))
    } else {
        logRat.pval <- pnorm(0, logRat.bar, mad(logRat))
    }
    logRat.pval

    logRat.ci <- c(logRat.ci, logRat.pval)
    names(logRat.ci)[1] <- "logRat"
    names(logRat.ci)[4] <- "p-value"

    m <- summary(b2p.fit, pars=c("p"), probs=c(0.025,0.975))$summary
    p.ci <- m[,c(1,4,5)]

    p1 <- m.dat$y[2]/m.dat$n[2]
    p2 <- m.dat$y[1]/m.dat$n[1]

    list(logRat.ci=logRat.ci, fn.pval=fn.pval, est.rat=exp(logRat.ci[1]), ML.rat=p1/p2, p.ci, p1=p1, p2=p2, m.dat=m.dat) # , m.fit=m.fit
}


m.Fn <- function(n, n.sPTB, all.n, all.n.sPTB, nChains=4, nItr=5000, thin=5)
{
    m.dat <- list(y=c(all.n.sPTB, n.sPTB),
                 n=c(all.n, n))
    str(m.dat)

    r <- m.2prop.Fn(m.dat, nChains=nChains, nItr=nItr, thin=thin)
    r
}


## All subjects at V1

m.dat <- list(y=c(102, 11),
             n=c(521, 68))
str(m.dat)

(r <- m.2prop.Fn(m.dat, nChains=nChains, nItr=nItr, thin=thin))
## $logRat.ci
##     logRat       2.5%      97.5%    p-value
## -0.1682949 -0.7693194  0.3508604  0.2983264

## $fn.pval
## [1] 0.288366

## $est.rat
##    logRat
## 0.8451046

## $ML.rat
## [1] 0.8262687

## [[5]]
##           mean       2.5%     97.5%
## p[1] 0.1972283 0.16423251 0.2333074
## p[2] 0.1719663 0.09367827 0.2689715

## $p1
## [1] 0.1617647

## $p2
## [1] 0.1957774

## $m.dat
## $m.dat$y
## [1] 102  11

## $m.dat$n
## [1] 521  68
