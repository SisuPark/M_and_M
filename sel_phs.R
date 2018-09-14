
sel.phs2 <- function(idx2, phs, ref.ph, minCaseSubjs=10)
{
    selPhs <- c()
    for ( ph in phs )
    {
        r <- n.stats2(ph, ref.ph, idx1)
        ##print(r$nCaseSubjs)
        if ( r$nCaseSubjs > minCaseSubjs )
        {
            selPhs <- c(selPhs, ph)
        }
    }
    selPhs
}

sel.phs <- function(idx1, thld=0.25)
{
    y <- mt$delStatus.n
    s <- mt$subjID
    #idx1 <- mt$delStatus.n!=2 & !is.na(mt$bLoad)
    y <- y[idx1]
    s <- s[idx1]
    ##(tt <- table(y))

    n.TERM <- sum(y==0)
    n.sPTB <- sum(y==1)

    n.subj.TERM <- length(unique(s[y==0]))
    n.subj.sPTB <- length(unique(s[y==1]))

    min.n.case.subjs <- thld * n.subj.sPTB

    selPhs <- c()
    for ( ph in selPh2.10 )
    {
        x <- pt2[,ph]
        y <- mt$delStatus.n
        s <- mt$subjID
        #idx1 <- mt$delStatus.n!=2 & !is.na(mt$bLoad)
        x <- x[idx1]
        y <- y[idx1]
        s <- s[idx1]
        idx <- x > 0
        x <- log10(x[idx])
        y <- y[idx]
        s <- s[idx]
        q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
        idx <- x > q[1] & x < q[2]
        x <- x[idx]
        y <- y[idx]
        s <- s[idx]
        n.case.subjs <- length(unique(s[y==1]))
        if ( n.case.subjs >= min.n.case.subjs )
        {
            selPhs <- c(selPhs, ph)
        }
    }

    list(selPhs=selPhs,
         n.subj.TERM=n.subj.TERM,
         n.subj.sPTB=n.subj.sPTB,
         n.TERM=n.TERM,
         n.sPTB=n.sPTB,
         min.n.case.subjs=min.n.case.subjs)
}


sel.phs.silva <- function(idx1, thld=0.25)
{
    y <- mt$delStatus.n
    s <- mt$subjID
    #idx1 <- mt$delStatus.n!=2 & !is.na(mt$bLoad)
    y <- y[idx1]
    s <- s[idx1]
    ##(tt <- table(y))

    n.TERM <- sum(y==0)
    n.sPTB <- sum(y==1)

    n.subj.TERM <- length(unique(s[y==0]))
    n.subj.sPTB <- length(unique(s[y==1]))

    min.n.case.subjs <- thld * n.subj.sPTB

    selPhs <- c()
    for ( ph in colnames(pt.silva) )
    {
        x <- pt.silva[,ph]
        y <- mt$delStatus.n
        s <- mt$subjID
        #idx1 <- mt$delStatus.n!=2 & !is.na(mt$bLoad)
        x <- x[idx1]
        y <- y[idx1]
        s <- s[idx1]
        idx <- x > 0
        x <- log10(x[idx])
        y <- y[idx]
        s <- s[idx]
        q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
        idx <- x > q[1] & x < q[2]
        x <- x[idx]
        y <- y[idx]
        s <- s[idx]
        n.case.subjs <- length(unique(s[y==1]))
        if ( n.case.subjs >= min.n.case.subjs )
        {
            selPhs <- c(selPhs, ph)
        }
    }

    list(selPhs=selPhs,
         n.subj.TERM=n.subj.TERM,
         n.subj.sPTB=n.subj.sPTB,
         n.TERM=n.TERM,
         n.sPTB=n.sPTB,
         min.n.case.subjs=min.n.case.subjs)
}


sel.phs.relman <- function(thld=0.25)
{
    y <- mt.relman$preterm.n
    s <- mt.relman$subjID
    ##(tt <- table(y))

    n.TERM <- sum(y==0)
    n.sPTB <- sum(y==1)

    n.subj.TERM <- length(unique(s[y==0]))
    n.subj.sPTB <- length(unique(s[y==1]))

    min.n.case.subjs <- thld * n.subj.sPTB

    selPhs <- c()
    for ( ph in colnames(pt.relman) )
    {
        x <- pt.relman[,ph]
        y <- mt.relman$preterm.n
        s <- mt.relman$subjID
        idx <- x > 0
        x <- log10(x[idx])
        y <- y[idx]
        s <- s[idx]
        q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
        idx <- x > q[1] & x < q[2]
        x <- x[idx]
        y <- y[idx]
        s <- s[idx]
        n.case.subjs <- length(unique(s[y==1]))
        if ( n.case.subjs >= min.n.case.subjs )
        {
            selPhs <- c(selPhs, ph)
        }
    }

    list(selPhs=selPhs,
         n.subj.TERM=n.subj.TERM,
         n.subj.sPTB=n.subj.sPTB,
         n.TERM=n.TERM,
         n.sPTB=n.sPTB,
         min.n.case.subjs=min.n.case.subjs)
}

## (n.TERM.10 <- as.integer(n.TERM*thld)) #
## (n.sPTB.10 <- as.integer(n.sPTB*thld)) #

## sel.10 <- c()
## sel.10.list <- list()
## for ( ph in selPh2.10 )
## {
##     x <- pt2[,ph]
##     y <- mt$delStatus.n
##     idx1 <- mt$delStatus.n!=2
##     x <- x[idx1]
##     y <- y[idx1]
##     idx <- x > 0
##     x <- log10(x[idx])
##     y <- y[idx]
##     s <- s[idx]
##     q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
##     idx <- x > q[1] & x < q[2]
##     x <- x[idx]
##     y <- y[idx]
##     s <- s[idx]
##     if ( length(x[y==0]) > n.TERM.10 && length(x[y==1]) > n.sPTB.10 )
##     {
##         sel.10 <- c(sel.10, ph)
##         sel.10.list[[ph]] <- c(length(x[y==0]), length(x[y==1]))
##     }
## }

## length(selPhs)
