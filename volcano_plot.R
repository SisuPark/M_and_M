
cexFn <- function(x){ -0.75*(x-2) + 2}

volcano.plot <- function(signal.mat, q.thld=0.05, e.thld=0.1,
                        legend.title="Median Relative Abundance",
                        xlab=expression(paste("max ", Delta, " pr( sPTB )")), ylab="-log10( q-values )",
                        xlim=c(-1.1,1.1), ylim=c(0,16), lfactor=1, rfactor=1.4,
                        yjust=-0.15, xjust=-0.15, adj=c(-0.07, -0.6), text.cex=0.8)
{
    op <- par(mar=c(4, 4, 4, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3) #cex.axis=0.7)
    x <- signal.mat[,1]
    idx <- signal.mat[,3] < 1e-14
    if ( sum(idx) )
    {
        signal.mat[idx,3] <- runif(sum(idx), 1e-16, 1e-14)
    }
    y <- -log10(signal.mat[,3])
    if ( any(is.na(xlim)) )
    {
        r <- range(x)
        xlim <- c(lfactor*r[1], rfactor*r[2])
    }
    if ( any(is.na(ylim)) )
    {
        ylim <- range(y[is.finite(y)])
    }
    plot(x,y, las=1, xlab=xlab, ylab=ylab, pch=19, xlim=xlim, ylim=ylim, type='n')
    abline(h=-log10(q.thld), col='gray80')
    abline(v=0, col='gray80')
    for ( i in seq(nrow(signal.mat)) )
    {
        if ( signal.mat[i,3] < q.thld & abs(signal.mat[i,1]) > e.thld )
        {
            if ( x[i] > 0 ){
                col <- "red"
            } else {
                col <- "green"
            }
            text(x[i], y[i], labels=gsub("_"," ",rownames(signal.mat)[i]), adj=adj, cex=text.cex, font=3)
            points(x[i], y[i], pch=19, col=col, cex=cexFn(signal.mat[i,ncol(signal.mat)]))
        } else {
            points(x[i], y[i], pch=20)
        }
    }
    legend(par('usr')[1], par('usr')[4], yjust=yjust, xjust=xjust, xpd=NA,
           legend=c("0.001","0.01","0.1"), pch=19, pt.cex=cexFn(c(3,2,1)),ncol=3,
           title=legend.title)
    par(op)
}


poly2.volcano.plot <- function(signal.mat, q.thld=0.05,
                              xlab="log OR of pr( sPTB )", ylab="-log10( q-values )",
                              xlim1=NA, xlim2=NA, ylim=c(0,16), lfactor=1, rfactor=2.4,
                              yjust=-0.15, xjust=-0.15, adj=c(-0.07, -0.6), text.cex=0.8)
{
    op <- par(mfrow=c(1,2), mar=c(4, 4, 4, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3) #cex.axis=0.7)
    ## linear term
    x <- signal.mat[,1]
    idx <- signal.mat[,3] < 1e-14
    if ( sum(idx) )
    {
        signal.mat[idx,3] <- runif(sum(idx), 1e-16, 1e-14)
    }
    y <- -log10(signal.mat[,3])
    if ( any(is.na(xlim1)) )
    {
        r <- range(x)
        xlim1 <- c(lfactor*r[1], rfactor*r[2])
    }
    if ( any(is.na(ylim)) )
    {
        ylim <- range(y[is.finite(y)])
    }
    plot(x,y, las=1, xlab=xlab, ylab=ylab, pch=19, xlim=xlim1, ylim=ylim, type='n')
    abline(h=-log10(q.thld), col='gray80')
    abline(v=0, col='gray80')
    for ( i in seq(nrow(signal.mat)) )
    {
        if ( signal.mat[i,3] < q.thld )
        {
            if ( x[i] > 0 ){
                col <- "red"
            } else {
                col <- "green"
            }
            text(x[i], y[i], labels=gsub("_"," ",rownames(signal.mat)[i]), adj=adj, cex=text.cex, font=3)
            points(x[i], y[i], pch=19, col=col, cex=cexFn(signal.mat[i,ncol(signal.mat)]))
        } else {
            points(x[i], y[i], pch=20)
        }
    }
    legend(par('usr')[1], par('usr')[4], yjust=yjust, xjust=xjust, xpd=NA,
           legend=c("0.001","0.01","0.1"), pch=19, pt.cex=cexFn(c(3,2,1)),ncol=3,
           title="Median Relative Abundance")
    ## quadratic term
    x <- signal.mat[,1+3]
    idx <- signal.mat[,3+3] < 1e-14
    if ( sum(idx) )
    {
        signal.mat[idx,3+3] <- runif(sum(idx), 1e-16, 1e-14)
    }
    y <- -log10(signal.mat[,3+3])
    if ( any(is.na(xlim2)) )
    {
        r <- range(x)
        xlim2 <- c(lfactor*r[1], rfactor*r[2])
    }
    if ( any(is.na(ylim)) )
    {
        ylim <- range(y[is.finite(y)])
    }
    plot(x,y, las=1, xlab=xlab, ylab=ylab, pch=19, xlim=xlim2, ylim=ylim, type='n')
    abline(h=-log10(q.thld), col='gray80')
    abline(v=0, col='gray80')
    for ( i in seq(nrow(signal.mat)) )
    {
        if ( signal.mat[i,3+3] < q.thld )
        {
            if ( x[i] > 0 ){
                col <- "red"
            } else {
                col <- "green"
            }
            text(x[i], y[i], labels=gsub("_"," ",rownames(signal.mat)[i]), adj=adj, cex=text.cex, font=3)
            points(x[i], y[i], pch=19, col=col, cex=cexFn(signal.mat[i,ncol(signal.mat)]))
        } else {
            points(x[i], y[i], pch=20)
        }
    }
    legend(par('usr')[1], par('usr')[4], yjust=yjust, xjust=xjust, xpd=NA,
           legend=c("0.001","0.01","0.1"), pch=19, pt.cex=cexFn(c(3,2,1)),ncol=3,
           title="Median Relative Abundance")
    par(op)
}
