
library(rstan)
source("stan_utils.R")

## mt.v1, ct.v1, bdt.v1,
## mt.v1.AA, ct.v1.AA, bdt.v1.AA,
## mt.v1.nonAA, ct.v1.nonAA, bdt.v1.nonAA,
load("bd_v1.rda")

load("mm_no_Thermodesulfatator_v3.rda") # mt, ct, pt, selPh2.10


## V1 CSTs

cst.q <- matrix(NA,nrow=length(CSTs), ncol=5)
rownames(cst.q) <- CSTs
colnames(cst.q) <- c("0\\%","25\\%", "50\\%", "75\\%","100\\%")
for ( cst in CSTs )
{
    idx <- mt.v1$cst==cst
    y <- log10(bdt.v1$BD_v1[idx])
    q <- quantile(y, probs = c(0,0.25,0.5, 0.75,1))
    cst.q[cst,] <- q
}

cst.q

write.csv(cst.q, file="cst_bd_quartiles.csv", quote=F)


## AA V1 phylotypes

phs.q <- matrix(NA,nrow=length(phs), ncol=5)
rownames(phs.q) <- phs
colnames(phs.q) <- c("0\\%","25\\%", "50\\%", "75\\%","100\\%")
for ( ph in phs )
{
    idx <- ct.v1[,ph] > 0
    y <- log10(bdt.v1$BD_v1[idx])
    q <- quantile(y, probs = c(0,0.25,0.5, 0.75,1))
    phs.q[ph,] <- q
}

phs.q

write.csv(phs.q, file="phs_bd_quartiles.csv", quote=F)



v1.cst.fn <- function(cst, nItr=50000, thin=5, nChains=4, nCores=4)
{
    idx <- mt.v1$cst==cst
    y <- log10(bdt.v1$BD_v1[idx])
    x <- mt.v1$sPTB[idx]
    idx <- !is.na(y) & !is.na(x)
    y <- y[idx]
    x <- x[idx]
    x <- ifelse(x==0,"TERM","sPTB")
    x <- factor(x, levels=c("sPTB","TERM"))
    qq <- get.qq(y)

    (f <- table(qq,x))
    (p <- 100*f[,1]/rowSums(f))
    (n <- rowSums(f))
    m.dat <- list(y=as.vector(f[c(1,4),1]),
                 n=as.vector(n[c(1,4)]))
    str(m.dat)

    r <- m.2prop.fn(m.dat, nItr=nItr, thin=thin, nChains=nChains, nCores=nCores)
    s <- c(paste0(f[1,1], "/", sum(f[1,]), " (", sprintf("%.1f",p[1]),")"),
          paste0(f[4,1], "/", sum(f[4,]), " (", sprintf("%.1f",p[4]),")"),
          sprintf("%.3f",r$logRat.ci[4]))
}

AA.bd.v1.cst.tbl <- matrix(NA,nrow=length(names(table(mt.v1$cst))), ncol=3)
colnames(AA.bd.v1.cst.tbl) <- c("Q1 sPTB/TERM (\\%)", "Q4 sPTB/TERM (\\%)", "p-val")
rownames(AA.bd.v1.cst.tbl) <- names(table(mt.v1$cst))
for ( cst in rownames(AA.bd.v1.cst.tbl) )
{
    print(cst)
    r <- v1.cst.fn(cst, nCores=1)
    AA.bd.v1.cst.tbl[cst,] <- r
}

AA.bd.v1.cst.tbl


tblTexFile <- "../docs/Tables/q_bd_V1_CST.tex"
tblLabel <- "q.bd.v1.cst:tbl"
lt <- latex(AA.bd.v1.cst.tbl, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)
