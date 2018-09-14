

## mt.v1, ct.v1, bdt.v1,
## mt.v1.AA, ct.v1.AA, bdt.v1.AA,
## mt.v1.nonAA, ct.v1.nonAA, bdt.v1.nonAA,
load("bd_v1.rda")


## Is the mean log10 abundance of $\beta$-defensin in AA at visit V1 and
## specific CST samples significantly different between TERM and sPTB subjects?

mean.log10.bd.AA.v1.cst <- function(cst)
{
    idx <- mt.v1.AA$cst==cst
    y <- log10(bdt.v1.AA$BD_v1[idx])
    x <- mt.v1.AA$sPTB[idx]
    r <- t.test( y[x==0], y[x==1] )
    s <- c(paste0(sprintf("%.2f",r$estimate[1]), " (", sum(x==0), ")"),
          paste0(sprintf("%.2f",r$estimate[2]), " (", sum(x==1), ")"),
          sprintf("%.2f",r$p.value))
    s
}

cst <- "I"
(r <- mean.log10.bd.AA.v1.cst(cst))

mean.log10.bd.AA.v1.cst.tbl <- matrix(NA,nrow=length(names(table(mt.v1$cst))), ncol=3)
colnames(mean.log10.bd.AA.v1.cst.tbl) <- c("TERM (n)", "sPTB (n)", "p-val")
rownames(mean.log10.bd.AA.v1.cst.tbl) <- names(table(mt.v1$cst))
for ( cst in rownames(mean.log10.bd.AA.v1.cst.tbl) )
{
    print(cst)
    r <- mean.log10.bd.AA.v1.cst(cst)
    mean.log10.bd.AA.v1.cst.tbl[cst,] <- r
}

mean.log10.bd.AA.v1.cst.tbl


tblTexFile <- "../docs/Tables/mean_log10_bd_AA_V1_CST.tex"
tblLabel <- "mean.log10.bd.AA.v1.cst:tbl"
lt <- latex(mean.log10.bd.AA.v1.cst.tbl, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



##
## plotting mean beta-defensin values in TERM and sPTB AA subjects
##

pdf("../pics/mean_log10_bd_V1_AA.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(c(1-0.5,2+0.75), range(log10(bdt.v1.AA$BD_v1)), type='n', xlab="", ylab="log10 ( BD )", axes=F)
box()
axis(2, las=1)
axis(1, at=1, labels=c("TERM"), las=2)
axis(1, at=2, labels=c("sPTB"), las=2)
a <- 0.1
c <- 1.5
lwd <- 3
y <- log10(bdt.v1.AA$BD_v1)
x <- mt.v1.AA$sPTB
dx <- 0.05
r <- t.test( y[x==0], y[x==1] )
m.TERM <- r$estimate[1]
m.sPTB <- r$estimate[2]
points(jitter(rep(1, length(y[x==0])), amount=a), y[x==0])
points(jitter(rep(2, length(y[x==1])), amount=a), y[x==1])
segments(1-c*a, m.TERM, 1+c*a, m.TERM, col='red', lwd=lwd)
segments(2-c*a, m.sPTB, 2+c*a, m.sPTB, col='red', lwd=lwd)
## significance bracket
segments(1+c*a+dx, m.TERM, 2.5, m.TERM, col=1, lwd=0.5)
segments(2+c*a+dx, m.sPTB, 2.5, m.sPTB, col=1, lwd=0.5)
segments(2.5, m.TERM, 2.5, m.sPTB, col=1, lwd=0.5)
text(2.5+dx, (m.TERM + m.sPTB)/2, labels="*", cex=1.5)
par(op)
r
##t = 3.7439, df = 117.25, p-value = 0.0002823
## mean of x mean of y
##  4.476802  4.100683
dev.off()
