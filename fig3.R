
## mt.v1, ct.v1, bdt.v1,
## mt.v1.AA, ct.v1.AA, bdt.v1.AA,
## mt.v1.nonAA, ct.v1.nonAA, bdt.v1.nonAA,
load("bd_v1.rda")

load("mm_no_Thermodesulfatator_v3.rda") # mt, ct, pt, selPh2.10

##
## Mb :: Lb :: BD interaction in the context of sPTB in AA at visit V1
##

ph <- "Mobiluncus_curtisii_mulieris"
idx <- ct.v1.AA[, ph] > 0
y <- mt.v1.AA$sPTB[idx]
## colnames(ct)[c(1,2,4,6)]
## [1] "Lactobacillus_iners"     "Lactobacillus_crispatus" "Lactobacillus_jensenii"
## [4] "Lactobacillus_gasseri"
x <- rowSums(ct.v1.AA[idx, c(1,2,4,6)])
n <- rowSums(ct.v1.AA[idx,])
x.Lb <- log10(x/n)
x.Mb <- log10(pt.v1.AA[idx,ph])
x.bd <- log10(bdt.v1.AA$BD_v1[idx])
x.cst <- mt.v1.AA$cst[idx]
x.subjID <- mt.v1.AA$subjID[idx]
idx <- is.finite(x.Lb) & is.finite(x.bd)
x.Lb <- x.Lb[idx]
x.Mb <- x.Mb[idx]
x.bd <- x.bd[idx]
x.cst <- x.cst[idx]
x.subjID <- x.subjID[idx]
y <- y[idx]
o <- order(x.Mb)
y <- y[o]
x.Lb <- x.Lb[o]
x.Mb <- x.Mb[o]
x.bd <- x.bd[o]
x.cst <- x.cst[o]
x.subjID <- x.subjID[o]
length(y) # 166

idx <- x.Lb>-3
x.Lb <- x.Lb[idx]
x.Mb <- x.Mb[idx]
x.bd <- x.bd[idx]
x.cst <- x.cst[idx]
x.subjID <- x.subjID[idx]
y <- y[idx]

T.Lb <- quantile(x.Lb, probs = c(0.3333,0.6666))
T.bd <- quantile(x.bd, probs = c(0.3333,0.6666))

t.Lb <- x.Lb
t.Lb[x.Lb<=T.Lb[1]] <- "T1"
t.Lb[x.Lb<=T.Lb[2] & x.Lb>T.Lb[1]] <- "T2"
t.Lb[x.Lb>T.Lb[2]] <- "T3"

t.bd <- x.bd
t.bd[x.bd<=T.bd[1]] <- "T1"
t.bd[x.bd<=T.bd[2] & x.bd>T.bd[1]] <- "T2"
t.bd[x.bd>T.bd[2]] <- "T3"


f.LbBD <- list()
p.LbBD <- list()
for ( i in 1:3 )
{
    tLbStr <- paste0("T",i)
    for ( j in 1:3 )
    {
        tBDStr <- paste0("T",j)
        idx <- t.Lb==tLbStr & t.bd==tBDStr
        yT <- y[idx]
        f <- table(yT)
        p <- 100*f/sum(f)
        id <- paste0(tLbStr,"::",tBDStr)
        f.LbBD[[id]] <- f
        p.LbBD[[id]] <- p
    }
}

myHist(x.Mb)

(f <- table(y))
##   0   1
## 125  40

(p <- 100*f/sum(f))
##  0        1
## 75.75758 24.24242


file <- "../pics/Lb_vs_BD_vs_Mb_in_AA_with_CSTs_with_no_CSTI_outlier.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(3.75, 3.75, 4.5, 3.5), mgp=c(2.5,0.5,0),tcl = -0.3)
plot(x.Lb, x.bd, xlab="log10( relAb of Lactobacillus )",
     ylab="log10( BD )",
     cex=ifelse(x.Mb> -2, 2, 1), type='n',
     las=1 )
points(x.Lb, x.bd, pch=ifelse(y==1, 17, 19), cex=ifelse(x.Mb> -2, 1.7, 0.7), col=cstColTbl[x.cst])
cf <- 1.1
points(x.Lb, x.bd, pch=ifelse(y==1, 2, 22), cex=ifelse(x.Mb> -2, cf*1.7, cf*1.1), col=1)#cstColTbl[x.cst])
abline(v=T.Lb[1], lwd=1)
abline(v=T.Lb[2], lwd=1)
mtext("Lb 33%", side=3, at=T.Lb[1], line=0.25)
mtext("Lb 66%", side=3, at=T.Lb[2], line=0.25)
abline(h=T.bd[1], lwd=1)
abline(h=T.bd[2], lwd=1)
mtext("bd 33%", side=4, at=T.bd[1], las=1, line=0.25)
mtext("bd 66%", side=4, at=T.bd[2], las=1, line=0.25)
legend(par('usr')[1], par('usr')[4], xjust=-6, yjust=-0.5 , xpd=NA, cex=0.8, legend=c("sPTB", "TERM"), pch=c(17,19), ncol=1, pt.cex=1.25)
 legend(par('usr')[1], par('usr')[4], xjust=0, yjust=-1.25 , xpd=NA, cex=0.8, legend=names(cstColTbl), ncol=length(cstColTbl), fill=cstColTbl)  #title="Community State Type", cex=1
par(op)
dev.off()



file <- "../pics/Lb_vs_BD_vs_Mb_in_AA_with_CSTs_with_no_CSTI_outlier_version_with_squares.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(3.75, 3.75, 4.5, 3.5), mgp=c(2.5,0.5,0),tcl = -0.3)
plot(x.Lb, x.bd, xlab="log10( relAb of Lactobacillus )",
     ylab="log10( BD )",
     cex=ifelse(x.Mb> -2, 2, 1), type='n',
     las=1 )
points(x.Lb, x.bd, pch=ifelse(y==1, 17, 15), cex=ifelse(x.Mb> -2, 1.7, 1.1), col=cstColTbl[x.cst])
cf <- 1.1
points(x.Lb, x.bd, pch=ifelse(y==1, 2, NA), cex=ifelse(x.Mb> -2, cf*1.7, cf*1.1), col=1)#cstColTbl[x.cst])
abline(v=T.Lb[1], lwd=1)
abline(v=T.Lb[2], lwd=1)
mtext("Lb 33%", side=3, at=T.Lb[1], line=0.25)
mtext("Lb 66%", side=3, at=T.Lb[2], line=0.25)
abline(h=T.bd[1], lwd=1)
abline(h=T.bd[2], lwd=1)
mtext("bd 33%", side=4, at=T.bd[1], las=1, line=0.25)
mtext("bd 66%", side=4, at=T.bd[2], las=1, line=0.25)
legend(par('usr')[1], par('usr')[4], xjust=-6, yjust=-0.55 , xpd=NA, cex=0.8, legend=c("sPTB", "TERM"), pch=c(17,15), ncol=1, pt.cex=1.25)
legend(par('usr')[1], par('usr')[4], xjust=0, yjust=-1.35 , xpd=NA, cex=0.65, legend=names(cstColTbl), ncol=length(cstColTbl), fill=cstColTbl)  #title="Community State Type", cex=1
par(op)
dev.off()
