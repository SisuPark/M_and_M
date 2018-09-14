##
## Kaplan-Meier survival models/curves
##

library(survival)

source("stan_utils.R")

load("mm_no_mPTB_june25_2017.rda")     # ct.no.mPTB, pt.no.mPTB, mt.no.mPTB, selPh.25
load("spmrf_bernoulli_o2_all_res.rda") # all.res2

##
## all samples; all subjects
##

outDir10 <- "/Users/pgajer/projects/M_and_M/pics/all_10perc/"
unlink(outDir10, recursive = T)
dir.create(outDir10)

idx <- mt$delStatus.n != 2
mt.no.mPTB <- mt[idx,]
ct.no.mPTB <- ct2[idx,]
pt.no.mPTB <- pt2[idx,]

phs <- c("Sneathia_sanguinegens", "Mobiluncus_curtisii_mulieris", "g_Megasphaera", "BVAB3", "g_Atopobium", "Porphyromonas_asaccharolytica", "Prevotella_buccalis")

pval <- c()
sig.thld <- c()
n.gr <- c()
n.less <- c()
hr <- c()
for ( ph in phs )
{
    i <- grep(ph, selPhs)
    r <- all.res2[[i]]
    ## 10% section
    sig.thld <- sig.thld.fn3(r,p=0.1)
    sig.thld[ph] <- sig.thld
    idx <- ct.no.mPTB[,ph]>0
    X <- mt.no.mPTB[idx,]
    nrow(X)
    x <- log10(pt.no.mPTB[idx,ph])
    X$is.above.thld <- ifelse(x>=sig.thld,"Yes","No")
    m <- survfit(Surv(delGA, rep(1,nrow(X))) ~ is.above.thld, data = X)
    file <- paste0(outDir10,ph,"_all_samples_p10_survival.pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    plot(m, xlim=c(15,43), col=1:2, las=1, xlab="Gestation Age [weeks]", main=ph)
    thld <- formatC(10^sig.thld,digits=3)
    legend("bottomleft",legend=c(paste0("RA < ",thld), paste0("RA >= ",thld)), col=1:2, lwd=1,inset=0.05)
    par(op)
    dev.off()
    m <- coxph(Surv(delGA, rep(1,nrow(X))) ~ is.above.thld, data = X)
    s <- summary(m)
    pval[ph] <- s$logtest[3]
    xx <- table(X$is.above.thld)
    n.gr[ph] <- xx[2][[1]]
    n.less[ph] <- xx[1][[1]]
    hr[ph] <- s$coefficients[2]
}

x <- cbind(sig.thld, 10^sig.thld,n.less,n.gr, hr, pval)
qval <- p.adjust(pval, method="fdr")
x <- cbind(x, qval)
colnames(x) <- c("log10 thld", "thld"," n.less","n.gr", "HR", "p-val", "q-val")
x

for ( ph in phs )
{
    i <- grep(ph, selPhs)
    r <- all.res2[[i]]
    sig.thld <- sig.thld.fn3(r,p=0.1)
    file <- paste0(outDir10,ph,"_risk_of_sPTB_vs_relAb_all_visits_all_subjects_p10.pdf")
    pdf(file, width=6, height=6)
    plot.logit.incr(ph, r, sig.thld)
    dev.off()
}



##
## all visits; African-American subjects
##

load("spmrf_bernoulli_o2_AA_res.rda") # AA.res

idx1 <- mt$delStatus.n!=2 & !is.na(mt$race) & mt$raceChar=="black"
s <- sel.phs(idx1)
AA.selPhs <- s$selPhs
length(AA.selPhs) # 74

## AA.selPhs <- c()
## for ( i in seq(AA.res) )
## {
##     r <- AA.res[[i]]
##     AA.selPhs[i] <- r$ph
## }
## length(AA.selPhs)


##(AA.sig.phs <- rownames(AA.gEff.sig))
AA.sig.phs <- c("Mobiluncus_curtisii_mulieris",
               "BVAB3",
               "Sneathia_sanguinegens",
               "Porphyromonas_asaccharolytica",
               "g_Megasphaera")


outDir10 <- "/Users/pgajer/projects/M_and_M/pics/AA_10perc/"
unlink(outDir10, recursive = T)
dir.create(outDir10)

idx <- mt$delStatus.n != 2 & !is.na(mt$raceChar) & mt$raceChar=="black"
mt.AA <- mt[idx,]
ct.AA <- ct[idx,]
pt.AA <- ct.AA/rowSums(ct.AA)


pval <- c()
sig.thld <- c()
n.gr <- c()
n.less <- c()
hr <- c()
for ( ph in AA.sig.phs )
{
    i <- grep(ph, AA.selPhs)
    r <- AA.res[[i]]
    sig.thld <- sig.thld.fn3(r,p=0.1)
    sig.thld[ph] <- sig.thld
    idx <- ct.AA[,ph]>0
    X <- mt.AA[idx,]
    nrow(X)
    x <- log10(pt.AA[idx,ph])
    X$is.above.thld <- ifelse(x>=sig.thld,"Yes","No")
    m <- survfit(Surv(delGA, rep(1,nrow(X))) ~ is.above.thld, data = X)
    file <- paste0(outDir10,ph,"_AA_p10_survival.pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    plot(m, xlim=c(15,43), col=1:2, las=1, xlab="Gestation Age [weeks]", main=ph)
    thld <- formatC(10^sig.thld,digits=3)
    legend("bottomleft",legend=c(paste0("RA < ",thld), paste0("RA >= ",thld)), col=1:2, lwd=1,inset=0.05)
    par(op)
    dev.off()
    m <- coxph(Surv(delGA, rep(1,nrow(X))) ~ is.above.thld, data = X)
    s <- summary(m)
    pval[ph] <- s$logtest[3]
    xx <- table(X$is.above.thld)
    n.gr[ph] <- xx[2][[1]]
    n.less[ph] <- xx[1][[1]]
    hr[ph] <- s$coefficients[2]
}

x <- cbind(sig.thld, 10^sig.thld,n.less,n.gr, hr, pval)
qval <- p.adjust(pval, method="fdr")
x <- cbind(x, qval)
colnames(x) <- c("log10 thld", "thld"," n.less","n.gr", "HR", "p-val", "q-val")
x

for ( ph in AA.sig.phs )
{
    i <- grep(ph, AA.selPhs)
    r <- AA.res[[i]]
    sig.thld <- sig.thld.fn3(r,p=0.1)
    file <- paste0(outDir10,ph,"_risk_of_sPTB_vs_relAb_AA_p10.pdf")
    pdf(file, width=6, height=6)
    plot.logit.incr(ph, r, sig.thld)
    dev.off()
}
