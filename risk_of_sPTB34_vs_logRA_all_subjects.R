##
## risk of sPTB34 vs log10 relAb's in all subjects over all visits
##

library(parallel)
library(rstan)

source("stan_utils.R")
source("volcano_plot.R")

load("mad_log10RA_vs_medianlog10nReads_plot_exp_model.rda") # m.dat, m.fit, nexp.m, nexpFn
load("mm_no_Thermodesulfatator_v3.rda") # ct, pt, selPh.10, selPh.25

## spmrf.bernoulli.o2.model <- stan_model(file="/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order2.stan")
## save(spmrf.bernoulli.o2.model, file="/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order2.rda")
load("/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order2.rda")

y <- mt$delStatus34
s <- mt$subjID
idx1 <- mt$delStatus34!=2
y <- y[idx1]
s <- s[idx1]
(tt <- table(y))

n.TERM <- tt[1][[1]]
n.sPTB <- tt[2][[1]]

(n.subj.TERM <- length(unique(s[y==0])))
(n.subj.sPTB <- length(unique(s[y==1])))

thld <- 0.25
(n.TERM.25 <- as.integer(n.TERM*thld)) #
(n.sPTB.25 <- as.integer(n.sPTB*thld)) #

sel.25 <- c()
sel.25.list <- list()
for ( ph in selPh.25 )
{
    x <- pt[,ph]
    y <- mt$delStatus34
    idx1 <- mt$delStatus34!=2
    x <- x[idx1]
    y <- y[idx1]
    idx <- x > 0
    x <- log10(x[idx])
    y <- y[idx]
    q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x > q[1] & x < q[2]
    x <- x[idx]
    y <- y[idx]
    if ( length(x[y==0]) > n.TERM.25 && length(x[y==1]) > n.sPTB.25 )
    {
        sel.25 <- c(sel.25, ph)
        sel.25.list[[ph]] <- c(length(x[y==0]), length(x[y==1]))
    }
}

length(sel.25) # 43


##
##
##

pics.dir <- "../pics/sPTB43_dir/"
dir.create(pics.dir)

pr.sPTB34.all.subjects.fn <- function(ph, nModels=5, nChains=3, nItr=2000, thin=1, nCores=3, showMessages=FALSE, picsDir=pics.dir)
{
    x <- ct[,ph]
    n <- rowSums(ct)
    y <- mt$delStatus34
    s <- mt$subjID
    idx1 <- mt$delStatus34!=2
    x <- x[idx1]
    y <- y[idx1]
    n <- n[idx1]
    s <- s[idx1]
    idx <- x > 0
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    s <- s[idx]
    m <- length(x)
    n.case.subjs <- length(unique(s[y==1]))
    n.case.samples <- length(s[y==1])
    n.ctrk.subjs <- length(unique(s[y==0]))
    n.ctrk.samples <- length(s[y==0])
    q <- quantile(x/n, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x/n > q[1] & x/n < q[2]
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    s <- s[idx]
    x <- log10(x/n)
    o <- order(x)
    x <- x[o]
    y <- y[o]

    r <- spmrf.bernoulli.o2.fn.xy(x, y, nModels=nModels, nItr=nItr, thin=thin, nChains=nChains, nCores=nCores, showMessages=showMessages)

    idx <- seq(r$y.med)

    file <- paste0(picsDir, ph, ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",m,")")
    plot.logit(x[idx], y[idx], expit(r$y.med), expit(r$y.l), expit(r$y.u), title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    par(op)
    dev.off()

    r$gEff["n(TERM)"] <- n.ctrk.samples
    r$gEff["n(subj TERM)"] <- n.ctrk.subjs
    r$gEff["n(sPTB)"] <- n.case.samples
    r$gEff["n(subj sPTB)"] <- n.case.subjs
    dx <- density(x)
    dx.mode <- denMode(dx)
    r$gEff["mode($\\log_{10}$(rat))"] <- -dx.mode

    list(ph=ph,
         file=file,
         r=r,
         x=x,
         y=y,
         m.dat=m.dat,
         m=m,
         m.dat=m.dat,
         n.case.subjs=n.case.subjs,
         n.case.samples=n.case.samples,
         n.ctrk.subjs=n.ctrk.subjs,
         n.ctrk.samples=n.ctrk.samples,
         minus.median.log10.RA=-median(x),
         gEff=r$gEff)
}


all34.res <- mclapply(sel.25, pr.sPTB34.all.subjects.fn, mc.cores=7)


all34.gEff <- matrix(nrow=length(sel.25), ncol=8)
colnames(all34.gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","-mode($\\log_{10}$(RA))")
rownames(all34.gEff) <- sel.25
for ( i in seq(sel.25) )
{
    r <- all34.res[[i]]
    ph <- r$ph
    all34.gEff[ph,] <- r$gEff
}

all34.gEff[,3] <- p.adjust(all34.gEff[,2], method="fdr")

o <- order(all34.gEff[,2])
all34.gEff <- all34.gEff[o,]

q.thld <- 0.05
e.thld <- 0.1
idx <- all34.gEff[,3] < q.thld & abs(all34.gEff[,1]) > e.thld
(all34.gEff.sig <- as.data.frame(all34.gEff)[idx,])
