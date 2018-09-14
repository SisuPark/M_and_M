##
##
##

## mt.v1, ct.v1, bdt.v1,
## mt.v1.AA, ct.v1.AA, bdt.v1.AA,
## mt.v1.nonAA, ct.v1.nonAA, bdt.v1.nonAA,
load("bd_v1.rda")

load("mm_no_Thermodesulfatator_v3.rda") # mt, ct, pt, selPh2.10

library(rstan)
source("stan_utils.R")

phs <- c("Lactobacillus_iners",
        "Lactobacillus_crispatus",
        "Lactobacillus_jensenii",
        "Lactobacillus_gasseri",
        "Mobiluncus_curtisii_mulieris",
        "Sneathia_sanguinegens",
        "g_Megasphaera",
        "BVAB3",
        "Porphyromonas_asaccharolytica",
        "g_Atopobium",
        "Prevotella_buccalis")


v1.ph.fn <- function(ph, nItr=50000, thin=5, nChains=4, nCores=4)
{
    idx <- ct.v1[, ph] > 0
    y <- log10(bdt.v1$BD_v1[idx])
    x <- mt.v1$sPTB[idx]
    idx <- !is.na(y) & !is.na(x)
    y <- y[idx]
    x <- x[idx]
    x <- ifelse(x==0,"TERM","sPTB")
    x <- factor(x, levels=c("sPTB","TERM"))
    qq <- get.qq(y)

    (f <- table(qq,x))
    (p <- 100*f[,1]/rowSums(f)) ## prop.test(f)
    (n <- rowSums(f))

    m.dat <- list(y=as.vector(f[c(1,4),1]),
                 n=as.vector(n[c(1,4)])) ## str(m.dat)

    r <- m.2prop.fn(m.dat, nItr=nItr, thin=thin, nChains=nChains, nCores=nCores) r$logRat.ci r$p.ci ## qqPlot(r$logRat)
    s <- c(paste0(f[1,1], "/", sum(f[1,]), " (", sprintf("%.1f",p[1]),")"), paste0(f[4,1], "/", sum(f[4,]), " (", sprintf("%.1f",p[4]),")"), sprintf("%.4f",r$logRat.ci[4]))
}


bd.v1.ph.tbl <- matrix(NA,nrow=length(phs), ncol=3)
colnames(bd.v1.ph.tbl) <- c("Q1 sPTB/TERM (\\%)", "Q4 sPTB/TERM (\\%)", "p-val")
rownames(bd.v1.ph.tbl) <- phs
for ( ph in phs )
{
    r <- v1.ph.fn(ph, nCores=1)
    bd.v1.ph.tbl[ph,] <- r
}
bd.v1.ph.tbl

q <- p.adjust(as.numeric(bd.v1.ph.tbl[,3]), method="fdr")
bd.v1.ph.tbl <- cbind(bd.v1.ph.tbl, sprintf("%.5f",q))
colnames(bd.v1.ph.tbl)[4] <- "q-val"


tblTexFile <- "../docs/Tables/q_bd_V1_ph.tex"
tblLabel <- "q.bd.v1.ph:tbl"
lt <- latex(bd.v1.ph.tbl, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)
