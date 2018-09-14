
library(rstan)
source("stan_utils.R")
load("mm_no_Thermodesulfatator_v3.rda") # ct, pt, selPh.10, selPh.25

## =======================
##   Lb in all non.mPTB
## =======================

r.Mb.LbT.list <- list()

for ( Tstr in c("T1","T2","T3") )
{
    print(Tstr)
    r.Mb.LbT <- spmrf.ph.T.fn(T=Tstr, phs=c("Lactobacillus_iners", "Lactobacillus_crispatus", "Lactobacillus_jensenii", "Lactobacillus_gasseri"), ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
    r.Mb.LbT.list[[Tstr]] <- r.Mb.LbT
    file <- paste0("../pics/risk_of_sPTB_vs_Mb_in_Lb",Tstr,".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0("Mb+ & Lb",Tstr," samples"," (n=",length(r.Mb.LbT$y),")")
    plot.spmrf.logit(r.Mb.LbT$m, r.Mb.LbT$x, r.Mb.LbT$y, title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    par(op)
    dev.off()
}
