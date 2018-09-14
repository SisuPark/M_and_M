
library(parallel)
library(rstan)
library(tools) # for file_path_as_absolute()

source("/Users/pgajer/organizer/programming/R/libs/stan_utils.R") # from https://github.com/betanalpha/knitr_case_studies/tree/master/rstan_workflow

setwd("/Users/pgajer/projects/M_and_M/data")
source("volcano_plot.R")
source("sel_phs.R")

load("mm_no_Thermodesulfatator_v3.rda")
load("mad_log10RA_vs_medianlog10nReads_plot_exp_model.rda") # m.dat, m.fit, nexp.m, nexpFn


## -------------
##    all
## -------------

idx1 <- !is.na(mt$bLoad) & mt$delStatus.n!=2
s <- sel.phs(idx1)
all.selPhs <- s$selPhs
length(all.selPhs) #
s$n.subj.TERM  #
s$n.subj.sPTB  #
s$n.TERM  #
s$n.sPTB  #
s$min.n.case.subjs #

all.picsDir <- "../pics/aa_all_dir/"
dir.create(all.picsDir)

system.time( all.res <- mclapply(all.selPhs, mc.cores=7, aa.spmrf.bernoulli.o2.fn, picsDir=all.picsDir) )

all.gEff <- matrix(nrow=length(all.selPhs), ncol=8)
colnames(all.gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","-mode($\\log_{10}$(AA))")
rownames(all.gEff) <- all.selPhs
for ( i in seq(all.selPhs) )
{
    ph <- all.selPhs[i]
    r <- all.res[[i]]
    all.gEff[ph,] <- r$gEff
    all.gEff[ph,2] <- r$r2$pval.min
}

all.gEff[,3] <- p.adjust(all.gEff[,2], method="fdr")

o <- order(all.gEff[,2])
all.gEff <- all.gEff[o,]

q.thld <- 0.05
e.thld <- 0.1
idx <- all.gEff[,3] < q.thld & abs(all.gEff[,1]) > e.thld
(all.gEff.sig <- as.data.frame(all.gEff)[idx,])
nrow(all.gEff.sig)

if ( nrow(all.gEff.sig) ){
    tblTexFile <- "../docs/sPTB_aa_spmrfo2/all_sPTB_aa_sig_gEff.tex"
    tblLabel <- "sPTB.aa.all:tbl"
    lt <- latex(format.gEff.tbl(all.gEff.sig), file=tblTexFile, label=tblLabel,
               caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
               where='!htp',vbar=FALSE)
}

pdf("../pics/sPTB_aa_spmrfo2/all_gEff_volcano_plot.pdf", width=6, height=6)
volcano.plot(all.gEff, q.thld=0.05) #, xlim=c(-1.1, 1.8))
dev.off()

## creating symbolic links to significant phylotypes
i <- 1
for ( ph in rownames(all.gEff.sig) )
{
    from <- file_path_as_absolute(paste0(all.picsDir, ph, ".pdf"))
    to <- paste0(all.picsDir, i, ".pdf")
    unlink(to)
    file.symlink(from, to)
    i <- i + 1
}
