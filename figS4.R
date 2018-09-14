
##
## Shannon diversity vs GA
##

library(rms)

source("stan_utils.R")
source("volcano_plot.R")

load("mad_log10RA_vs_medianlog10nReads_plot_exp_model.rda") # m.dat, m.fit, nexp.m, nexpFn
load("mm_no_Thermodesulfatator_v3.rda") # ct, pt, selPh.10, selPh.25


##
## all subjects
##

Sdiv <- apply(pt, 1, entropy)

x <- mt$visitGA
y <- Sdiv
idx <- !is.na(mt$visitGA)
x <- x[idx]
y <- y[idx]
gr <- mt$sPTB[idx]
gr.char <- ifelse(gr==1, "sPTB","TERM")
gr.f <- factor(gr.char, levels=c("TERM", "sPTB"))
gr.n <- as.integer(gr.f)-1
gr.n1 <- as.integer(gr.f)
o <- order(x)
x <- x[o]
y <- y[o]
gr.f <- gr.f[o]
gr.char <- gr.char[o]


ddist <- datadist(x, y, gr, gr.f, gr.char)
options("datadist" = "ddist")

m.ols3 <- ols( y ~ gr.f*rcs(x,3) ) #  penalty= list(simple =if( penalize ==1) pens[i] else 0, nonlinear= pens[i]))
summary(m.ols3)

sapply(3:7, getInfCrit, fit=m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=AIC)

sapply(3:7, getInfCrit, fit=m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=BIC)

## so number of knots = 3 fits the data best


m.ols3.pred.TERM <- predict(m.ols3, data.frame(x=x[gr.char=="TERM"], gr.f=gr.f[gr.char=="TERM"]), conf.int=0.95)
m.ols3.pred.sPTB <- predict(m.ols3, data.frame(x=x[gr.char=="sPTB"], gr.f=gr.f[gr.char=="sPTB"]), conf.int=0.95)

r <- range(x[gr.char=="TERM"])

pdf("../pics/shannon_diversity_vs_GA_ols3_v1_all.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="log10( Shannon Diversity )", xlab="Gestation Age [weeks]")
points(jitter(x), y, col=ifelse(gr.char=="TERM","green","red"),pch=ifelse(gr.char=="TERM",20,19), cex=0.7)
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("topright", legend=c("TERM","sPTB"), pch=19, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()

pdf("../pics/shannon_diversity_vs_GA_ols3_v2_all.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, ylim=c(0.75, 2.25), pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="log10( Shannon Diversity )", xlab="Gestation Age [weeks]")
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("bottomleft", legend=c("TERM","sPTB"), lty=1, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()




##
## diversity for AA
##

idx <- mt$delStatus.n != 2 & !is.na(mt$raceChar) & mt$raceChar=="black"
mt.AA <- mt[idx,]
ct.AA <- ct[idx,]
pt.AA <- pt[idx,]

entropy <- function(x){
    x <- x[x>0]
    x <- x / sum(x)
    -sum( x * log2(x) )
}

Sdiv.AA <- apply(pt.AA, 1, entropy)

x <- mt.AA$visitGA
y <- Sdiv.AA
idx <- !is.na(mt.AA$visitGA)
x <- x[idx]
y <- y[idx]
gr <- mt.AA$sPTB[idx]
gr.char <- ifelse(gr==1, "sPTB","TERM")
gr.f <- factor(gr.char, levels=c("TERM", "sPTB"))
gr.n <- as.integer(gr.f)-1
gr.n1 <- as.integer(gr.f)
o <- order(x)
x <- x[o]
y <- y[o]
gr.f <- gr.f[o]
gr.char <- gr.char[o]

plot(jitter(x), y, las=1, pch=19, col=ifelse(gr.f=="TERM","green","red"), cex=ifelse(gr.f=="TERM",0.5,1)) # pch=ifelse(gr.f=="TERM",1,19)

ddist <- datadist(x, y, gr, gr.f, gr.char)
options("datadist" = "ddist")

m.ols3 <- ols( y ~ gr.f*rcs(x,3) ) #  penalty= list(simple =if( penalize ==1) pens[i] else 0, nonlinear= pens[i]))
summary(m.ols3)

sapply(3:7, getInfCrit, fit=m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=AIC)

sapply(3:7, getInfCrit, fit=m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=BIC)

## so number of knots = 3 fits the data best


m.ols3.pred.TERM <- predict(m.ols3, data.frame(x=x[gr.char=="TERM"], gr.f=gr.f[gr.char=="TERM"]), conf.int=0.95)
m.ols3.pred.sPTB <- predict(m.ols3, data.frame(x=x[gr.char=="sPTB"], gr.f=gr.f[gr.char=="sPTB"]), conf.int=0.95)

r <- range(x[gr.char=="TERM"])

pdf("../pics/shannon_diversity_vs_GA_ols3_v1_AA.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="log10( Shannon Diversity )", xlab="Gestation Age [weeks]")
points(jitter(x), y, col=ifelse(gr.char=="TERM","green","red"),pch=ifelse(gr.char=="TERM",20,19), cex=0.7)
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("topright", legend=c("TERM","sPTB"), pch=19, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()

pdf("../pics/shannon_diversity_vs_GA_ols3_AA.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, ylim=c(0.75, 2.25), pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="Shannon Diversity", xlab="Gestation Age [weeks]")
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("bottomleft", legend=c("TERM","sPTB"), lty=1, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()
