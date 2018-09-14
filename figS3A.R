##
## bacterial load vs GA in TERM and sPTB samples
##

library(rms)

load("mm_no_mPTB_june25_2017.rda") # ct.no.mPTB, pt.no.mPTB, mt.no.mPTB, selPh.25

y <- log10(mt.no.mPTB$bLoad)
idx <- is.finite(y)
y <- y[idx]
subjID <- mt.no.mPTB$subjID[idx]
gr <- mt.no.mPTB$sPTB[idx]
gr.char <- ifelse(gr==1, "sPTB","TERM")
gr.f <- factor(gr.char, levels=c("TERM", "sPTB"))
gr.n <- as.integer(gr.f)-1
gr.n1 <- as.integer(gr.f)


ddist <- datadist(x, y, gr, gr.f, gr.char)
options("datadist" = "ddist")

m.ols3 <- ols( y ~ gr.f*rcs(x,3) ) #  penalty= list(simple =if( penalize ==1) pens[i] else 0, nonlinear= pens[i]))
summary(m.ols3)

## "In order to avoid overfitting we will select a model based on the AIC/BIC
## criterion. The selection is simply finding the lowest value where in general AIC
## allows slightly more complex models compared to BIC."
##
## from
## https://www.r-bloggers.com/an-exercise-in-non-linear-modeling/

source("getInfCrit.R")

## getInfCrit(4, m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=AIC)

## AIC
sapply(3:7, getInfCrit, fit=m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=AIC)
## gr*rcs(x, 3) gr*rcs(x, 4) gr*rcs(x, 5) gr*rcs(x, 6) gr*rcs(x, 7)
##     2733.649     2736.779     2735.113     2738.008     2733.824

## BIC
sapply(3:7, getInfCrit, fit=m.ols3, rm_var="gr.f*rcs(x,3)", add_var_str="gr*rcs(x, %d)", ic_fn=BIC)
## gr*rcs(x, 3) gr*rcs(x, 4) gr*rcs(x, 5) gr*rcs(x, 6) gr*rcs(x, 7)
##     2785.315     2800.147     2806.040     2821.257     2824.943

## so number of knots = 3 fits the data best


m.ols3.pred.TERM <- predict(m.ols3, data.frame(x=x[gr.char=="TERM"], gr.f=gr.f[gr.char=="TERM"]), conf.int=0.95)
m.ols3.pred.sPTB <- predict(m.ols3, data.frame(x=x[gr.char=="sPTB"], gr.f=gr.f[gr.char=="sPTB"]), conf.int=0.95)

r <- range(x[gr.char=="TERM"])

pdf("../pics/bacterial_load_vs_GA_ols3_v1_all.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="log10( Bacterial Load )", xlab="Gestation Age [weeks]")
points(jitter(x), y, col=ifelse(gr.char=="TERM","green","red"),pch=ifelse(gr.char=="TERM",20,19), cex=0.7)
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("bottomright", legend=c("TERM","sPTB"), pch=19, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()

pdf("../pics/bacterial_load_vs_GA_ols3_all.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, ylim=c(6.5, 8.5), pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="log10( Bacterial Load )", xlab="Gestation Age [weeks]")
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("bottomleft", legend=c("TERM","sPTB"), lty=1, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()


## AA

x <- mt.AA$visitGA
y <- log10(mt.AA$bLoad)
idx <- is.finite(y) & is.finite(x)
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

ddist <- datadist(x, y, gr, gr.f, gr.char)
options("datadist" = "ddist")

m.ols3 <- ols( y ~ gr.f*rcs(x,3) ) #  penalty= list(simple =if( penalize ==1) pens[i] else 0, nonlinear= pens[i]))
summary(m.ols3)

m.ols3.pred.TERM <- predict(m.ols3, data.frame(x=x[gr.char=="TERM"], gr.f=gr.f[gr.char=="TERM"]), conf.int=0.95)
m.ols3.pred.sPTB <- predict(m.ols3, data.frame(x=x[gr.char=="sPTB"], gr.f=gr.f[gr.char=="sPTB"]), conf.int=0.95)

r <- range(x[gr.char=="TERM"])

pdf("../pics/bacterial_load_vs_GA_ols3_v1_AA.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab="log10( Bacterial Load )", xlab="Gestation Age [weeks]")
points(jitter(x), y, col=ifelse(gr.char=="TERM","green","red"),pch=ifelse(gr.char=="TERM",20,19), cex=0.7)
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("bottomright", legend=c("TERM","sPTB"), pch=19, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()

pdf("../pics/bacterial_load_vs_GA_ols3_AA.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(x, y, las=1, xlim=r, ylim=c(6.5, 8.5), pch=19, col=ifelse(gr.char=="TERM","green","red"), type='n', ylab=expression(paste(log[10], "( Bacterial Absolute Abundance )")), xlab="Gestation Age [weeks]")
polygon(c(x[gr.char=="sPTB"], rev(x[gr.char=="sPTB"])), c(m.ols3.pred.sPTB$lower, rev(m.ols3.pred.sPTB$upper)), col='gray90', border=NA)
lines(x[gr.char=="sPTB"], m.ols3.pred.sPTB$linear.predictors, col='red',lwd=2)
polygon(c(x[gr.char=="TERM"], rev(x[gr.char=="TERM"])), c(m.ols3.pred.TERM$lower, rev(m.ols3.pred.TERM$upper)), col='gray80', border=NA)
lines(x[gr.char=="TERM"], m.ols3.pred.TERM$linear.predictors, col='green',lwd=2)
legend("bottomleft", legend=c("TERM","sPTB"), lty=1, col=c("green","red"), inset=0.05, cex=0.9)
par(op)
dev.off()
