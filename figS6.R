

## mt.v1, ct.v1, bdt.v1,
## mt.v1.AA, ct.v1.AA, bdt.v1.AA,
## mt.v1.nonAA, ct.v1.nonAA, bdt.v1.nonAA,
load("bd_v1.rda")


## Is the mean log10 abundance of $\beta$-defensin in all, AA and nonAA subjects
## at visit V1

## All
idx <- mt.v1$cst==cst
y <- log10(bdt.v1$BD_v1[idx])
x <- mt.v1$sPTB[idx]
r <- t.test( y[x==0], y[x==1] )
s <- c(paste0(sprintf("%.2f",r$estimate[1]), " (", sum(x==0), ")"),
      paste0(sprintf("%.2f",r$estimate[2]), " (", sum(x==1), ")"),
      sprintf("%.2f",r$p.value))


## AA
idx <- mt.v1.AA$cst==cst
y <- log10(bdt.v1.AA$BD_v1[idx])
x <- mt.v1.AA$sPTB[idx]
r <- t.test( y[x==0], y[x==1] )
s <- c(paste0(sprintf("%.2f",r$estimate[1]), " (", sum(x==0), ")"),
      paste0(sprintf("%.2f",r$estimate[2]), " (", sum(x==1), ")"),
      sprintf("%.2f",r$p.value))


## nonAA
idx <- mt.v1.nonAA$cst==cst
y <- log10(bdt.v1.nonAA$BD_v1[idx])
x <- mt.v1.nonAA$sPTB[idx]
r <- t.test( y[x==0], y[x==1] )
s <- c(paste0(sprintf("%.2f",r$estimate[1]), " (", sum(x==0), ")"),
      paste0(sprintf("%.2f",r$estimate[2]), " (", sum(x==1), ")"),
      sprintf("%.2f",r$p.value))





##
## plotting mean beta-defensin values in TERM and sPTB in all, AA and nonAA
## subjects at V1
##

pdf("../pics/mean_log10_bd_V1.pdf", width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot(c(1-0.5,6+0.75), range(log10(c(bdt.v1$BD_v1, bdt.v1.AA$BD_v1, bdt.v1.nonAA$BD_v1))), type='n', xlab="", ylab="log10 ( BD )", axes=F)
box()
axis(2, las=1)
axis(1, at=1, labels=c("TERM"), las=2)
axis(1, at=2, labels=c("sPTB"), las=2)
axis(1, at=3, labels=c("TERM"), las=2)
axis(1, at=4, labels=c("sPTB"), las=2)
axis(1, at=5, labels=c("TERM"), las=2)
axis(1, at=6, labels=c("sPTB"), las=2)
a <- 0.1
c <- 1.5
lwd <- 3
## All
y <- log10(bdt.v1$BD_v1)
x <- mt.v1$sPTB
dx <- 0.05
r <- t.test( y[x==0], y[x==1] )
m.TERM <- r$estimate[1]
m.sPTB <- r$estimate[2]
points(jitter(rep(1, length(y[x==0])), amount=a), y[x==0])
points(jitter(rep(2, length(y[x==1])), amount=a), y[x==1])
segments(1-c*a, m.TERM, 1+c*a, m.TERM, col='red', lwd=lwd)
segments(2-c*a, m.sPTB, 2+c*a, m.sPTB, col='red', lwd=lwd)
## AA
y <- log10(bdt.v1.AA$BD_v1)
x <- mt.v1.AA$sPTB
dx <- 0.05
r <- t.test( y[x==0], y[x==1] )
m.TERM <- r$estimate[1]
m.sPTB <- r$estimate[2]
points(jitter(rep(3, length(y[x==0])), amount=a), y[x==0])
points(jitter(rep(4, length(y[x==1])), amount=a), y[x==1])
segments(3-c*a, m.TERM, 3+c*a, m.TERM, col='red', lwd=lwd)
segments(4-c*a, m.sPTB, 4+c*a, m.sPTB, col='red', lwd=lwd)
## AA
y <- log10(bdt.v1.AA$BD_v1)
x <- mt.v1.AA$sPTB
dx <- 0.05
r <- t.test( y[x==0], y[x==1] )
m.TERM <- r$estimate[1]
m.sPTB <- r$estimate[2]
points(jitter(rep(5, length(y[x==0])), amount=a), y[x==0])
points(jitter(rep(6, length(y[x==1])), amount=a), y[x==1])
segments(5-c*a, m.TERM, 5+c*a, m.TERM, col='red', lwd=lwd)
segments(6-c*a, m.sPTB, 6+c*a, m.sPTB, col='red', lwd=lwd)
par(op)
dev.off()
