
##
## Comparing CST proportions between sPTB and TERM subjects at visit 1 and over
## all visits
##

load("mm_dada2_ct_pt.rda") # ct, pt, mt


## ---------------------------------
##     All subjects - all visits
## ---------------------------------

idx <- mt$delStatus.n != 2
(cst.freq <- table(mt$delStatus[idx], mt$cst[idx]))
(cst.perc <- 100 * cst.freq / rowSums(cst.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 27.272727  7.438017 17.768595 26.859504 15.289256  5.371901
## TERM 31.512272 10.134600 19.081552 19.477435 13.618369  6.175772

table(mt$cst)
##    I   II  III IV-A IV-B    V
##  498  157  312  361  224   96


cst.pvals <- c()

##
## CST I
##
cstI <- ifelse(mt$cst[idx]=="I", 1, 0)
length(cstI) # 1505

mI.dat <- list(y=cstI,
              delSt=mt$delStatus[idx],
              subjID=factor(mt$subjID[idx]))

str(mI.dat)

m.I <- glmer( y ~ delSt + (1|subjID), data=mI.dat, family=poisson )
summary(m.I)

cst.pvals[1] <- coef(summary(m.I))[2,4]



##
## CST II
##
cstII <- ifelse(mt$cst[idx]=="II", 1, 0)
length(cstII) # 1505

mII.dat <- list(y=cstII,
              delSt=mt$delStatus[idx],
              subjID=factor(mt$subjID[idx]))

str(mII.dat)

m.II <- glmer( y ~ delSt + (1|subjID), data=mII.dat, family=poisson )
summary(m.II)

cst.pvals[2] <- coef(summary(m.II))[2,4]

##
## CST III
##
cstIII <- ifelse(mt$cst[idx]=="III", 1, 0)
length(cstIII) # 1505

mIII.dat <- list(y=cstIII,
              delSt=mt$delStatus[idx],
              subjID=factor(mt$subjID[idx]))

str(mIII.dat)

m.III <- glmer( y ~ delSt + (1|subjID), data=mIII.dat, family=poisson )
summary(m.III)

cst.pvals[3] <- coef(summary(m.III))[2,4]



##
## CST IV-A
##
cstIVA <- ifelse(mt$cst[idx]=="IV-A", 1, 0)
length(cstIVA) # 1505

mIVA.dat <- list(y=cstIVA,
              delSt=mt$delStatus[idx],
              subjID=factor(mt$subjID[idx]))

str(mIVA.dat)

m.IVA <- glmer( y ~ delSt + (1|subjID), data=mIVA.dat, family=poisson )
summary(m.IVA)

cst.pvals[4] <- coef(summary(m.IVA))[2,4]



##
## CST IV-B
##
cstIVB <- ifelse(mt$cst[idx]=="IV-B", 1, 0)
length(cstIVB) # 1505

mIVB.dat <- list(y=cstIVB,
              delSt=mt$delStatus[idx],
              subjID=factor(mt$subjID[idx]))

str(mIVB.dat)

m.IVB <- glmer( y ~ delSt + (1|subjID), data=mIVB.dat, family=poisson )
summary(m.IVB)

cst.pvals[5] <- coef(summary(m.IVB))[2,4]

##
## CST V
##
cstV <- ifelse(mt$cst[idx]=="V", 1, 0)
length(cstV) # 1505

mV.dat <- list(y=cstV,
              delSt=mt$delStatus[idx],
              subjID=factor(mt$subjID[idx]))

str(mV.dat)

m.V <- glmer( y ~ delSt + (1|subjID), data=mV.dat, family=poisson )
summary(m.V)

cst.pvals[6] <- coef(summary(m.V))[2,4]



## -----------------------------------
##     All subjects - Visit V1
## -----------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V1"
mt.V1 <- mt[idx,]

(cst.v1.freq <- table(mt.V1$delStatus, mt.V1$cst))
(cst.v1.perc <- 100 * cst.v1.freq / rowSums(cst.v1.freq))

file <- "../pics/v1_upperLim60.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.v1.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


cst.v1.pvals <- c()

##
## CST I at V1
##
cstI <- ifelse(mt.V1$cst=="I", 1, 0)
length(cstI) # 524

mI.dat <- list(y=cstI,
              delSt=mt.V1$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.v1.pvals[1] <- coefs.I[2]

##
## CST II at V1
##
cstII <- ifelse(mt.V1$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.V1$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.v1.pvals[2] <- coefs.II[2]


##
## CST III at V1
##
cstIII <- ifelse(mt.V1$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.V1$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.v1.pvals[3] <- coefs.III[2]

##
## CST IV-A at V1
##
cstIVA <- ifelse(mt.V1$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.V1$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.v1.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V1
##

cstIVB <- ifelse(mt.V1$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.V1$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.v1.pvals[5] <- coefs.IVB[2]

##
## CST V at V1
##

cstV <- ifelse(mt.V1$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.V1$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.v1.pvals[6] <- coefs.V[2]


cst.v1.perc <- rbind(cst.v1.perc, cst.v1.pvals)

rownames(cst.v1.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_V1.tex"
tblLabel <- "cst.delSt.v1:tbl"
lt <- latex(formatC(cst.v1.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.v1.perc, digits=3), file="../docs/cst_delSt_perc_v1.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_v1.tex"
tblLabel <- "cst.delSt.freq.v1:tbl"
lt <- latex(cst.v1.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.v1.freq, file="../docs/cst_delSt_freq_v1.csv", quote=F)



## ---------------------------------------------
##     African-American - all visits
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & !is.na(mt$raceChar) & mt$raceChar=="black"
mt.AA <- mt[idx,]

(cst.AA.freq <- table(mt.AA$delStatus, mt.AA$cst))
(cst.AA.perc <- 100 * cst.AA.freq / rowSums(cst.AA.freq))

write.csv(cst.AA.freq, file="../docs/cst_delSt_freq_AA.csv", quote=F)

cst.AA.pvals <- c()

##
## CST I - all visits
##
cstI <- ifelse(mt.AA$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.AA$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.AA.pvals[1] <- coefs.I[2]

##
## CST II - all visits
##
cstII <- ifelse(mt.AA$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.AA$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.AA.pvals[2] <- coefs.II[2]


##
## CST III - all visits
##
cstIII <- ifelse(mt.AA$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.AA$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.AA.pvals[3] <- coefs.III[2]

##
## CST IV-A - all visits
##
cstIVA <- ifelse(mt.AA$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.AA$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.AA.pvals[4] <- coefs.IVA[2]


##
## CST IV-B - all visits
##

cstIVB <- ifelse(mt.AA$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.AA$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.AA.pvals[5] <- coefs.IVB[2]

##
## CST V - all visits
##

cstV <- ifelse(mt.AA$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.AA$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.AA.pvals[6] <- coefs.V[2]


cst.AA.perc <- rbind(cst.AA.perc, cst.AA.pvals)

rownames(cst.AA.perc)[3] <- "p-val"

write.csv(formatC(cst.AA.perc, digits=3), file="../docs/cst_delSt_perc_AA.csv", quote=F)


file <- "../pics/AA_all_visits_upperLim60.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.AA.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


## ---------------------------------------------
##     African-American at visit V1
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V1" & !is.na(mt$raceChar) & mt$raceChar=="black"
mt.AA.v1 <- mt[idx,]

(cst.AA.v1.freq <- table(mt.AA.v1$delStatus, mt.AA.v1$cst))
##       I II III IV-A IV-B  V
## sPTB 15  3  11   27   16  5
## TERM 64 32  70   81   57  9
(cst.AA.v1.perc <- 100 * cst.AA.v1.freq / rowSums(cst.AA.v1.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 19.480519  3.896104 14.285714 35.064935 20.779221  6.493506
## TERM 20.447284 10.223642 22.364217 25.878594 18.210863  2.875399


file <- "../pics/AA_v1_upperLim60.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.AA.v1.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


cst.AA.v1.pvals <- c()

##
## CST I at V1
##
cstI <- ifelse(mt.AA.v1$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.AA.v1$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.AA.v1.pvals[1] <- coefs.I[2]

##
## CST II at V1
##
cstII <- ifelse(mt.AA.v1$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.AA.v1$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.AA.v1.pvals[2] <- coefs.II[2]


##
## CST III at V1
##
cstIII <- ifelse(mt.AA.v1$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.AA.v1$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.AA.v1.pvals[3] <- coefs.III[2]

##
## CST IV-A at V1
##
cstIVA <- ifelse(mt.AA.v1$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.AA.v1$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.AA.v1.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V1
##

cstIVB <- ifelse(mt.AA.v1$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.AA.v1$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.AA.v1.pvals[5] <- coefs.IVB[2]

##
## CST V at V1
##

cstV <- ifelse(mt.AA.v1$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.AA.v1$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.AA.v1.pvals[6] <- coefs.V[2]


cst.AA.v1.perc <- rbind(cst.AA.v1.perc, cst.AA.v1.pvals)

rownames(cst.AA.v1.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_AA_V1.tex"
tblLabel <- "cst.delSt.AA.v1:tbl"
lt <- latex(formatC(cst.AA.v1.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.AA.v1.perc, digits=3), file="../docs/cst_delSt_perc_AA_v1.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_AA_v1.tex"
tblLabel <- "cst.delSt.freq.AA.v1:tbl"
lt <- latex(cst.AA.v1.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.AA.v1.freq, file="../docs/cst_delSt_freq_AA_v1.csv", quote=F)



## ---------------------------------------------
##     non-African-American - all visits
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & !is.na(mt$raceChar) & mt$raceChar!="black"
mt.nonAA <- mt[idx,]

(cst.nonAA.freq <- table(mt.nonAA$delStatus, mt.nonAA$cst))
(cst.nonAA.perc <- 100 * cst.nonAA.freq / rowSums(cst.nonAA.freq))

write.csv(cst.nonAA.freq, file="../docs/cst_delSt_freq_nonAA.csv", quote=F)

cst.nonAA.pvals <- c()

##
## CST I - all visits
##
cstI <- ifelse(mt.nonAA$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.nonAA$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.nonAA.pvals[1] <- coefs.I[2]

##
## CST II - all visits
##
cstII <- ifelse(mt.nonAA$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.nonAA$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.nonAA.pvals[2] <- coefs.II[2]


##
## CST III - all visits
##
cstIII <- ifelse(mt.nonAA$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.nonAA$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.nonAA.pvals[3] <- coefs.III[2]

##
## CST IV-A - all visits
##
cstIVA <- ifelse(mt.nonAA$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.nonAA$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.nonAA.pvals[4] <- coefs.IVA[2]


##
## CST IV-B - all visits
##

cstIVB <- ifelse(mt.nonAA$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.nonAA$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.nonAA.pvals[5] <- coefs.IVB[2]

##
## CST V - all visits
##

cstV <- ifelse(mt.nonAA$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.nonAA$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.nonAA.pvals[6] <- coefs.V[2]


cst.nonAA.perc <- rbind(cst.nonAA.perc, cst.nonAA.pvals)

rownames(cst.nonAA.perc)[3] <- "p-val"

write.csv(formatC(cst.nonAA.perc, digits=3), file="../docs/cst_delSt_perc_nonAA.csv", quote=F)



file <- "../pics/nonAA_all_visits.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.nonAA.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


## ---------------------------------------------
##     non-African-American at visit V1
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V1" & !is.na(mt$raceChar) & mt$raceChar!="black"
mt.nonAA.v1 <- mt[idx,]

(cst.nonAA.v1.freq <- table(mt.nonAA.v1$delStatus, mt.nonAA.v1$cst))
##       I II III IV-A IV-B  V
## sPTB 11  3   4    6    0  1
## TERM 57 10  14    7    3 17
(cst.nonAA.v1.perc <- 100 * cst.nonAA.v1.freq / rowSums(cst.nonAA.v1.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 44.000000 12.000000 16.000000 24.000000  0.000000  4.000000
## TERM 52.777778  9.259259 12.962963  6.481481  2.777778 15.740741

file <- "../pics/nonAA_v1.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.nonAA.v1.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()

cst.nonAA.v1.pvals <- c()

##
## CST I at V1
##
cstI <- ifelse(mt.nonAA.v1$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.nonAA.v1$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.nonAA.v1.pvals[1] <- coefs.I[2]

##
## CST II at V1
##
cstII <- ifelse(mt.nonAA.v1$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.nonAA.v1$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.nonAA.v1.pvals[2] <- coefs.II[2]


##
## CST III at V1
##
cstIII <- ifelse(mt.nonAA.v1$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.nonAA.v1$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.nonAA.v1.pvals[3] <- coefs.III[2]

##
## CST IV-A at V1
##
cstIVA <- ifelse(mt.nonAA.v1$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.nonAA.v1$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.nonAA.v1.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V1
##

cstIVB <- ifelse(mt.nonAA.v1$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.nonAA.v1$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.nonAA.v1.pvals[5] <- coefs.IVB[2]

##
## CST V at V1
##

cstV <- ifelse(mt.nonAA.v1$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.nonAA.v1$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.nonAA.v1.pvals[6] <- coefs.V[2]


cst.nonAA.v1.perc <- rbind(cst.nonAA.v1.perc, cst.nonAA.v1.pvals)

rownames(cst.nonAA.v1.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_nonAA_V1.tex"
tblLabel <- "cst.delSt.nonAA.v1:tbl"
lt <- latex(formatC(cst.nonAA.v1.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.nonAA.v1.perc, digits=3), file="../docs/cst_delSt_perc_nonAA_v1.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_nonAA_v1.tex"
tblLabel <- "cst.delSt.freq.nonAA.v1:tbl"
lt <- latex(cst.nonAA.v1.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.nonAA.v1.freq, file="../docs/cst_delSt_freq_nonAA_v1.csv", quote=F)
