
load("mm_dada2_ct_pt.rda") # ct, pt, mt


## -----------------------------------
##     All subjects - Visit V2
## -----------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V2"
mt.V2 <- mt[idx,]

(cst.v2.freq <- table(mt.V2$delStatus, mt.V2$cst))
(cst.v2.perc <- 100 * cst.v2.freq / rowSums(cst.v2.freq))

file <- "../pics/v2_upperLim60.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.v2.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


cst.v2.pvals <- c()

##
## CST I at V2
##
cstI <- ifelse(mt.V2$cst=="I", 1, 0)
length(cstI) # 524

mI.dat <- list(y=cstI,
              delSt=mt.V2$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.v2.pvals[1] <- coefs.I[2]

##
## CST II at V2
##
cstII <- ifelse(mt.V2$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.V2$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.v2.pvals[2] <- coefs.II[2]


##
## CST III at V2
##
cstIII <- ifelse(mt.V2$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.V2$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.v2.pvals[3] <- coefs.III[2]

##
## CST IV-A at V2
##
cstIVA <- ifelse(mt.V2$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.V2$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.v2.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V2
##

cstIVB <- ifelse(mt.V2$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.V2$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.v2.pvals[5] <- coefs.IVB[2]

##
## CST V at V2
##

cstV <- ifelse(mt.V2$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.V2$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.v2.pvals[6] <- coefs.V[2]


cst.v2.perc <- rbind(cst.v2.perc, cst.v2.pvals)

rownames(cst.v2.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_V2.tex"
tblLabel <- "cst.delSt.v2:tbl"
lt <- latex(formatC(cst.v2.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.v2.perc, digits=3), file="../docs/cst_delSt_perc_v2.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_v2.tex"
tblLabel <- "cst.delSt.freq.v2:tbl"
lt <- latex(cst.v2.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.v2.freq, file="../docs/cst_delSt_freq_v2.csv", quote=F)




## ---------------------------------------------
##     African-American at visit V2
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V2" & !is.na(mt$raceChar) & mt$raceChar=="black"
mt.AA.v2 <- mt[idx,]

(cst.AA.v2.freq <- table(mt.AA.v2$delStatus, mt.AA.v2$cst))
##       I II III IV-A IV-B  V
## sPTB 15  3  11   27   16  5
## TERM 64 32  70   81   57  9
(cst.AA.v2.perc <- 100 * cst.AA.v2.freq / rowSums(cst.AA.v2.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 19.480519  3.896104 14.285714 35.064935 20.779221  6.493506
## TERM 20.447284 10.223642 22.364217 25.878594 18.210863  2.875399


file <- "../pics/AA_v2_upperLim60.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.AA.v2.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


cst.AA.v2.pvals <- c()

##
## CST I at V2
##
cstI <- ifelse(mt.AA.v2$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.AA.v2$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.AA.v2.pvals[1] <- coefs.I[2]

##
## CST II at V2
##
cstII <- ifelse(mt.AA.v2$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.AA.v2$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.AA.v2.pvals[2] <- coefs.II[2]


##
## CST III at V2
##
cstIII <- ifelse(mt.AA.v2$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.AA.v2$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.AA.v2.pvals[3] <- coefs.III[2]

##
## CST IV-A at V2
##
cstIVA <- ifelse(mt.AA.v2$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.AA.v2$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.AA.v2.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V2
##

cstIVB <- ifelse(mt.AA.v2$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.AA.v2$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.AA.v2.pvals[5] <- coefs.IVB[2]

##
## CST V at V2
##

cstV <- ifelse(mt.AA.v2$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.AA.v2$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.AA.v2.pvals[6] <- coefs.V[2]


cst.AA.v2.perc <- rbind(cst.AA.v2.perc, cst.AA.v2.pvals)

rownames(cst.AA.v2.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_AA_V2.tex"
tblLabel <- "cst.delSt.AA.v2:tbl"
lt <- latex(formatC(cst.AA.v2.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.AA.v2.perc, digits=3), file="../docs/cst_delSt_perc_AA_v2.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_AA_v2.tex"
tblLabel <- "cst.delSt.freq.AA.v2:tbl"
lt <- latex(cst.AA.v2.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.AA.v2.freq, file="../docs/cst_delSt_freq_AA_v2.csv", quote=F)



## ---------------------------------------------
##     non-African-American at visit V2
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V2" & !is.na(mt$raceChar) & mt$raceChar!="black"
mt.nonAA.v2 <- mt[idx,]

(cst.nonAA.v2.freq <- table(mt.nonAA.v2$delStatus, mt.nonAA.v2$cst))
##       I II III IV-A IV-B  V
## sPTB 11  3   4    6    0  1
## TERM 57 10  14    7    3 17
(cst.nonAA.v2.perc <- 100 * cst.nonAA.v2.freq / rowSums(cst.nonAA.v2.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 44.000000 12.000000 16.000000 24.000000  0.000000  4.000000
## TERM 52.777778  9.259259 12.962963  6.481481  2.777778 15.740741

file <- "../pics/nonAA_v2.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.nonAA.v2.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()

cst.nonAA.v2.pvals <- c()

##
## CST I at V2
##
cstI <- ifelse(mt.nonAA.v2$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.nonAA.v2$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.nonAA.v2.pvals[1] <- coefs.I[2]

##
## CST II at V2
##
cstII <- ifelse(mt.nonAA.v2$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.nonAA.v2$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.nonAA.v2.pvals[2] <- coefs.II[2]


##
## CST III at V2
##
cstIII <- ifelse(mt.nonAA.v2$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.nonAA.v2$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.nonAA.v2.pvals[3] <- coefs.III[2]

##
## CST IV-A at V2
##
cstIVA <- ifelse(mt.nonAA.v2$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.nonAA.v2$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.nonAA.v2.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V2
##

cstIVB <- ifelse(mt.nonAA.v2$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.nonAA.v2$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.nonAA.v2.pvals[5] <- coefs.IVB[2]

##
## CST V at V2
##

cstV <- ifelse(mt.nonAA.v2$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.nonAA.v2$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.nonAA.v2.pvals[6] <- coefs.V[2]


cst.nonAA.v2.perc <- rbind(cst.nonAA.v2.perc, cst.nonAA.v2.pvals)

rownames(cst.nonAA.v2.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_nonAA_V2.tex"
tblLabel <- "cst.delSt.nonAA.v2:tbl"
lt <- latex(formatC(cst.nonAA.v2.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.nonAA.v2.perc, digits=3), file="../docs/cst_delSt_perc_nonAA_v2.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_nonAA_v2.tex"
tblLabel <- "cst.delSt.freq.nonAA.v2:tbl"
lt <- latex(cst.nonAA.v2.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.nonAA.v2.freq, file="../docs/cst_delSt_freq_nonAA_v2.csv", quote=F)






## -----------------------------------
##     All subjects - Visit V3
## -----------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V3"
mt.V3 <- mt[idx,]

(cst.v3.freq <- table(mt.V3$delStatus, mt.V3$cst))
(cst.v3.perc <- 100 * cst.v3.freq / rowSums(cst.v3.freq))

file <- "../pics/v3_upperLim60.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.v3.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


cst.v3.pvals <- c()

##
## CST I at V3
##
cstI <- ifelse(mt.V3$cst=="I", 1, 0)
length(cstI) # 524

mI.dat <- list(y=cstI,
              delSt=mt.V3$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.v3.pvals[1] <- coefs.I[2]

##
## CST II at V3
##
cstII <- ifelse(mt.V3$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.V3$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.v3.pvals[2] <- coefs.II[2]


##
## CST III at V3
##
cstIII <- ifelse(mt.V3$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.V3$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.v3.pvals[3] <- coefs.III[2]

##
## CST IV-A at V3
##
cstIVA <- ifelse(mt.V3$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.V3$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.v3.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V3
##

cstIVB <- ifelse(mt.V3$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.V3$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.v3.pvals[5] <- coefs.IVB[2]

##
## CST V at V3
##

cstV <- ifelse(mt.V3$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.V3$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.v3.pvals[6] <- coefs.V[2]


cst.v3.perc <- rbind(cst.v3.perc, cst.v3.pvals)

rownames(cst.v3.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_V3.tex"
tblLabel <- "cst.delSt.v3:tbl"
lt <- latex(formatC(cst.v3.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.v3.perc, digits=3), file="../docs/cst_delSt_perc_v3.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_v3.tex"
tblLabel <- "cst.delSt.freq.v3:tbl"
lt <- latex(cst.v3.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.v3.freq, file="../docs/cst_delSt_freq_v3.csv", quote=F)




## ---------------------------------------------
##     African-American at visit V3
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V3" & !is.na(mt$raceChar) & mt$raceChar=="black"
mt.AA.v3 <- mt[idx,]

(cst.AA.v3.freq <- table(mt.AA.v3$delStatus, mt.AA.v3$cst))
##       I II III IV-A IV-B  V
## sPTB 15  3  11   27   16  5
## TERM 64 32  70   81   57  9
(cst.AA.v3.perc <- 100 * cst.AA.v3.freq / rowSums(cst.AA.v3.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 19.480519  3.896104 14.285714 35.064935 20.779221  6.493506
## TERM 20.447284 10.223642 22.364217 25.878594 18.210863  2.875399


file <- "../pics/AA_v3_upperLim60.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.AA.v3.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()


cst.AA.v3.pvals <- c()

##
## CST I at V3
##
cstI <- ifelse(mt.AA.v3$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.AA.v3$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.AA.v3.pvals[1] <- coefs.I[2]

##
## CST II at V3
##
cstII <- ifelse(mt.AA.v3$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.AA.v3$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.AA.v3.pvals[2] <- coefs.II[2]


##
## CST III at V3
##
cstIII <- ifelse(mt.AA.v3$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.AA.v3$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.AA.v3.pvals[3] <- coefs.III[2]

##
## CST IV-A at V3
##
cstIVA <- ifelse(mt.AA.v3$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.AA.v3$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.AA.v3.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V3
##

cstIVB <- ifelse(mt.AA.v3$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.AA.v3$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.AA.v3.pvals[5] <- coefs.IVB[2]

##
## CST V at V3
##

cstV <- ifelse(mt.AA.v3$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.AA.v3$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.AA.v3.pvals[6] <- coefs.V[2]


cst.AA.v3.perc <- rbind(cst.AA.v3.perc, cst.AA.v3.pvals)

rownames(cst.AA.v3.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_AA_V3.tex"
tblLabel <- "cst.delSt.AA.v3:tbl"
lt <- latex(formatC(cst.AA.v3.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.AA.v3.perc, digits=3), file="../docs/cst_delSt_perc_AA_v3.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_AA_v3.tex"
tblLabel <- "cst.delSt.freq.AA.v3:tbl"
lt <- latex(cst.AA.v3.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.AA.v3.freq, file="../docs/cst_delSt_freq_AA_v3.csv", quote=F)



## ---------------------------------------------
##     non-African-American at visit V3
## ---------------------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V3" & !is.na(mt$raceChar) & mt$raceChar!="black"
mt.nonAA.v3 <- mt[idx,]

(cst.nonAA.v3.freq <- table(mt.nonAA.v3$delStatus, mt.nonAA.v3$cst))
##       I II III IV-A IV-B  V
## sPTB 11  3   4    6    0  1
## TERM 57 10  14    7    3 17
(cst.nonAA.v3.perc <- 100 * cst.nonAA.v3.freq / rowSums(cst.nonAA.v3.freq))
##              I        II       III      IV-A      IV-B         V
## sPTB 44.000000 12.000000 16.000000 24.000000  0.000000  4.000000
## TERM 52.777778  9.259259 12.962963  6.481481  2.777778 15.740741

file <- "../pics/nonAA_v3.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(cst.nonAA.v3.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("sPTB","TERM"))
box()
par(op)
dev.off()

cst.nonAA.v3.pvals <- c()

##
## CST I at V3
##
cstI <- ifelse(mt.nonAA.v3$cst=="I", 1, 0)
length(cstI) # 390

mI.dat <- list(y=cstI,
              delSt=mt.nonAA.v3$delStatus)
str(mI.dat)

m.I <- glm( y ~ delSt, data=mI.dat, family=binomial )
summary(m.I)

(coefs.I <- coef(summary(m.I))[2,c(1,4)])
cst.nonAA.v3.pvals[1] <- coefs.I[2]

##
## CST II at V3
##
cstII <- ifelse(mt.nonAA.v3$cst=="II", 1, 0)
length(cstII) # 524

mII.dat <- list(y=cstII,
              delSt=mt.nonAA.v3$delStatus)
str(mII.dat)

m.II <- glm( y ~ delSt, data=mII.dat, family=binomial )
summary(m.II)

(coefs.II <- coef(summary(m.II))[2,c(1,4)])
cst.nonAA.v3.pvals[2] <- coefs.II[2]


##
## CST III at V3
##
cstIII <- ifelse(mt.nonAA.v3$cst=="III", 1, 0)
length(cstIII) # 524

mIII.dat <- list(y=cstIII,
              delSt=mt.nonAA.v3$delStatus)
str(mIII.dat)

m.III <- glm( y ~ delSt, data=mIII.dat, family=binomial )
summary(m.III)

(coefs.III <- coef(summary(m.III))[2,c(1,4)])
cst.nonAA.v3.pvals[3] <- coefs.III[2]

##
## CST IV-A at V3
##
cstIVA <- ifelse(mt.nonAA.v3$cst=="IV-A", 1, 0)
length(cstIVA) # 524

mIVA.dat <- list(y=cstIVA,
              delSt=mt.nonAA.v3$delStatus)
str(mIVA.dat)

m.IVA <- glm( y ~ delSt, data=mIVA.dat, family=binomial )
summary(m.IVA)

(coefs.IVA <- coef(summary(m.IVA))[2,c(1,4)])
cst.nonAA.v3.pvals[4] <- coefs.IVA[2]


##
## CST IV-B at V3
##

cstIVB <- ifelse(mt.nonAA.v3$cst=="IV-B", 1, 0)
length(cstIVB) # 524

mIVB.dat <- list(y=cstIVB,
              delSt=mt.nonAA.v3$delStatus)
str(mIVB.dat)

m.IVB <- glm( y ~ delSt, data=mIVB.dat, family=binomial )
summary(m.IVB)

(coefs.IVB <- coef(summary(m.IVB))[2,c(1,4)])
cst.nonAA.v3.pvals[5] <- coefs.IVB[2]

##
## CST V at V3
##

cstV <- ifelse(mt.nonAA.v3$cst=="V", 1, 0)
length(cstV) # 524

mV.dat <- list(y=cstV,
              delSt=mt.nonAA.v3$delStatus)
str(mV.dat)

m.V <- glm( y ~ delSt, data=mV.dat, family=binomial )
summary(m.V)

(coefs.V <- coef(summary(m.V))[2,c(1,4)])
cst.nonAA.v3.pvals[6] <- coefs.V[2]


cst.nonAA.v3.perc <- rbind(cst.nonAA.v3.perc, cst.nonAA.v3.pvals)

rownames(cst.nonAA.v3.perc)[3] <- "p-val"


tblTexFile <- "../docs/cst_delSt_nonAA_V3.tex"
tblLabel <- "cst.delSt.nonAA.v3:tbl"
lt <- latex(formatC(cst.nonAA.v3.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(cst.nonAA.v3.perc, digits=3), file="../docs/cst_delSt_perc_nonAA_v3.csv", quote=F)


tblTexFile <- "../docs/cst_delSt_freq_nonAA_v3.tex"
tblLabel <- "cst.delSt.freq.nonAA.v3:tbl"
lt <- latex(cst.nonAA.v3.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(cst.nonAA.v3.freq, file="../docs/cst_delSt_freq_nonAA_v3.csv", quote=F)
