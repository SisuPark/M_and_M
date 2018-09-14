
## ------------------------------------------------
##    Comparing CST proportions in AA vs nonAA
## ------------------------------------------------

load("mm_no_mPTB_june25_2017.rda") # ct.no.mPTB, pt.no.mPTB, mt.no.mPTB, selPh.25

(race.cst.freq <- table(mt.no.mPTB$race2, mt.no.mPTB$cst))
(race.cst.perc <- 100 * race.cst.freq / rowSums(race.cst.freq))
##               I        II       III      IV-A      IV-B         V
## AA    23.215899  9.485095 20.957543 24.841915 18.066847  3.432701
## nonAA 52.020202 10.353535 12.878788  9.090909  2.272727 13.383838


race.cst.pvals <- c()

##
## CST I - all visits
##
cstI <- ifelse(mt.no.mPTB$cst=="I", 1, 0)
length(cstI) # 524
mI.dat <- list(y=cstI,
              race=mt.no.mPTB$race2)
str(mI.dat)

m.I <- glm( y ~ race, data=mI.dat, family=binomial )
summary(m.I)

race.cst.pvals[1] <- coef(summary(m.I))[2,4]


##
## CST II - all visits
##
cstII <- ifelse(mt.no.mPTB$cst=="II", 1, 0)
length(cstII) # 524
mII.dat <- list(y=cstII,
              race=mt.no.mPTB$race2)
str(mII.dat)

m.II <- glm( y ~ race, data=mII.dat, family=binomial )
summary(m.II)

race.cst.pvals[2] <- coef(summary(m.II))[2,4]


##
## CST III - all visits
##
cstIII <- ifelse(mt.no.mPTB$cst=="III", 1, 0)
length(cstIII) # 524
mIII.dat <- list(y=cstIII,
              race=mt.no.mPTB$race2)
str(mIII.dat)

m.III <- glm( y ~ race, data=mIII.dat, family=binomial )
summary(m.III)

race.cst.pvals[3] <- coef(summary(m.III))[2,4]


##
## CST IV-A - all visits
##
cstIVA <- ifelse(mt.no.mPTB$cst=="IV-A", 1, 0)
length(cstIVA) # 524
mIVA.dat <- list(y=cstIVA,
              race=mt.no.mPTB$race2)
str(mIVA.dat)

m.IVA <- glm( y ~ race, data=mIVA.dat, family=binomial )
summary(m.IVA)

race.cst.pvals[4] <- coef(summary(m.IVA))[2,4]


##
## CST IV-B - all visits
##
cstIVB <- ifelse(mt.no.mPTB$cst=="IV-B", 1, 0)
length(cstIVB) # 524
mIVB.dat <- list(y=cstIVB,
              race=mt.no.mPTB$race2)
str(mIVB.dat)

m.IVB <- glm( y ~ race, data=mIVB.dat, family=binomial )
summary(m.IVB)

race.cst.pvals[5] <- coef(summary(m.IVB))[2,4]


##
## CST V - all visits
##
cstV <- ifelse(mt.no.mPTB$cst=="V", 1, 0)
length(cstV) # 524
mV.dat <- list(y=cstV,
              race=mt.no.mPTB$race2)
str(mV.dat)

m.V <- glm( y ~ race, data=mV.dat, family=binomial )
summary(m.V)

race.cst.pvals[6] <- coef(summary(m.V))[2,4]


race.cst.perc <- rbind(race.cst.perc, race.cst.pvals)
rownames(race.cst.perc)[3] <- "p-vals"

write.csv(race.cst.freq, file="../docs/cst_race_freq.csv", quote=F)
write.csv(formatC(race.cst.perc, digits=3), file="../docs/cst_race_perc.csv", quote=F)


file <- "../pics/race_all_visits.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(race.cst.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("AA","non-AA"))
box()
par(op)
dev.off()



## ---------------------------------
##    CSTs in AA and nonAA at V1
## ---------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V1"
mt.V1 <- mt[idx,]

(race.cst.v1.freq <- table(mt.V1$race2, mt.V1$cst))
##         I  II III IV-A IV-B   V
## AA     79  35  81  108   73  14
## nonAA  68  13  18   13    3  18

(race.cst.v1.perc <- 100 * race.cst.v1.freq / rowSums(race.cst.v1.freq))
##               I        II       III      IV-A      IV-B         V
## AA    20.256410  8.974359 20.769231 27.692308 18.717949  3.589744
## nonAA 51.127820  9.774436 13.533835  9.774436  2.255639 13.533835

file <- "../pics/race_v1.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(race.cst.v1.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("AA","non-AA"))
box()
par(op)
dev.off()


table(mt.V1$race2)
##  AA nonAA
## 390   133

## using AA as a baseline class



race.cst.v1.pvals <- c()

##
## CST I at V1
##
cstI <- ifelse(mt.V1$cst=="I", 1, 0)
length(cstI) # 524
mI.dat <- list(y=cstI,
              race=mt.V1$race2)
str(mI.dat)

m.I <- glm( y ~ race, data=mI.dat, family=binomial )
summary(m.I)

race.cst.v1.pvals[1] <- coef(summary(m.I))[2,4]


##
## CST II at V1
##
cstII <- ifelse(mt.V1$cst=="II", 1, 0)
length(cstII) # 524
mII.dat <- list(y=cstII,
              race=mt.V1$race2)
str(mII.dat)

m.II <- glm( y ~ race, data=mII.dat, family=binomial )
summary(m.II)

race.cst.v1.pvals[2] <- coef(summary(m.II))[2,4]


##
## CST III at V1
##
cstIII <- ifelse(mt.V1$cst=="III", 1, 0)
length(cstIII) # 524
mIII.dat <- list(y=cstIII,
              race=mt.V1$race2)
str(mIII.dat)

m.III <- glm( y ~ race, data=mIII.dat, family=binomial )
summary(m.III)

race.cst.v1.pvals[3] <- coef(summary(m.III))[2,4]


##
## CST IV-A at V1
##
cstIVA <- ifelse(mt.V1$cst=="IV-A", 1, 0)
length(cstIVA) # 524
mIVA.dat <- list(y=cstIVA,
              race=mt.V1$race2)
str(mIVA.dat)

m.IVA <- glm( y ~ race, data=mIVA.dat, family=binomial )
summary(m.IVA)

race.cst.v1.pvals[4] <- coef(summary(m.IVA))[2,4]


##
## CST IV-B at V1
##
cstIVB <- ifelse(mt.V1$cst=="IV-B", 1, 0)
length(cstIVB) # 524
mIVB.dat <- list(y=cstIVB,
              race=mt.V1$race2)
str(mIVB.dat)

m.IVB <- glm( y ~ race, data=mIVB.dat, family=binomial )
summary(m.IVB)

race.cst.v1.pvals[5] <- coef(summary(m.IVB))[2,4]


##
## CST V at V1
##
cstV <- ifelse(mt.V1$cst=="V", 1, 0)
length(cstV) # 524
mV.dat <- list(y=cstV,
              race=mt.V1$race2)
str(mV.dat)

m.V <- glm( y ~ race, data=mV.dat, family=binomial )
summary(m.V)

race.cst.v1.pvals[6] <- coef(summary(m.V))[2,4]


race.cst.v1.perc <- rbind(race.cst.v1.perc, race.cst.v1.pvals)
rownames(race.cst.v1.perc)[3] <- "p-vals"


tblTexFile <- "../docs/cst_v1_race_freq.tex"
tblLabel <- "cst.v1.race.freq:tbl"
lt <- latex(race.cst.v1.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)

write.csv(race.cst.v1.freq, file="../docs/cst_v1_race_freq.csv", quote=F)


tblTexFile <- "../docs/cst_v1_race.tex"
tblLabel <- "cst.v1.race:tbl"
lt <- latex(formatC(race.cst.v1.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(race.cst.v1.perc, digits=3), file="../docs/cst_v1_race_perc.csv", quote=F)
