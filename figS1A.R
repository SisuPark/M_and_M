## ---------------------------------
##    CSTs in AA vs nonAA at V2
## ---------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V2"
mt.V2 <- mt[idx,]

(race.cst.v2.freq <- table(mt.V2$race2, mt.V2$cst))
(race.cst.v2.perc <- 100 * race.cst.v2.freq / rowSums(race.cst.v2.freq))

file <- "../pics/race_v2.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(race.cst.v2.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("AA","non-AA"))
box()
par(op)
dev.off()


table(mt.V2$race2)
##  AA nonAA
## 390   133

## using AA as a baseline class



race.cst.v2.pvals <- c()

##
## CST I at V2
##
cstI <- ifelse(mt.V2$cst=="I", 1, 0)
length(cstI) # 524
mI.dat <- list(y=cstI,
              race=mt.V2$race2)
str(mI.dat)

m.I <- glm( y ~ race, data=mI.dat, family=binomial )
summary(m.I)

race.cst.v2.pvals[1] <- coef(summary(m.I))[2,4]


##
## CST II at V2
##
cstII <- ifelse(mt.V2$cst=="II", 1, 0)
length(cstII) # 524
mII.dat <- list(y=cstII,
              race=mt.V2$race2)
str(mII.dat)

m.II <- glm( y ~ race, data=mII.dat, family=binomial )
summary(m.II)

race.cst.v2.pvals[2] <- coef(summary(m.II))[2,4]


##
## CST III at V2
##
cstIII <- ifelse(mt.V2$cst=="III", 1, 0)
length(cstIII) # 524
mIII.dat <- list(y=cstIII,
              race=mt.V2$race2)
str(mIII.dat)

m.III <- glm( y ~ race, data=mIII.dat, family=binomial )
summary(m.III)

race.cst.v2.pvals[3] <- coef(summary(m.III))[2,4]


##
## CST IV-A at V2
##
cstIVA <- ifelse(mt.V2$cst=="IV-A", 1, 0)
length(cstIVA) # 524
mIVA.dat <- list(y=cstIVA,
              race=mt.V2$race2)
str(mIVA.dat)

m.IVA <- glm( y ~ race, data=mIVA.dat, family=binomial )
summary(m.IVA)

race.cst.v2.pvals[4] <- coef(summary(m.IVA))[2,4]


##
## CST IV-B at V2
##
cstIVB <- ifelse(mt.V2$cst=="IV-B", 1, 0)
length(cstIVB) # 524
mIVB.dat <- list(y=cstIVB,
              race=mt.V2$race2)
str(mIVB.dat)

m.IVB <- glm( y ~ race, data=mIVB.dat, family=binomial )
summary(m.IVB)

race.cst.v2.pvals[5] <- coef(summary(m.IVB))[2,4]


##
## CST V at V2
##
cstV <- ifelse(mt.V2$cst=="V", 1, 0)
length(cstV) # 524
mV.dat <- list(y=cstV,
              race=mt.V2$race2)
str(mV.dat)

m.V <- glm( y ~ race, data=mV.dat, family=binomial )
summary(m.V)

race.cst.v2.pvals[6] <- coef(summary(m.V))[2,4]


race.cst.v2.perc <- rbind(race.cst.v2.perc, race.cst.v2.pvals)
rownames(race.cst.v2.perc)[3] <- "p-vals"


tblTexFile <- "../docs/cst_v2_race_freq.tex"
tblLabel <- "cst.v2.race.freq:tbl"
lt <- latex(race.cst.v2.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)

write.csv(race.cst.v2.freq, file="../docs/cst_v2_race_freq.csv", quote=F)


tblTexFile <- "../docs/cst_v2_race.tex"
tblLabel <- "cst.v2.race:tbl"
lt <- latex(formatC(race.cst.v2.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(race.cst.v2.perc, digits=3), file="../docs/cst_v2_race_perc.csv", quote=F)



## ---------------------------------
##    CSTs in AA vs nonAA at V3
## ---------------------------------

idx <- mt$delStatus.n != 2 & mt$visit=="V3"
mt.V3 <- mt[idx,]

(race.cst.v3.freq <- table(mt.V3$race2, mt.V3$cst))
(race.cst.v3.perc <- 100 * race.cst.v3.freq / rowSums(race.cst.v3.freq))

file <- "../pics/race_v3.pdf"
pdf(file, , width=6, height=6)
op <- par(mar=c(2.5, 3, 0.5, 1.5), mgp=c(2.5,0.6,0),tcl = -0.3)
barplot(race.cst.v3.perc, beside=T, las=1, ylim=c(0,60), legend.text=c("AA","non-AA"))
box()
par(op)
dev.off()


table(mt.V3$race2)
##  AA nonAA
## 390   133

## using AA as a baseline class



race.cst.v3.pvals <- c()

##
## CST I at V3
##
cstI <- ifelse(mt.V3$cst=="I", 1, 0)
length(cstI) # 524
mI.dat <- list(y=cstI,
              race=mt.V3$race2)
str(mI.dat)

m.I <- glm( y ~ race, data=mI.dat, family=binomial )
summary(m.I)

race.cst.v3.pvals[1] <- coef(summary(m.I))[2,4]


##
## CST II at V3
##
cstII <- ifelse(mt.V3$cst=="II", 1, 0)
length(cstII) # 524
mII.dat <- list(y=cstII,
              race=mt.V3$race2)
str(mII.dat)

m.II <- glm( y ~ race, data=mII.dat, family=binomial )
summary(m.II)

race.cst.v3.pvals[2] <- coef(summary(m.II))[2,4]


##
## CST III at V3
##
cstIII <- ifelse(mt.V3$cst=="III", 1, 0)
length(cstIII) # 524
mIII.dat <- list(y=cstIII,
              race=mt.V3$race2)
str(mIII.dat)

m.III <- glm( y ~ race, data=mIII.dat, family=binomial )
summary(m.III)

race.cst.v3.pvals[3] <- coef(summary(m.III))[2,4]


##
## CST IV-A at V3
##
cstIVA <- ifelse(mt.V3$cst=="IV-A", 1, 0)
length(cstIVA) # 524
mIVA.dat <- list(y=cstIVA,
              race=mt.V3$race2)
str(mIVA.dat)

m.IVA <- glm( y ~ race, data=mIVA.dat, family=binomial )
summary(m.IVA)

race.cst.v3.pvals[4] <- coef(summary(m.IVA))[2,4]


##
## CST IV-B at V3
##
cstIVB <- ifelse(mt.V3$cst=="IV-B", 1, 0)
length(cstIVB) # 524
mIVB.dat <- list(y=cstIVB,
              race=mt.V3$race2)
str(mIVB.dat)

m.IVB <- glm( y ~ race, data=mIVB.dat, family=binomial )
summary(m.IVB)

race.cst.v3.pvals[5] <- coef(summary(m.IVB))[2,4]


##
## CST V at V3
##
cstV <- ifelse(mt.V3$cst=="V", 1, 0)
length(cstV) # 524
mV.dat <- list(y=cstV,
              race=mt.V3$race2)
str(mV.dat)

m.V <- glm( y ~ race, data=mV.dat, family=binomial )
summary(m.V)

race.cst.v3.pvals[6] <- coef(summary(m.V))[2,4]


race.cst.v3.perc <- rbind(race.cst.v3.perc, race.cst.v3.pvals)
rownames(race.cst.v3.perc)[3] <- "p-vals"


tblTexFile <- "../docs/cst_v3_race_freq.tex"
tblLabel <- "cst.v3.race.freq:tbl"
lt <- latex(race.cst.v3.freq, file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)

write.csv(race.cst.v3.freq, file="../docs/cst_v3_race_freq.csv", quote=F)


tblTexFile <- "../docs/cst_v3_race.tex"
tblLabel <- "cst.v3.race:tbl"
lt <- latex(formatC(race.cst.v3.perc, digits=3), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


write.csv(formatC(race.cst.v3.perc, digits=3), file="../docs/cst_v3_race_perc.csv", quote=F)
