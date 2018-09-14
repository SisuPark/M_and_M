##
## defining CSTs
##

load("mm_no_Thermodesulfatator_v3.rda") # ct, pt, selPh.10

dyn.load("vegdist/vegdist.so")
source("vegdist/vegdist.R")

JSdist <- vegdist(ct, method="jensen-shannon",useShrinkage=T)
hc <- hclust(as.dist(JSdist), method="ward.D2")

nClrs <- 6
cst <- cutree(hc, k=nClrs)
table(cst)

cstTbl <- c("IV-B","I","II","IV-A","III","V")
cst <- cstTbl[cst]
