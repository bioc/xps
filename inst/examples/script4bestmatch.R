#------------------------------------------------------------------------------#
# Script: probesets with best match between HuExon, HuGene and HGU133P2
#         step-by-step functions to demonstrate how to compare results obtained
#         with xps to results obtained with apt (Affymetrix Power Tools)
#         using the Affymetrix human tissue-mixture exon array dataset
#
# Note: this script assumes that you have already produced normalized expression
#       levels using the accompanying script "script4exon.R"
#
# Copyright (c) 2008-2008 Dr. Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Functions 
#------------------------------------------------------------------------------#

## remove rows with duplicate probesets: keep only probeset with highest percentage
uniqueframe <- function(ma) {
   maxunique <- function(id, m) {
      m <- m[which(m[,1] == id),];
      m <- m[which(m[,2] == max(m[,2])),];
      return(m[1,]);
   }

   dup <- duplicated(ma[,1])
   uni <- unique(ma[dup,1])

   ds <- NULL
   for (i in uni) {ds <- rbind(ds, maxunique(i,ma))}

   tmp <- ma[dup==F,]
   ds <- rbind(ds, tmp)
   ds <- ds[order(rownames(ds)),]
   return(ds)
}#uniqueframe


## return subframe ma1 containing only rownames common to ma1 and ma2
intersectframes <- function(ma1, ma2) {
   id  <- intersect(rownames(ma1), rownames(ma2))
   ok  <- match(id, rownames(ma1))
   ok  <- ok[!is.na(ok)]
   ma1 <- ma1[ok,]
   ma1 <- ma1[order(rownames(ma1)),]
   return(ma1)
}#intersectframes


## return subframe ma1 containing only rownames common to ma1 and ma2
intersectrows <- function(ma1, ma2, col1=NULL, col2=NULL, ord=NULL) {
   if (is.null(col1)) {
      colA <- rownames(ma1)
   } else {
      colA <- ma1[,col1]
   }#if

   if (is.null(col2)) {
      col2 <- rownames(ma2)
   } else {
      col2 <- ma2[,col2]
   }#if

   id  <- intersect(colA, col2)
   ok  <- match(id, colA)
   ok  <- ok[!is.na(ok)]
   ma1 <- ma1[ok,]

   if (is.null(ord)) {
      ma1 <- ma1[order(rownames(ma1)),]
   } else if (ord > 0) {
      ma1 <- ma1[order(ma1[,ord]),]
   }#if

   return(ma1)
}#intersectrowss


#------------------------------------------------------------------------------#
# Affymetrix BestMatch tables 
#------------------------------------------------------------------------------#

## read BestMatch tables
hx2hg <- read.delim("HuExVsHuGene_BestMatch.txt", row.names=3, comment.char="")
hx2hg <- hx2hg[,5:6]
hx2hg <- uniqueframe(hx2hg)
colnames(hx2hg) <- c("HuGene", "PercentX2G")

up2hx <- read.delim("U133PlusVsHuEx_BestMatch.txt", row.names=3, comment.char="")
up2hx <- up2hx[,5:6]
up2hx <- uniqueframe(up2hx)
colnames(up2hx) <- c("HuExon", "PercentU2X")

up2hg <- read.delim("U133PlusVsHuGene_BestMatch.txt", row.names=3, comment.char="")
up2hg <- up2hg[,5:6]
up2hg <- uniqueframe(up2hg)
colnames(up2hg) <- c("HuGene", "PercentU2G")


## get unique exon and genome identical to u133p2
u2x  <- intersectframes(up2hx, up2hg)
u2g  <- intersectframes(up2hg, up2hx)
u2gx <- cbind(u2g[,1,drop=F],u2x[,1,drop=F])
write.table(u2x,"U133P2_Exon.txt",sep="\t",col.names=NA)
write.table(u2g,"U133P2_Genome.txt",sep="\t",col.names=NA)
write.table(u2gx,"U133P2_Genome_Exon.txt",sep="\t",col.names=NA)


#------------------------------------------------------------------------------#
# RMA normalized expression levels 
#------------------------------------------------------------------------------#

## import RMA data from xps (from text file created when running rma)
xps.urma <- read.delim("MixU133P2RMA.txt", row.names=2, comment.char="")
xps.grma <- read.delim("HuGeneMixRMAMetacore.txt", row.names=2, comment.char="")
xps.xrma <- read.delim("HuExonMixRMAMetacore.txt", row.names=2, comment.char="")
xps.grma.bq16 <- read.delim("HuGeneMixRMAbgqu16mp8.txt", row.names=2, comment.char="")
xps.xrma.bq16 <- read.delim("HuExonMixRMAbgqu16mp8.txt", row.names=2, comment.char="")
xps.urma <- xps.urma[,-1]
xps.grma <- xps.grma[,-1]
xps.xrma <- xps.xrma[,-1]
xps.grma.bq16 <- xps.grma.bq16[,-1]
xps.xrma.bq16 <- xps.xrma.bq16[,-1]

## read RMA data from apt
apt.urma <- read.delim("apt-u133p2-rma-summary.txt", row.names=1, comment.char="", skip=50)
apt.grma <- read.delim("apt-hugene-metacore-rma-summary.txt", row.names=1, comment.char="", skip=50)
apt.xrma <- read.delim("apt-huexon-metacore-rma-summary.txt", row.names=1, comment.char="", skip=50)


#------------------------------------------------------------------------------#
# Common unique metacore probesets 
#------------------------------------------------------------------------------#

## get metacore unique exon and genome identical to u133p2
tmp  <- intersectrows(u2gx, xps.grma, 1, NULL)
meta <- intersectrows(tmp,  xps.xrma, 2, NULL)


#------------------------------------------------------------------------------#
# RMA subsets with common unique metacore probesets 
#------------------------------------------------------------------------------#

## U133P2 metacore
xpsu <- intersectrows(xps.urma, meta, NULL, NULL)


## HuGene vs U133P2 xps metacore
xpsg <- intersectrows(xps.grma, meta, NULL, 1, -1)
xpsg <- xpsg[match(meta[,1],rownames(xpsg)),]

## HuGene vs U133P2 xps-bq16 metacore
bq16g <- intersectrows(xps.grma.bq16, meta, NULL, 1, -1)
bq16g <- bq16g[match(meta[,1],rownames(bq16g)),]

## HuGene vs U133P2 apt metacore
aptg <- intersectrows(apt.grma, meta, NULL, 1, -1)
aptg <- aptg[match(meta[,1],rownames(aptg)),]


## HuExon vs U133P2 xps metacore
xpsx <- intersectrows(xps.xrma, meta, NULL, 2, -1)
xpsx <- xpsx[match(meta[,2],rownames(xpsx)),]

## HuExon vs U133P2 xps-bq16 metacore
bq16x <- intersectrows(xps.xrma.bq16, meta, NULL, 2, -1)
bq16x <- bq16x[match(meta[,2],rownames(bq16x)),]

## HuExon vs U133P2 apt metacore
aptx <- intersectrows(apt.xrma, meta, NULL, 2, -1)
aptx <- aptx[match(meta[,2],rownames(aptx)),]


## export common meta probesets
tmp <- cbind(rownames(xpsu),rownames(xpsx),rownames(xpsg))
colnames(tmp) <- c("U133P2","HuExon","HuGene")
write.table(tmp,"U133P2_Genome_Exon_meta.txt",sep="\t",col.names=NA)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## HuGene vs U133P2 xps metacore
plot(xpsu[,1], xpsg[,1], log="xy")
plot(log2(xpsu[,1] * xpsg[,1])/2, log2(xpsg[,1]/xpsu[,1]), main="XPS: U133P2 vs HuGene", xlab="A = Log2(Gene*U133)", ylab="M = Log2(Gene/U133)")

## HuGene vs U133P2 xps-bq16 metacore
plot(xpsu[,1], bq16g[,1], log="xy")
plot(log2(xpsu[,1] * bq16g[,1])/2, log2(bq16g[,1]/xpsu[,1]), main="XPS-BQ16: U133P2 vs HuGene", xlab="A = Log2(Gene*U133)", ylab="M = Log2(Gene/U133)")

## HuGene vs U133P2 apt metacore
plot(xpsu[,1], aptg[,1], log="xy")
plot(log2(xpsu[,1] * aptg[,1])/2, log2(aptg[,1]/xpsu[,1]), main="APT: U133P2 vs HuGene", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",)


## HuExon vs U133P2 xps metacore
plot(xpsu[,1], xpsx[,1], log="xy")
plot(log2(xpsu[,1] * xpsx[,1])/2, log2(xpsx[,1]/xpsu[,1]), main="XPS: U133P2 vs HuExon", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)")

## HuExon vs U133P2 xps-bq16 metacore
plot(xpsu[,1], bq16x[,1], log="xy")
plot(log2(xpsu[,1] * bq16x[,1])/2, log2(bq16x[,1]/xpsu[,1]), main="XPS-BQ16: U133P2 vs HuExon", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)")

## HuExon vs U133P2 apt metacore
plot(xpsu[,1], aptx[,1], log="xy")
plot(log2(xpsu[,1] * aptx[,1])/2, log2(aptx[,1]/xpsu[,1]), main="APT: U133P2 vs HuExon", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",)

plot(aptx[,1], xpsx[,1], log="xy")
plot(log2(aptx[,1] * xpsx[,1])/2, log2(xpsx[,1]/aptx[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)")


## HuExon vs HuGene xps metacore
plot(xpsx[,1], xpsg[,1], log="xy")
plot(log2(xpsx[,1] * xpsg[,1])/2, log2(xpsg[,1]/xpsx[,1]), main="XPS: HuExon vs HuGene", xlab="A = Log2(Exon*Gene)", ylab="M = Log2(Gene/Exon)")
plot(log2(bg16x[,1] * bg16g[,1])/2, log2(bg16g[,1]/bg16x[,1]), main="XPS-BQ16: HuExon vs HuGene", xlab="A = Log2(Exon*Gene)", ylab="M = Log2(Gene/Exon)")
plot(log2(xpsg[,1] * bq16g[,1])/2, log2(bq16g[,1]/xpsg[,1]), main="HuGene: XPS-BQ16 vs XPS", xlab="A = Log2(XPS-BQ16*XPS)", ylab="M = Log2(XPS-BQ16/XPS)")


## use mean of triplicates!!! see paper M.Robinson & T.Speed
plot(xpsu[,1], bg16g[,1], log="xy")
plot(rowMeans(bg16g), log="xy")

## pairs xps
tmp <- cbind(rowMeans(xpsu), rowMeans(xpsg), rowMeans(xpsx))
colnames(tmp) <- c("XPS.U133P2","XPS.HuGene","XPS.HuExon")
png(file="XPS_U133P2_Gene_Exon_mean.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()
png(file="XPS_U133P2_Gene_Exon_mean_dot.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp),pch='.')
dev.off()

## pairs apt  (xpsu is identical to aptu)
tmp <- cbind(rowMeans(xpsu), rowMeans(aptg), rowMeans(aptx))
colnames(tmp) <- c("APT.U133P2","APT.HuGene","APT.HuExon")
png(file="APT_U133P2_Gene_Exon_mean.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()

## pairs xps vs apt
tmp <- cbind(rowMeans(xpsg), rowMeans(aptg), rowMeans(xpsx), rowMeans(aptx))
colnames(tmp) <- c("XPS.HuGene","APT.HuGene","XPS.HuExon","APT.HuExon")
png(file="XPS_APT_Gene_Exon_mean.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()


#------------------------------------------------------------------------------#
# Limma: compare breast vs prostate
#------------------------------------------------------------------------------#

library(limma)

## create design matrix
tissue <- c("Breast","Breast","Breast","Prostate","Prostate","Prostate") 
design <- model.matrix(~factor(tissue)) 
colnames(design) <- c("Breast","BreastvsProstate") 


## limma U133P2: xpsu
tmp <- as.matrix(log2(xpsu))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
xpsu.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
xpsu.lm <- xpsu.lm[order(xpsu.lm[,"ID"]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(xpsu.lm) <- c("xpsu.ID","xpsu.logFC","xpsu.P.Value","xpsu.adj.P.Val")


## limma HuGene vs U133P2 xps metacore: xpsg
tmp <- as.matrix(log2(xpsg))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
xpsg.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
xpsg.lm <- xpsg.lm[match(meta[,1],xpsg.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(xpsg.lm) <- c("xpsg.ID","xpsg.logFC","xpsg.P.Value","xpsg.adj.P.Val")

## limma HuGene vs U133P2 xps-bq16 metacore: bq16g
tmp <- as.matrix(log2(bq16g))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
bq16g.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
bq16g.lm <- bq16g.lm[match(meta[,1],bq16g.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(bq16g.lm) <- c("bq16g.ID","bq16g.logFC","bq16g.P.Value","bq16g.adj.P.Val")

## limma HuGene vs U133P2 apt metacore: aptg
tmp <- as.matrix(log2(aptg))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
aptg.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
aptg.lm <- aptg.lm[match(meta[,1],aptg.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(aptg.lm) <- c("aptg.ID","aptg.logFC","aptg.P.Value","aptg.adj.P.Val")


## limma HuExon vs U133P2 xps metacore: xpsx
tmp <- as.matrix(log2(xpsx))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
xpsx.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
xpsx.lm <- xpsx.lm[match(meta[,2],xpsx.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(xpsx.lm) <- c("xpsx.ID","xpsx.logFC","xpsx.P.Value","xpsx.adj.P.Val")

## limma HuExon vs U133P2 xps-bq16 metacore: bq16x
tmp <- as.matrix(log2(bq16x))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
bq16x.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
bq16x.lm <- bq16x.lm[match(meta[,2],bq16x.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(bq16x.lm) <- c("bq16x.ID","bq16x.logFC","bq16x.P.Value","bq16x.adj.P.Val")

## limma HuExon vs U133P2 apt metacore: aptx
tmp <- as.matrix(log2(aptx))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
aptx.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
aptx.lm <- aptx.lm[match(meta[,2],aptx.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(aptx.lm) <- c("aptx.ID","aptx.logFC","aptx.P.Value","aptx.adj.P.Val")


## export all
meta.lm <- cbind(xpsu.lm, xpsg.lm, bq16g.lm, aptg.lm, xpsx.lm, bq16x.lm, aptx.lm)
write.table(meta.lm,"Limma_U133P2_Genome_Exon_meta.txt",sep="\t",col.names=NA)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## plot HuGene P-Value for xps vs apt
png(file="XPS_APT_HuGene_PValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(meta.lm[,"xpsg.P.Value"], meta.lm[,"aptg.P.Value"],xlab="XPS: P-Value",ylab="APT: P-Value")
plot(meta.lm[,"xpsg.P.Value"], meta.lm[,"aptg.P.Value"],xlab="XPS: P-Value",ylab="APT: P-Value",xlim=c(0,0.001),ylim=c(0,0.001))
plot(meta.lm[,"bq16g.P.Value"], meta.lm[,"aptg.P.Value"],xlab="XPS-BQ16: P-Value",ylab="APT: P-Value")
plot(meta.lm[,"bq16g.P.Value"], meta.lm[,"aptg.P.Value"],xlab="XPS-BQ16: P-Value",ylab="APT: P-Value",xlim=c(0,0.001),ylim=c(0,0.001))
par(oldpar);
dev.off()

## plot HuGene adjusted P-Value for xps vs apt
png(file="XPS_APT_HuGene_AdjPValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(meta.lm[,"xpsg.adj.P.Val"], meta.lm[,"aptg.adj.P.Val"],xlab="XPS: Adjusted P-Value",ylab="APT: Adjusted P-Value")
plot(meta.lm[,"xpsg.adj.P.Val"], meta.lm[,"aptg.adj.P.Val"],xlab="XPS: Adjusted P-Value",ylab="APT: Adjusted P-Value",xlim=c(0,0.01),ylim=c(0,0.01))
plot(meta.lm[,"bq16g.adj.P.Val"], meta.lm[,"aptg.adj.P.Val"],xlab="XPS-BQ16: Adjusted P-Value",ylab="APT: Adjusted P-Value")
plot(meta.lm[,"bq16g.adj.P.Val"], meta.lm[,"aptg.adj.P.Val"],xlab="XPS-BQ16: Adjusted P-Value",ylab="APT: Adjusted P-Value",xlim=c(0,0.01),ylim=c(0,0.01))
par(oldpar);
dev.off()


## plot HuExon P-Value for xps vs apt
png(file="XPS_APT_HuExon_PValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(meta.lm[,"xpsx.P.Value"], meta.lm[,"aptx.P.Value"],xlab="XPS: P-Value",ylab="APT: P-Value")
plot(meta.lm[,"xpsx.P.Value"], meta.lm[,"aptx.P.Value"],xlab="XPS: P-Value",ylab="APT: P-Value",xlim=c(0,0.000001),ylim=c(0,0.000001))
plot(meta.lm[,"bq16x.P.Value"], meta.lm[,"aptx.P.Value"],xlab="XPS-BQ16: P-Value",ylab="APT: P-Value")
plot(meta.lm[,"bq16x.P.Value"], meta.lm[,"aptx.P.Value"],xlab="XPS-BQ16: P-Value",ylab="APT: P-Value",xlim=c(0,0.000001),ylim=c(0,0.000001))
par(oldpar);
dev.off()

## plot HuExon adjusted P-Value for xps vs apt
png(file="XPS_APT_HuExon_AdjPValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(meta.lm[,"xpsx.adj.P.Val"], meta.lm[,"aptx.adj.P.Val"],xlab="XPS: Adjusted P-Value",ylab="APT: Adjusted P-Value")
plot(meta.lm[,"xpsx.adj.P.Val"], meta.lm[,"aptx.adj.P.Val"],xlab="XPS: Adjusted P-Value",ylab="APT: Adjusted P-Value",xlim=c(0,0.0001),ylim=c(0,0.0001))
plot(meta.lm[,"bq16x.adj.P.Val"], meta.lm[,"aptx.adj.P.Val"],xlab="XPS-BQ16: Adjusted P-Value",ylab="APT: Adjusted P-Value")
plot(meta.lm[,"bq16x.adj.P.Val"], meta.lm[,"aptx.adj.P.Val"],xlab="XPS-BQ16: Adjusted P-Value",ylab="APT: Adjusted P-Value",xlim=c(0,0.0001),ylim=c(0,0.0001))
par(oldpar);
dev.off()

## plot logFC for xps
tmp <- cbind(meta.lm[,"xpsu.logFC"], meta.lm[,"xpsg.logFC"], meta.lm[,"xpsx.logFC"])
colnames(tmp) <- c("U133P2.logFC","HuGene.logFC","HuExon.logFC")
png(file="XPS_U133P2_HuGene_HuExon_logFC.png",width=440,height=440)
pairs(tmp, labels=colnames(tmp))
dev.off()

## plot logFC for apt
tmp <- cbind(meta.lm[,"xpsu.logFC"], meta.lm[,"aptg.logFC"], meta.lm[,"aptx.logFC"])
colnames(tmp) <- c("U133P2.logFC","HuGene.logFC","HuExon.logFC")
png(file="APT_U133P2_HuGene_HuExon_logFC.png",width=440,height=440)
pairs(tmp, labels=colnames(tmp))
dev.off()


#------------------------------------------------------------------------------#

