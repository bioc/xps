#------------------------------------------------------------------------------#
# Script: step-by-step functions to demonstrate how to compare results obtained
#         with xps to results obtained with apt (Affymetrix Power Tools)
#         using the Affymetrix human tissue-mixture exon array dataset
#
# Note: this script assumes that you have already produced normalized expression
#       levels using the accompanying script "script4exon.R"
#
# Copyright (c) 2008-2008 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
#------------------------------------------------------------------------------#

### open R session containing processed data for HG-U133_Plus_2
#   see 3.step of accompanying script "script4exon.R" how to create the data
#   R session should contain objects: data.rma, data.mas5, call.mas5

### load library xps
library(xps)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. RMA: comparison xps vs affy vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## apt rma: sketch-quantile
apt-probeset-summarize -a rma-bg,quant-norm.sketch=-1.usepm=true.bioc=true,pm-only,med-polish.expon=true -d HG-U133_Plus_2.cdf *.CEL

# import data.frame
apt.rma.sk <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)
apt.rma.sk <- apt.rma.sk[order(rownames(apt.rma.sk)),]

## apt rma: use quantile and not sketch-quantile
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -d HG-U133_Plus_2.cdf *.CEL

# import data.frame
apt.rma <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)
apt.rma <- apt.rma[order(rownames(apt.rma)),]

## xps rma: see "script4exon.R" how to get "data.u133p2"
data.rma <- rma(data.u133p2,"MixU133P2RMA",filedir=datdir,tmpdir="", background="pmonly",normalize=TRUE)
xps.rma  <- validData(data.rma)
xps.rma  <- xps.rma[order(rownames(xps.rma)),]

## affy rma:
affy.rma <- justRMA()
affy.rma <- 2^exprs(affy.rma)


## compare expression levels obtained with quantile-sketch vs quantile
png(file="U133P2_APT_RMA_RMAsketch_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * apt.rma.sk[,1])/2, log2(apt.rma.sk[,1]/apt.rma[,1]), main="APT: RMA vs RMA-sketch", xlab="A = Log2(SketchRMA*RMA)", ylab="M = Log2(SketchRMA/RMA)",log="",ylim=c(-0.1,0.1))
dev.off()

## compare expression levels obtained with xps vs apt
png(file="U133P2_XPS_APT_RMA_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * xps.rma[,1])/2, log2(xps.rma[,1]/apt.rma[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.1,0.1))
dev.off()

## plot xps vs affy vs apt
tmp <- cbind(xps.rma[,1],affy.rma[,1],apt.rma[,1])
colnames(tmp) <- c("xps.rma","affy.rma","apt.rma")
png(file="U133P2_XPS_AFFY_APT_RMA.png",width=400,height=400)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()

# plot percent difference
plot(1:54675, 100*(apt.rma[,1] - xps.rma[,1])/apt.rma[,1])
plot(1:54675, 100*(affy.rma[,1] - xps.rma[,1])/affy.rma[,1])
plot(1:54675, 100*(affy.rma[,1] - apt.rma[,1])/affy.rma[,1])



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. MAS5: comparison xps vs affy vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## apt mas5 signal
apt-probeset-summarize -a mas5-bg,pm-mm,mas5-signal -d HG-U133_Plus_2.cdf *.CEL

# apt rma: import data.frame
apt.mas5 <- read.delim("mas5-bg.pm-mm.mas5-signal.summary.txt", row.names=1, comment.char="", skip=50)
# note: as mentioned in FAQ, apt does not allow for signal level normalization across the chip.
#       since xps mas5 scales to sc=500, we need to normalize the apt results first
apt.mas5 <- apply(apt.mas5, 2, function(x){x*(500/mean(x, trim=0.02))})
apt.mas5 <- apt.mas5[order(rownames(apt.mas5)),]

## xps mas5: get data.frame
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5All",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500, update=TRUE)
xps.mas5 <- validData(data.mas5)
xps.mas5 <- xps.mas5[order(rownames(xps.mas5)),]

## affy mas5:
affy <- ReadAffy()
affy.mas5 <- mas5(affy, normalize=T, sc=500)
affy.mas5 <- exprs(affy.mas5)

## ExpressionConsole mas5:
ec.mas5 <- read.delim("EC_MAS500_signal.txt", row.names=1, comment.char="")
ec.mas5 <- ec.mas5[order(rownames(ec.mas5)),]

## compare expression levels obtained with xps, apt, ExpressionConsole, affy
png(file="U133P2_XPS_APT_EC_AFFY_MAS5_MvA_3x2.png",width=400,height=600)
oldpar <- par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(log2(xps.mas5[,1] * apt.mas5[,1])/2, log2(xps.mas5[,1]/apt.mas5[,1]), main="MAS5: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-1,1))
plot(log2(affy.mas5[,1] * apt.mas5[,1])/2, log2(affy.mas5[,1]/apt.mas5[,1]), main="MAS5: AFFY vs APT", xlab="A = Log2(AFFY*APT)", ylab="M = Log2(AFFY/APT)",log="",ylim=c(-1,1))
plot(log2(xps.mas5[,1] * ec.mas5[,1])/2, log2(xps.mas5[,1]/ec.mas5[,1]), main="MAS5: XPS vs EC", xlab="A = Log2(XPS*EC)", ylab="M = Log2(XPS/EC)",log="",ylim=c(-1,1))
plot(log2(affy.mas5[,1] * ec.mas5[,1])/2, log2(affy.mas5[,1]/ec.mas5[,1]), main="MAS5: AFFY vs EC", xlab="A = Log2(AFFY*EC)", ylab="M = Log2(AFFY/EC)",log="",ylim=c(-1,1))
plot(log2(apt.mas5[,1] * ec.mas5[,1])/2, log2(apt.mas5[,1]/ec.mas5[,1]), main="MAS5: APT vs EC", xlab="A = Log2(APT*EC)", ylab="M = Log2(APT/EC)",log="",ylim=c(-1,1))
plot(log2(xps.mas5[,1] * affy.mas5[,1])/2, log2(xps.mas5[,1]/affy.mas5[,1]), main="MAS5: XPS vs AFFY", xlab="A = Log2(XPS*AFFY)", ylab="M = Log2(XPS/AFFY)",log="",ylim=c(-1,1))
par(oldpar);
dev.off()

# compare expression levels obtained with xps, apt, ExpressionConsole, affy
png(file="U133P2_Diff_XPS_APT_EC_AFFY_MAS5_2x3.png",width=600,height=440)
oldpar <- par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(xps.mas5[,1], 100*(apt.mas5[,1] - xps.mas5[,1])/apt.mas5[,1], xlab="Intensity XPS", ylab="(APT - XPS) / APT [%]", log="x", ylim=c(-100,100))
plot(affy.mas5[,1], 100*(apt.mas5[,1] - affy.mas5[,1])/apt.mas5[,1], xlab="Intensity AFFY", ylab="(APT - AFFY) / APT [%]", log="x", ylim=c(-100,100))
plot(ec.mas5[,1], 100*(apt.mas5[,1] - ec.mas5[,1])/apt.mas5[,1], xlab="Intensity EC", ylab="(APT - EC) / APT [%]", log="x", ylim=c(-100,100))
plot(xps.mas5[,1], 100*(ec.mas5[,1] - xps.mas5[,1])/ec.mas5[,1], xlab="Intensity XPS", ylab="(EC - XPS) / EC [%]", log="x", ylim=c(-100,100))
plot(affy.mas5[,1], 100*(ec.mas5[,1] - affy.mas5[,1])/ec.mas5[,1], xlab="Intensity AFFY", ylab="(AFFY - EC) / EC [%]", log="x", ylim=c(-100,100))
plot(xps.mas5[,1], 100*(affy.mas5[,1] - xps.mas5[,1])/affy.mas5[,1], xlab="Intensity XPS", ylab="(AFFY - XPS) / AFFY [%]", log="x", ylim=c(-100,100))
par(oldpar);
dev.off()

## plot xps vs affy vs EC vs apt
tmp <- cbind(xps.mas5[,1],affy.mas5[,1],ec.mas5[,1],apt.mas5[,1])
colnames(tmp) <- c("xps.mas5","affy.mas5","ec.mas5","apt.mas5")
png(file="U133P2_XPS_AFFY_EC_APT_MAS5.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. MAS5 Call:  comparison xps vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### apt mas5 detection call
apt-probeset-summarize -a pm-mm,mas5-detect.calls=1.pairs=1 -d HG-U133_Plus_2.cdf -x 10 *.CEL

# apt mas5.call: import data.frame
apt.pval <- read.delim("pm-mm.mas5-detect.summary.txt", row.names=1, comment.char="", skip=50)
apt.pval <- apt.pval[order(rownames(apt.pval)),]

# xps mas5.call: get data.frame
call.mas5 <- mas5.call(data.u133p2,"MixU133P2Call")
xps.pval  <- pvalData(call.mas5)
rownames(xps.pval) <- xps.pval[,2]
xps.pval <- xps.pval[,-c(1:2)]
xps.pval <- xps.pval[order(rownames(xps.pval)),]

## affy mas5calls
affy <- ReadAffy()
affy.dc5  <- mas5calls(affy)
affy.pval <- assayData(affy.dc5)[["se.exprs"]]
#write.table(affy.pval,"BrPr_pval.txt",sep="\t",col.names=NA)
#affy.pval <- read.delim("BrPr_pval.txt", row.names=1, comment.char="")

## ExpressionConsole mas5 call:
ec.pval <- read.delim("EC_MAS500_pval.txt", row.names=1, comment.char="")
ec.pval <- ec.pval[order(rownames(ec.pval)),]

# compare detection calls obtained with xps vs affy
png(file="U133P2_Diff_XPS_AFFY_MAS5_pval.png",width=400,height=400)
plot(xps.pval[,1], 100*(affy.pval[,1] - xps.pval[,1])/affy.pval[,1], main="MAS5 P-Value: XPS vs AFFY", xlab="XPS P-Value", ylab="Difference [%]", ylim=c(-10,10))
dev.off()

# plot xps vs affy vs apt
tmp <- cbind(xps.pval[,1],affy.pval[,1],ec.pval[,1],apt.pval[,1])
colnames(tmp) <- c("xps.pval","affy.pval","ec.pval","apt.pval")
png(file="U133P2_XPS_AFFY_EC_APT_MAS5_pval.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp),xlim=c(-15,0),ylim=c(-15,0))
dev.off()


#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1 
#------------------------------------------------------------------------------#

### open R session containing processed data for HuGene-1_0-st-v1
#   see 3.step of accompanying script "script4exon.R" how to create the data
#   R session should contain objects: data.g.rma, call.g.dabg

### load library xps
library(xps)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. RMA: comparison xps vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## xps rma: see "script4exon.R" how to get "data.genome"
data.g.rma <- rma(data.genome,"HuGeneMixRMAMetacore",filedir=datdir,tmpdir="",
                  background="antigenomic",normalize=T,exonlevel="metacore+affx")
xps.rma  <- validData(data.g.rma)

## probesetList.txt
#  need to create probeset list for APT so that APT can use same probesets as XPS
tmp <- as.data.frame(rownames(xps.rma))
colnames(tmp) <- "probeset_id"
write.table(tmp, "probesetList.txt", quote=FALSE, row.names=FALSE)

## apt rma: use quantile and not sketch-quantile; use probesetList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r3.pgf -c HuGene-1_0-st-v1.r3.clf -s probesetList.txt *.CEL

# import data.frame
apt.rma <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)

# compare expression levels obtained with xps vs apt
png(file="HuGene_XPS_APT_RMA_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * xps.rma[,1])/2, log2(xps.rma[,1]/apt.rma[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,0.2))
dev.off()


## apt medpol only
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r3.pgf -c HuGene-1_0-st-v1.r3.clf -s probesetList.txt *.CEL

# import data.frame
apt.mp <- read.delim("pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)

## xps medpol only
expr.mp <- express(data.genome,"HuGeneMixMedPolMetacore",summarize.method="medianpolish",
           summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",
           summarize.params=c(10, 0.01, 1.0),exonlevel="metacore+affx")
xps.mp <- validData(expr.mp)

# compare medpol obtained with xps vs apt
png(file="HuGene_XPS_APT_MedPol_MvA.png",width=400,height=400)
plot(log2(apt.mp[,1] * xps.mp[,1])/2, log2(xps.mp[,1]/apt.mp[,1]), main="MedPol: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

## xps rma: bg="all", quantile="all", medpol="metacore+affx"
data.rma.bq16 <- rma(data.genome,"HuGeneMixRMAbgqu16mp8",filedir=datdir,tmpdir="",
                    background="antigenomic",normalize=T,exonlevel=c(16316,16316,8252))
xps.rma.bq16 <- validData(data.rma.bq16)

# compare expression levels obtained with xps vs apt
png(file="HuGene_XPS_APT_RMA_bg16qu16_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * xps.rma.bq16[,1])/2, log2(xps.rma.bq16[,1]/apt.rma[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,0.2))
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. DABG: comparison xps vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### apt detection above background call; probesetList.txt
apt-probeset-summarize -a dabg -p HuGene-1_0-st-v1.r3.pgf -c HuGene-1_0-st-v1.r3.clf -b HuGene-1_0-st-v1.r3.bgp -s probesetList.txt  -x 8 *.CEL

# apt dabg: import data.frame
apt.dabg <- read.delim("dabg.summary.txt", row.names=1, comment.char="", skip=50)

# xps dabg: see "script4exon.R" how to get "data.genome"
call.g.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore",filedir=datdir,
                         exonlevel="metacore+affx")
xps.pval <- pvalData(call.g.dabg)
rownames(xps.pval) <- xps.pval[,2]
xps.pval <- xps.pval[,-c(1:2)]

# compare detection calls obtained with xps vs affy
png(file="HuGene_Diff_XPS_APT_DABG_pval.png",width=400,height=400)
plot(xps.pval[,1], (xps.pval[,1] - apt.dabg[,1]), main="DABG P-Value: XPS vs APT", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.00001,0.00001))
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. MAS5 Call: comparison to DABG Call
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# xps mas5.call: see "script4exon.R" how to get "data.genome"
call.g.mas5 <- mas5.call(data.genome,"HuGeneMixCallMetacore",filedir=datdir,tmpdir="",
                         exonlevel="metacore+affx")
mas.pval <- pvalData(call.g.mas5)
rownames(mas.pval) <- mas.pval[,2]
mas.pval <- mas.pval[,-c(1:2)]

# compare DABG to MAS5 detection calls
png(file="HuGene_Diff_XPS_DABG_MAS5_pval.png",width=400,height=400)
plot(xps.pval[,1], (xps.pval[,1] - mas.pval[,1]), main="P-Value: DABG Call vs MAS5 Call", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.001,0.001))
dev.off()


#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
#------------------------------------------------------------------------------#

### open R session containing processed data for HuEx-1_0-st-v2
#   see 3.step of accompanying script "script4exon.R" how to create the data
#   R session should contain objects: data.x.rma.ps, call.x.dabg.ps

### load library xps
library(xps)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. RMA: comparison xps vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## xps rma: see "script4exon.R" how to get "data.exon"
data.x.rma.ps <- rma(data.exon,"MixRMAMetacorePS",filedir=datdir,tmpdir="",background="antigenomic",
                     normalize=T,option="probeset",exonlevel="metacore")
xps.rma.ps <- validData(data.x.rma.ps)

## metacorePSList.txt
#  need to create probeset list for APT so that APT can use same probesets as XPS
tmp <- as.data.frame(rownames(xps.rma.ps))
colnames(tmp) <- "probeset_id"
write.table(tmp, "metacorePSList.txt", quote=FALSE, row.names=FALSE)

## apt rma: use quantile and not sketch-quantile; use metacorePSList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -b HuEx-1_0-st-v2.r2.antigenomic.bgp -s metacorePSList.txt *.CEL

# import data.frame
apt.rma.ps <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)

# compare expression levels obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps[,1])/2, log2(xps.rma.ps[,1]/apt.rma.ps[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()


## apt medpol only
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -s metacorePSList.txt *.CEL

# import data.frame
apt.mp.ps <- read.delim("pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)

## xps medpol only
expr.mp.ps <- express(data.exon,"HuExonMixMedPolMetacore",summarize.method="medianpolish",
              summarize.select="pmonly",summarize.option="probeset",summarize.logbase="log2",
              summarize.params=c(10, 0.01, 1.0),exonlevel="metacore")
xps.mp.ps <- validData(expr.mp.ps)

# compare medpol obtained with xps vs apt
png(file="HuExon_XPS_APT_MedPol_MvA_ps.png",width=400,height=400)
plot(log2(apt.mp.ps[,1] * xps.mp.ps[,1])/2, log2(xps.mp.ps[,1]/apt.mp.ps[,1]), main="MedPol: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()


## xps rma: bg="all", quantile="all", medpol="metacore"
data.rma.ps.bq16 <- rma(data.exon,"HuExonMixRMAbgqu16mp8",filedir=datdir,tmpdir="",
                    background="antigenomic",normalize=T,option="probeset",exonlevel=c(16316,16316,8192))
xps.rma.ps.bq16 <- validData(data.rma.ps.bq16)

# compare expression levels obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_bg16qu16_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq16[,1])/2, log2(xps.rma.ps.bq16[,1]/apt.rma.ps[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()


## xps rma: bg="all", quantile="all", medpol="all"
data.rma.ps.all <- rma(data.exon,"MixRMAAllPS",filedir=datdir,tmpdir="",background="antigenomic",
                       normalize=T,option="probeset",exonlevel="all")
# create probeset list for apt
xps.rma.ps.all <- validData(data.rma.ps.all)
tmp <- as.data.frame(rownames(xps.rma.ps.all))
colnames(tmp) <- "probeset_id"
write.table(tmp, "allPSList.txt", quote=FALSE, row.names=FALSE)

## apt rma: use allPSList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -b HuEx-1_0-st-v2.r2.antigenomic.bgp -s allPSList.txt *.CEL

# import data.frame
apt.rma.ps.all <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)

# compare expression levels obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_MvA_ps_all.png",width=400,height=400)
plot(log2(apt.rma.ps.all[,1] * xps.rma.ps.all[,1])/2, log2(xps.rma.ps.all[,1]/apt.rma.ps.all[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()


## xps rma: see "script4exon.R" how to get "data.exon"
data.x.rma <- rma(data.exon,"MixRMAMetacore",filedir=datdir,tmpdir="",background="antigenomic",
                  normalize=T,option="transcript",exonlevel="metacore")
xps.rma  <- validData(data.x.rma)

### need to create probeset list for APT so that APT can use same probesets as XPS
# metacoreList.mps
xps.rma <- validData(data.x.rma)
writeLines(rownames(xps.rma), "metacore.txt")
metaProbesets(scheme.exon,"metacore.txt","metacoreList.mps","metacore")

## apt rma: use quantile and not sketch-quantile; use metacoreList.mps
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -m metacoreList.mps *.CEL

# import data.frame
apt.rma <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)


# compare expression levels obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * xps.rma[,1])/2, log2(xps.rma[,1]/apt.rma[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

## xps rma: bg="all", quantile="all", medpol="metacore+affx"
data.rma.bq16 <- rma(data.exon,"HuExonMixRMAbgqu16mp8",filedir=datdir,tmpdir="",
                     background="antigenomic",normalize=T,exonlevel=c(16316,16316,8192))
xps.rma.bq16 <- validData(data.rma.bq16)

png(file="HuExon_XPS_APT_RMA_bg16qu16_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * xps.rma.bq16[,1])/2, log2(xps.rma.bq16[,1]/apt.rma[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()


## apt medpol only
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -m metacoreList.mps *.CEL

# import data.frame
apt.mp <- read.delim("pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)

## xps medpol only
expr.mp <- express(data.exon,"HuExonMedPolMetacore",summarize.method="medianpolish",
           summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",
           summarize.params=c(10, 0.01, 1.0),exonlevel="metacore")
xps.mp <- validData(expr.mp)

# compare medpol obtained with xps vs apt
png(file="HuExon_XPS_APT_MedPol_MvA.png",width=400,height=400)
plot(log2(apt.mp[,1] * xps.mp[,1])/2, log2(xps.mp[,1]/apt.mp[,1]), main="MedPol: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()


## apt rma: sketch-quantile
apt-probeset-summarize -a rma-bg,quant-norm.sketch=-1.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -m metacoreList.mps *.CEL

# import data.frame
apt.rma.sk <- read.delim("rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=50)
apt.rma.sk <- apt.rma.sk[order(rownames(apt.rma.sk)),]

# compare expression levels obtained with xps vs apt
png(file="HuExon_APT_RMA_RMAsketch_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * apt.rma.sk[,1])/2, log2(apt.rma.sk[,1]/apt.rma[,1]), main="APT: RMA vs RMA-sketch", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.1,0.1))
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. DABG: comparison xps vs apt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### apt detection above background call; metacorePSList.txt
apt-probeset-summarize -a dabg -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -b HuEx-1_0-st-v2.r2.antigenomic.bgp -s metacorePSList.txt  -x 12 *.CEL

# apt dabg: import data.frame
apt.dabg.ps <- read.delim("dabg.summary.txt", row.names=1, comment.char="", skip=50)

# xps dabg: see "script4exon.R" how to get "data.exon"
call.x.dabg.ps <- dabg.call(data.exon,"MixDABGMetacorePS",filedir=datdir,option="probeset",exonlevel="metacore")
xps.pval.ps <- pvalData(call.x.dabg.ps)
rownames(xps.pval.ps) <- xps.pval.ps[,2]
xps.pval.ps <- xps.pval.ps[,-c(1:2)]

# compare detection calls obtained with xps vs affy
png(file="HuExon_Diff_XPS_APT_DABG_pval_ps.png",width=400,height=400)
plot(xps.pval.ps[,1], (xps.pval.ps[,1] - apt.dabg.ps[,1]), main="DABG P-Value: XPS vs APT", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.00001,0.00001))
dev.off()


# xps dabg at transcript level: see "script4exon.R" how to get "data.exon"
call.x.dabg <- dabg.call(data.exon,"MixDABGMetacore",filedir=datdir,option="transcript",exonlevel="metacore")
xps.pval <- pvalData(call.x.dabg)
rownames(xps.pval.ps) <- xps.pval[,2]
xps.pval <- xps.pval[,-c(1:2)]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Appendix: APT PLIER
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## apt plier: use quantile and not sketch-quantile
apt-probeset-summarize -a plier-mm -d HG-U133_Plus_2.cdf *.CEL

# import data.frame
apt.plier <- read.delim("plier-mm.summary.txt", row.names=1, comment.char="", skip=50)
apt.plier <- apt.plier[order(rownames(apt.plier)),]

## Scatter-plot for Breast_A vs Breast_B 
png(file="U133P2_APT_PLIER_Scatter.png",width=400,height=400)
plot(apt.plier[,1], apt.plier[,2], main="APT: PLIER", xlab="Breast_A", ylab="Breast_B",log="xy")
dev.off()

## MvA-plot for Breast_A vs Breast_B 
png(file="U133P2_APT_PLIER_MvA.png",width=400,height=400)
plot(log2(apt.plier[,1] * apt.plier[,2])/2, log2(apt.plier[,1]/apt.plier[,2]), main="APT: PLIER", xlab="A = Log2(BrA*BrB)", ylab="M = Log2(BrA/BrB)",log="")
dev.off()

## MvA-plot for Breast_A vs Breast_B using RMA
png(file="U133P2_APT_RMA_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * apt.rma[,2])/2, log2(apt.rma[,1]/apt.rma[,2]), main="APT: RMA", xlab="A = Log2(BrA*BrB)", ylab="M = Log2(BrA/BrB)",log="",ylim=c(-10,10))
dev.off()





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create metacoreList.mps using R (too slow: about 6-8hrs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### load library xps
library(xps)

### metacoreList.mps
#   need to create probeset list for APT so that APT can use same probesets as XPS
xps.rma <- validData(data.x.rma)
tc.id <- rownames(xps.rma)
writeLines(tc.id, "metacore.txt")

# get probeset_id, transcript_cluster_id, probe_count from probeset annotation tree
anp <- export.scheme(scheme.exon, treetype="anp", varlist="fProbesetID:fTranscriptID:fNProbes:fCrossHybType:fLevelID:fBounded", outfile="HuExon_anp.txt",as.dataframe=TRUE,verbose=FALSE)


# 1. trial: create metacoreList.mps using for-loop
mps <- vector(mode="character",length=(length(tc.id) + 1))
mps[1] <- paste("probeset_id", "transcript_cluster_id", "probeset_list", "probe_count", sep="\t")
Sys.time()
for (i in 1:100) {
   tmp <- anp[anp[,"TranscriptClusterID"] == tc.id[i],]
   tmp <- tmp[tmp[,"Level"] == "core",]
   tmp <- tmp[tmp[,"CrossHybType"] == 1,]
   tmp <- tmp[tmp[,"Bounded_NoBoundedEvidence"] == 0,]
   txt <- paste(tmp[,"ProbesetID"], collapse=" ")
   txt <- paste(tmp[1,"TranscriptClusterID"], tmp[1,"TranscriptClusterID"], txt, sum(tmp[,"ProbeCount"]), sep="\t")
   mps[i+1] <- txt
}
Sys.time()
writeLines(mps, "metacoreList.mps")


# 2. trial: create metacoreList.mps using sapply()
"mclist" <- function(x) {
   tmp <- subset(anp, anp[,"TranscriptClusterID"] == x)
   tmp <- subset(tmp, (tmp[,"Level"] == "core" & tmp[,"CrossHybType"] == 1 & tmp[,"Bounded_NoBoundedEvidence"] == 0))
   txt <- paste(tmp[,"ProbesetID"], collapse="\t")
   txt <- paste(tmp[1,"TranscriptClusterID"], tmp[1,"TranscriptClusterID"], txt, sum(tmp[,"ProbeCount"]), sep="\t")
   return(txt)
}

Sys.time()
mps <- sapply(tc.id, mclist)
Sys.time()
mps <- append(paste("probeset_id", "transcript_cluster_id", "probeset_list", "probe_count", sep="\t"), mps)
writeLines(mps, "metacoreList.mps")


#------------------------------------------------------------------------------#



