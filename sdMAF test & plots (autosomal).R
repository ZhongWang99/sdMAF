####################################################################################
# This file contains the process of running sdMAF testing between male and female,
#     which will generate a .csv file.
# Detailed plotting procedure also included.
# Build GRCh37, autosomal.
####################################################################################
####################################################################################
# Using simple linear regression, we regress Y(Sex) on 
# G(Genotype), when the test statistics is significant, 
# we say the difference shows
# Below are framework codes applicable for all autosomal data; One data are ready, Type 'Ctrl+F' and replace all 'chra' with 'chrb'
# replace all 'CHR=a' with 'CHR=b'
# replace all 'Chromosome a' with 'Chromosome b'
####################################################################################

library(qqman)
library(genio)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(cowplot)

chr7phase3.GSNP.bim = read_bim('data/chr7/chr7_SNP_only.recode.bim')
chr7phase3.GSNP.bim$chr = as.numeric(chr7phase3.GSNP.bim$chr)

chr7.AFtest = read.csv('chr7.AFtest.csv')

chr7.SNP.FreqTable = read.csv('chr7.SNP.FreqTable.csv')

# Initialize dataset & Simple linear regression
index1 = grep("F_A1A1", colnames(chr7.SNP.FreqTable))
afdiff.dt = chr7.SNP.FreqTable[,index1:(index1+5)]   

pr1x = apply(afdiff.dt,MARGIN = 1,FUN = wald.1df.hwd.auto)
pr2x = apply(afdiff.dt,MARGIN = 1,FUN = wald.1df.hwe.auto)


# Record Data
chr7.AFtest = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,SNP=chr7phase3.GSNP.bim$id,LOGP1df=pr1x,LOGP2df=pr2x)
chr7.AFtest = cbind.data.frame(chr7.AFtest, WALD1df.HWD=pr1x, WALD1df.HWE=pr2x)
write.csv(chr7.AFtest,file='chr7.AFtest.csv')



########################################################
# Two vertical Manhattan plots (First: p-value, Second: Allele Frequency Difference ) 
########################################################

# 1df WALD, assuming HWD (Combined RAF>=0.05)
m.data.temp = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,
                               SNP=chr7phase3.GSNP.bim$id,
                               LOGP=chr7.AFtest$WALD1df.HWD,
                               RMAF=chr7.SNP.FreqTable$F.M_RAF)
m.data.temp$SNP = as.character(m.data.temp$SNP)

lgt = length(which(m.data.temp$SNP=='.'))
snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
m.data.temp$SNP[which(m.data.temp$SNP=='.')] = snp.sp
m.data.temp$BP = m.data.temp$BP/1000000

m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
m.data.temp = m.data.temp[-which(chr7.SNP.FreqTable$A_RAF<0.05),]


tiff("chr7.aftest_pvalue_2Manhattan(RAF_0.05), 1dfwald_hwd.tiff",width=3000,height=5000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(2,1))
manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mb)',ylab=expression(-log[10](italic(p))),
                  cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="A",cex.main=cex,col="black",font=2,line=line)
manhattan1(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="RMAF", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mb)',ylim=c(-0.6,0.6),ylab='Female - Male MAF',
           cex.axis=1,cex.lab=1.2,cex=0.5,yaxt='n')
axis(side=2, at=c(-5:5)/10)
abline(h=0,col='red')
title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
dev.off()



# 1df WALD, assuming HWD (Combined 0.01<=RAF<0.05)
m.data.temp = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,
                               SNP=chr7phase3.GSNP.bim$id,
                               LOGP=chr7.AFtest$WALD1df.HWD,
                               RMAF=chr7.SNP.FreqTable$F.M_RAF)
m.data.temp$SNP = as.character(m.data.temp$SNP)

lgt = length(which(m.data.temp$SNP=='.'))
snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
m.data.temp$SNP[which(m.data.temp$SNP=='.')] = snp.sp
m.data.temp$BP = m.data.temp$BP/1000000

m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
m.data.temp = m.data.temp[which(chr7.SNP.FreqTable$A_RAF<0.05 & chr7.SNP.FreqTable$A_RAF>=0.01),]


tiff("chr7.aftest_pvalue_2Manhattan(0.01_RAF_0.05), 1dfwald_hwd.tiff",width=3000,height=5000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(2,1))
manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mbp)',ylab=expression(-log[10](italic(p))),
                  cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="A",cex.main=cex,col="black",font=2,line=line)
manhattan1(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="RMAF", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mbp)',ylim=c(-0.5,0.5),ylab='Female - Male RAF',
           cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=0,col='red')
title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
dev.off()



########################################################
# Histogram of p-value 
########################################################

# 1df WALD, assuming HWD (Combined RAF>=0.05)
m.data.temp = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,
                               SNP=chr7phase3.GSNP.bim$id,
                               LOGP=chr7.AFtest$WALD1df.HWD,
                               RMAF=chr7.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[-which(chr7.SNP.FreqTable$A_RAF<0.05),]
pvalues = 10^(-m.data.temp$LOGP)

tiff("chr7.aftest_pvalue_histogram(RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
datax = cbind.data.frame(pvalue=pvalues)
ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
dev.off()



# 1df WALD, assuming HWD (Combined 0.01<=RAF<0.05)
m.data.temp = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,
                               SNP=chr7phase3.GSNP.bim$id,
                               LOGP=chr7.AFtest$WALD1df.HWD,
                               RMAF=chr7.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[which(chr7.SNP.FreqTable$A_RAF<0.05 & chr7.SNP.FreqTable$A_RAF>=0.01),]
pvalues = 10^(-m.data.temp$LOGP)

tiff("chr7.aftest_pvalue_histogram(0.01_RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
datax = cbind.data.frame(pvalue=pvalues)
ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
dev.off()



########################################################
# QQ-plots of p-value 
########################################################

# 1df WALD, assuming HWD (Combined RAF>=0.05)
m.data.temp = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,
                               SNP=chr7phase3.GSNP.bim$id,
                               LOGP=chr7.AFtest$WALD1df.HWD,
                               RMAF=chr7.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[-which(chr7.SNP.FreqTable$A_RAF<0.05),]
pvalues = 10^(-m.data.temp$LOGP)

tiff("chr7.aftest_pvalue_QQ(RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
qqMargin(pvalues)
dev.off()



# 1df WALD, assuming HWD (Combined 0.01<=RAF<0.05)
m.data.temp = cbind.data.frame(CHR=7,BP=chr7phase3.GSNP.bim$pos,
                               SNP=chr7phase3.GSNP.bim$id,
                               LOGP=chr7.AFtest$WALD1df.HWD,
                               RMAF=chr7.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[which(chr7.SNP.FreqTable$A_RAF<0.05 & chr7.SNP.FreqTable$A_RAF>=0.01),]
pvalues = 10^(-m.data.temp$LOGP)

tiff("chr7.aftest_pvalue_QQ(0.01_RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
qqMargin(pvalues)
dev.off()




