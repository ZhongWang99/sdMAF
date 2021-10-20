####################################################################################
# This file contains the process of running sdMAF testing between male and female,
#     which will generate a .csv file.
# Detailed plotting procedure also included.
# Build GRCh37, chrX.
####################################################################################

library(qqman)
library(genio)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(cowplot)

chrX.SNP.FreqTable = read.csv('chrX.SNP.FreqTable.csv')

xphase3.GSNP.bim = read_bim('data/chrX/xphase3_with_sex.XY.GSNP.bim')
xphase3.GSNP.bim$chr = as.numeric(xphase3.GSNP.bim$chr)


# Initialize dataset & sdMAF test
ind.temp = grep("F_A1A1", colnames(chrX.SNP.FreqTable))
afdiff.m1.dt = as.matrix(chrX.SNP.FreqTable[which(xphase3.GSNP.bim$chr==23),ind.temp:(ind.temp+5)])    # NPR and PAR3, 'Xchr'
afdiff.m2.dt = as.matrix(chrX.SNP.FreqTable[which(xphase3.GSNP.bim$chr!=23),ind.temp:(ind.temp+5)])    # PAR1 and PAR2, 'Auto'

pr1x = apply(afdiff.m1.dt,MARGIN = 1,FUN = wald.1df.hwd.xchr)
pr2x = apply(afdiff.m1.dt,MARGIN = 1,FUN = wald.1df.hwe.xchr)
pr3x = apply(afdiff.m2.dt,MARGIN = 1,FUN = wald.1df.hwd.auto)
pr4x = apply(afdiff.m2.dt,MARGIN = 1,FUN = wald.1df.hwe.auto)


# Record Data of testing results
chrX.AFtest = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],LOGP1df=c(pr1x,pr2x),
                               LOGP2df=c(rep(NA,length(pr1x)),pr3x))
chrX.AFtest = cbind.data.frame(chrX.AFtest, WALD1df.HWD=c(pr1x,pr3x), WALD1df.HWE=c(pr2x,pr4x))
write.csv(chrX.AFtest,file='chrX.AFtest.csv')


#####################################################################
# Vertical Manhattan plots (First: p-value, Second: Allele 
#     Frequency Difference ) 
#####################################################################

# 1df WALD, assuming HWD (Combined MAF>=0.05)
m.data.temp = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               LOGP=chrX.AFtest$WALD1df.HWD,
                               RMAF=chrX.SNP.FreqTable$F.M_RAF)
m.data.temp$SNP = as.character(m.data.temp$SNP)

lgt = length(which(m.data.temp$SNP=='.'))
snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
m.data.temp$SNP[which(m.data.temp$SNP=='.')] = snp.sp
snphighlight = m.data.temp$SNP[which((m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)|(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)|
                                       (m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560))]
m.data.temp$BP = m.data.temp$BP/1000000

m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
m.data.temp = m.data.temp[-which(chrX.SNP.FreqTable$A_RAF<0.05),]

tiff("chrX.aftest_pvalue_2Manhattan(MAF_0.05), 1dfwald_hwd.tiff",width=3000,height=5000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(2,1))
manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylab=expression(-log[10](italic(p))),
                  highlight = snphighlight,xaxt="n",
                  cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="A",cex.main=cex,col="black",font=2,line=line)
manhattan1(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="RMAF", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylim=c(-0.6,0.6),ylab='Female - Male MAF',
           highlight = snphighlight, xaxt="n", yaxt="n",
           cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
axis(2, at = seq(-0.5, 0.5, by = 0.1), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=0,col='red')
text(c(2,90,155), c(0.45,0.45,0.45), c("PAR1","PAR3","PAR2"),
     cex = 1, col='grey')
title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
dev.off()


# 1df WALD, assuming HWD (Combined 0.01<=RAF<0.05)
m.data.temp = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               LOGP=chrX.AFtest$WALD1df.HWD,
                               RMAF=chrX.SNP.FreqTable$F.M_RAF)
m.data.temp$SNP = as.character(m.data.temp$SNP)

lgt = length(which(m.data.temp$SNP=='.'))
snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
m.data.temp$SNP[which(m.data.temp$SNP=='.')] = snp.sp
snphighlight = m.data.temp$SNP[which((m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)|(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)|
                                       (m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560))]
m.data.temp$BP = m.data.temp$BP/1000000

m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
m.data.temp = m.data.temp[which(chrX.SNP.FreqTable$A_RAF<0.05 & chrX.SNP.FreqTable$A_RAF>=0.01),]

tiff("chrX.aftest_pvalue_2Manhattan(0.01_RAF_0.05), 1dfwald_hwd.tiff",width=3000,height=5000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(2,1))
manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylab=expression(-log[10](italic(p))),
                  highlight = snphighlight, xaxt="n",
                  cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="A",cex.main=cex,col="black",font=2,line=line)
manhattan1(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="RMAF", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylim=c(-0.5,0.5),ylab='Female - Male RAF',
           highlight = snphighlight, xaxt="n",
           cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=0,col='red')
title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
dev.off()



#####################################################################
# Histogram of p-value (4 subgraphs: NPR,PAR1,PAR2,PAR3) 
#####################################################################

# 1df WALD, assuming HWD (Combined RAF>=0.05)
m.data.temp = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               LOGP=chrX.AFtest$WALD1df.HWD,
                               RMAF=chrX.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[-which(chrX.SNP.FreqTable$A_RAF<0.05),]
pvalues = 10^(-m.data.temp$LOGP)

pvalues1 = pvalues[-which((m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)|(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)|
                            (m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560))]
datax = cbind.data.frame(pvalue=pvalues1)
g1 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g1 = ggdraw(add_sub(g1, "NPR", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues2 = pvalues[which(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)]
datax = cbind.data.frame(pvalue=pvalues2)
g2 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g2 = ggdraw(add_sub(g2, "PAR1", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues3 = pvalues[which(m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560)]
datax = cbind.data.frame(pvalue=pvalues3)
g3 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g3 = ggdraw(add_sub(g3, "PAR2", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues4 = pvalues[which(m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)]
datax = cbind.data.frame(pvalue=pvalues4)
g4 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g4 = ggdraw(add_sub(g4, "PAR3", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

p = plot_grid(g1,g2,g3,g4,labels=c('A','B','C','D'),ncol=2,nrow=2)

tiff("chrX.aftest_pvalue_histogram(MAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
p
dev.off()



# 1df WALD, assuming HWD (Combined 0.01<=RAF<0.05)
m.data.temp = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               LOGP=chrX.AFtest$WALD1df.HWD,
                               RMAF=chrX.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[which(chrX.SNP.FreqTable$A_RAF<0.05 & chrX.SNP.FreqTable$A_RAF>=0.01),]
pvalues = 10^(-m.data.temp$LOGP)

pvalues1 = pvalues[-which((m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)|(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)|
                            (m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560))]
datax = cbind.data.frame(pvalue=pvalues1)
g1 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g1 = ggdraw(add_sub(g1, "NPR", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues2 = pvalues[which(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)]
datax = cbind.data.frame(pvalue=pvalues2)
g2 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g2 = ggdraw(add_sub(g2, "PAR1", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues3 = pvalues[which(m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560)]
datax = cbind.data.frame(pvalue=pvalues3)
g3 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g3 = ggdraw(add_sub(g3, "PAR2", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues4 = pvalues[which(m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)]
datax = cbind.data.frame(pvalue=pvalues4)
g4 = ggplot(datax, aes(x=pvalue)) + geom_histogram(color="black", fill="grey", boundary=0) + xlab('p value')
g4 = ggdraw(add_sub(g4, "PAR3", size=8, vpadding=grid::unit(0, "lines"),y = 38.5, x = 0.03, hjust = 0, fontface = "bold"))

p = plot_grid(g1,g2,g3,g4,labels=c('A','B','C','D'),ncol=2,nrow=2)

tiff("chrX.aftest_pvalue_histogram(0.01_RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
p
dev.off()



#####################################################################
# QQ-plots of p-value (4 subgraphs: NPR,PAR1,PAR2,PAR3) using 
#     pre-defined qqMargin function (qqMargin.R)
#####################################################################

# 1df WALD, assuming HWD (Combined RAF>=0.05)
m.data.temp = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               LOGP=chrX.AFtest$WALD1df.HWD,
                               RMAF=chrX.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[-which(chrX.SNP.FreqTable$A_RAF<0.05),]
pvalues = 10^(-m.data.temp$LOGP)

pvalues1 = pvalues[-which((m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)|(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)|
                            (m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560))]
g1 = qqMargin(pvalues1, x_breaks=seq(0,6,by=2), y_breaks=seq(0,300,by=100), x_lim = c(0,6), y_lim = c(0,300))
g1 = ggdraw(add_sub(g1, "NPR", size=8, vpadding=grid::unit(0, "lines"),y = 37, x = 0.03, hjust = 0, fontface = "bold"))

pvalues2 = pvalues[which(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)]
g2 = qqMargin(pvalues2, x_breaks=seq(0,6,by=2), y_breaks=seq(0,300,by=100), x_lim = c(0,6), y_lim = c(0,300))
g2 = ggdraw(add_sub(g2, "PAR1", size=8, vpadding=grid::unit(0, "lines"),y = 37, x = 0.03, hjust = 0, fontface = "bold"))

pvalues3 = pvalues[which(m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560)]
g3 = qqMargin(pvalues3, x_breaks=seq(0,6,by=2), y_breaks=seq(0,300,by=100), x_lim = c(0,6), y_lim = c(0,300))
g3 = ggdraw(add_sub(g3, "PAR2", size=8, vpadding=grid::unit(0, "lines"),y = 37, x = 0.03, hjust = 0, fontface = "bold"))

pvalues4 = pvalues[which(m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)]
g4 = qqMargin(pvalues4, x_breaks=seq(0,6,by=2), y_breaks=seq(0,300,by=100), x_lim = c(0,6), y_lim = c(0,300))
g4 = ggdraw(add_sub(g4, "PAR3", size=8, vpadding=grid::unit(0, "lines"),y = 37, x = 0.03, hjust = 0, fontface = "bold"))

p = plot_grid(g1,g2,g3,g4,labels=c('A','B','C','D'),ncol=2,nrow=2)

tiff("chrX.aftest_pvalue_QQ(RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
p
dev.off()



# 1df WALD, assuming HWD (Combined 0.01<=RAF<0.05)
m.data.temp = cbind.data.frame(CHR=23,BP=xphase3.GSNP.bim$pos[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               SNP=xphase3.GSNP.bim$id[c(which(xphase3.GSNP.bim$chr==23),which(xphase3.GSNP.bim$chr!=23))],
                               LOGP=chrX.AFtest$WALD1df.HWD,
                               RMAF=chrX.SNP.FreqTable$F.M_RAF)
m.data.temp = m.data.temp[which(chrX.SNP.FreqTable$A_RAF<0.05 & chrX.SNP.FreqTable$A_RAF>=0.01),]
pvalues = 10^(-m.data.temp$LOGP)

pvalues1 = pvalues[-which((m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)|(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)|
                            (m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560))]
g1 = qqMargin(pvalues1)
g1 = ggdraw(add_sub(g1, "NPR", size=8, vpadding=grid::unit(0, "lines"),y = 33.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues2 = pvalues[which(m.data.temp$BP>=60001 & m.data.temp$BP<=2699520)]
g2 = qqMargin(pvalues2)
g2 = ggdraw(add_sub(g2, "PAR1", size=8, vpadding=grid::unit(0, "lines"),y = 33.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues3 = pvalues[which(m.data.temp$BP>=154931044 & m.data.temp$BP<=155260560)]
g3 = qqMargin(pvalues3)
g3 = ggdraw(add_sub(g3, "PAR2", size=8, vpadding=grid::unit(0, "lines"),y = 33.5, x = 0.03, hjust = 0, fontface = "bold"))

pvalues4 = pvalues[which(m.data.temp$BP>=88400000 & m.data.temp$BP<=92000000)]
g4 = qqMargin(pvalues4)
g4 = ggdraw(add_sub(g4, "PAR3", size=8, vpadding=grid::unit(0, "lines"),y = 33.5, x = 0.03, hjust = 0, fontface = "bold"))

p = plot_grid(g1,g2,g3,g4,labels=c('A','B','C','D'),ncol=2,nrow=2)

tiff("chrX.aftest_pvalue_QQ(0.01_RAF_0.05), 1dfwald_hwd.tiff",width=2500,height=1980,res=300)
p
dev.off()









