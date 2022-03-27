######################################################
# Perform the similar data analysis to gnomAD data
#   Sample size of each SNP is not fixed
#   Use .INFO data generated from original gnomAD v3.0
######################################################

library(qqman)
library(genio)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(cowplot)
library(liftOver)

chrX_gnomad = read.table("gnomad_nfe.INFO",header=T,stringsAsFactors = F)

invalid_ind = c(which(as.numeric(chrX_gnomad$AN_nfe)<=50000),which(is.na(as.numeric(as.character(chrX_gnomad$AF_nfe)))),
                which(is.na(as.numeric(as.character(chrX_gnomad$AF_nfe_XX)))),which(is.na(as.numeric(as.character(chrX_gnomad$AF_nfe_XY)))))
psedo = which(chrX_gnomad$POS[-invalid_ind]<=2781479|chrX_gnomad$POS[-invalid_ind]>=155701383)
chrX_gnomad$AF_nfe = as.numeric(as.character(chrX_gnomad$AF_nfe))
chrX_gnomad$AF_nfe_XX = as.numeric(as.character(chrX_gnomad$AF_nfe_XX))
chrX_gnomad$AF_nfe_XY = as.numeric(as.character(chrX_gnomad$AF_nfe_XY))
# 2781479|m.data.temp$BP>=155701383

# correct one
afdiff.m1.dt = as.matrix(data.frame(aaF=as.numeric(chrX_gnomad$nhomalt_nfe_XX[-invalid_ind][-psedo])-as.numeric(chrX_gnomad$AC_nfe_XX[-invalid_ind][-psedo])+as.numeric(chrX_gnomad$AN_nfe_XX[-invalid_ind][-psedo])/2,
                                    AaF=as.numeric(chrX_gnomad$AC_nfe_XX[-invalid_ind][-psedo])-as.numeric(chrX_gnomad$nhomalt_nfe_XX[-invalid_ind][-psedo])*2,
                                    AAF=as.numeric(chrX_gnomad$nhomalt_nfe_XX[-invalid_ind][-psedo]),
                                    aM=as.numeric(chrX_gnomad$AN_nfe_XY[-invalid_ind][-psedo])-as.numeric(chrX_gnomad$AC_nfe_XY[-invalid_ind][-psedo]),
                                    AaM=0,
                                    AM=as.numeric(chrX_gnomad$AC_nfe_XY[-invalid_ind][-psedo])))

afdiff.m2.dt = as.matrix(data.frame(aaF=as.numeric(chrX_gnomad$nhomalt_nfe_XX[-invalid_ind][psedo])-as.numeric(chrX_gnomad$AC_nfe_XX[-invalid_ind][psedo])+as.numeric(chrX_gnomad$AN_nfe_XX[-invalid_ind][psedo])/2,
                                    AaF=as.numeric(chrX_gnomad$AC_nfe_XX[-invalid_ind][psedo])-as.numeric(chrX_gnomad$nhomalt_nfe_XX[-invalid_ind][psedo])*2,
                                    AAF=as.numeric(chrX_gnomad$nhomalt_nfe_XX[-invalid_ind][psedo]),
                                    aaM=as.numeric(chrX_gnomad$nhomalt_nfe_XY[-invalid_ind][psedo])-as.numeric(chrX_gnomad$AC_nfe_XY[-invalid_ind][psedo])+as.numeric(chrX_gnomad$AN_nfe_XY[-invalid_ind][psedo])/2,
                                    AaM=as.numeric(chrX_gnomad$AC_nfe_XY[-invalid_ind][psedo])-as.numeric(chrX_gnomad$nhomalt_nfe_XY[-invalid_ind][psedo])*2,
                                    AAM=as.numeric(chrX_gnomad$nhomalt_nfe_XY[-invalid_ind][psedo])))


t1 = Sys.time()
pr1x = apply(afdiff.m1.dt,MARGIN = 1,FUN = wald.1df.hwd.xchr)
pr2x = apply(afdiff.m1.dt,MARGIN = 1,FUN = wald.1df.hwe.xchr)
pr3x = apply(afdiff.m2.dt,MARGIN = 1,FUN = wald.1df.hwd.auto)
pr4x = apply(afdiff.m2.dt,MARGIN = 1,FUN = wald.1df.hwe.auto)
t2 = Sys.time()
t2 - t1


# Generate a SNP Table (A1 is not necessarily minor allele, A1 matches REF in original dataset)
#    the selection of A1 also relates to RAF
#    here RAF corresponds to A1
afdiff.m1.dt = as.data.frame(afdiff.m1.dt)
afdiff.m2.dt = as.data.frame(afdiff.m2.dt)
chrXgnomad.SNP.FreqTable = cbind.data.frame(CHR=23, POS=as.numeric(c(chrX_gnomad$POS[-invalid_ind][-psedo],chrX_gnomad$POS[-invalid_ind][psedo])),
                                            A1=c(chrX_gnomad$REF[-invalid_ind][-psedo],chrX_gnomad$REF[-invalid_ind][psedo]),
                                            A2=c(chrX_gnomad$ALT[-invalid_ind][-psedo],chrX_gnomad$ALT[-invalid_ind][psedo]),
                                            A_A1=c(chrX_gnomad$AN_nfe[-invalid_ind][-psedo],chrX_gnomad$AN_nfe[-invalid_ind][psedo])-
                                              c(chrX_gnomad$AC_nfe[-invalid_ind][-psedo],chrX_gnomad$AC_nfe[-invalid_ind][psedo]), 
                                            A_A2=c(chrX_gnomad$AC_nfe[-invalid_ind][-psedo],chrX_gnomad$AC_nfe[-invalid_ind][psedo]), 
                                            F_A1A1=c(afdiff.m1.dt$aaF,afdiff.m2.dt$aaF),
                                            F_A1A2=c(afdiff.m1.dt$AaF,afdiff.m2.dt$AaF),
                                            F_A2A2=c(afdiff.m1.dt$AAF,afdiff.m2.dt$AAF),
                                            M_A1A1.A1=c(afdiff.m1.dt$aM,afdiff.m2.dt$aaM),
                                            M_A1A2=c(afdiff.m1.dt$AaM,afdiff.m2.dt$AaM),
                                            M_A2A2.A2=c(afdiff.m1.dt$AM,afdiff.m2.dt$AAM),
                                            F_RAF=1-c(chrX_gnomad$AF_nfe_XX[-invalid_ind][-psedo],
                                                      chrX_gnomad$AF_nfe_XX[-invalid_ind][psedo]),
                                            M_RAF=1-c(chrX_gnomad$AF_nfe_XY[-invalid_ind][-psedo],
                                                      chrX_gnomad$AF_nfe_XY[-invalid_ind][psedo]), 
                                            A_RAF=1-c(chrX_gnomad$AF_nfe[-invalid_ind][-psedo],
                                                      chrX_gnomad$AF_nfe[-invalid_ind][psedo]), 
                                            'F-M_RAF'=-c(chrX_gnomad$AF_nfe_XX[-invalid_ind][-psedo]-
                                                           chrX_gnomad$AF_nfe_XY[-invalid_ind][-psedo],
                                                         chrX_gnomad$AF_nfe_XX[-invalid_ind][psedo]-
                                                           chrX_gnomad$AF_nfe_XY[-invalid_ind][psedo]))
write.csv(chrXgnomad.SNP.FreqTable,file='chrXgnomad.nfe.SNP.FreqTable.csv')


# Record Data
chrXgnomad.AFtest = cbind.data.frame(CHR=23,BP=as.numeric(c(chrX_gnomad$POS[-invalid_ind][-psedo],chrX_gnomad$POS[-invalid_ind][psedo])))
chrXgnomad.AFtest = cbind.data.frame(chrXgnomad.AFtest, WALD1df.HWD=c(pr1x,pr3x), WALD1df.HWE=c(pr2x,pr4x))
write.csv(chrXgnomad.AFtest,file='chrXgnomad.nfe.AFtest.csv')

########################################################
# Two vertical Manhattan plots (First: p-value, Second: Allele Frequency Difference ) Region Specific
########################################################

chrXgnomad.AFtest = read.csv('chrXgnomad.nfe.AFtest.csv')

# 1df WALD, assuming HWD (Combined MAF>=0.05)
#   here LOGP stays the same whether we choose A1 or A2 as the target
m.data.temp = cbind.data.frame(CHR=23,BP=as.numeric(chrXgnomad.AFtest$BP),
                               SNP="0",
                               AF = c(chrX_gnomad$AF_nfe[-invalid_ind][-psedo],chrX_gnomad$AF_nfe[-invalid_ind][psedo]),
                               LOGP = chrXgnomad.AFtest$WALD1df.HWD,
                               RMAF = c(chrX_gnomad$AF_nfe_XX[-invalid_ind][-psedo]-
                                          chrX_gnomad$AF_nfe_XY[-invalid_ind][-psedo],
                                        chrX_gnomad$AF_nfe_XX[-invalid_ind][psedo]-
                                          chrX_gnomad$AF_nfe_XY[-invalid_ind][psedo]))

snp.sp = as.character(sample(0:100000000,size=nrow(m.data.temp),replace = F))
m.data.temp$SNP = snp.sp
reverse_ind = which(m.data.temp$AF>0.5)          # Some SNPs have AF over 0.5
m.data.temp$AF[reverse_ind] = 1-m.data.temp$AF[reverse_ind]
m.data.temp$RMAF[reverse_ind] = -m.data.temp$RMAF[reverse_ind]

m.data.temp = m.data.temp[which(m.data.temp$AF>=0.05),]
snphighlight = m.data.temp$SNP[which(m.data.temp$BP<=2781479|m.data.temp$BP>=155701383|(m.data.temp$BP>=89145000 & m.data.temp$BP<=92745001))]
m.data.temp$BP = m.data.temp$BP/1000000

m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1

tiff("chrXgnomad.AFtest_pvalue_2Manhattan(RAF_0.05,nfe).tiff",width=3000,height=5000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(2,1))
manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylab=expression(-log[10](italic(p))),
                  highlight = snphighlight,xaxt="n",
                  ylim = c(min(na.omit(m.data.temp$LOGP)), 47982.36),
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