####################################################################################
# Form SNP table specifying basic ethnicity-stratified results
# Will need .GT.FORMAT file pre-obtained using vcftools commands,
# Build GRCh37, chrX.
####################################################################################

# Find extreme SNPs
chrX.SNP.FreqTable = read.csv('chrX.SNP.FreqTable.csv')
pos1 = chrX.SNP.FreqTable$POS[which(chrX.SNP.FreqTable$POS<=2699520)][order(abs(
       chrX.SNP.FreqTable$F.M_RAF[which(chrX.SNP.FreqTable$POS<=2699520)]), decreasing=TRUE)][1:2]
pos2 = chrX.SNP.FreqTable$POS[which(chrX.SNP.FreqTable$POS>=154931044)][order(abs(
       chrX.SNP.FreqTable$F.M_RAF[which(chrX.SNP.FreqTable$POS>=154931044)]), decreasing=TRUE)][1:2]
pos3 = chrX.SNP.FreqTable$POS[which(chrX.SNP.FreqTable$POS>=88400000 & chrX.SNP.FreqTable$POS<=92000000)][order(abs(
       chrX.SNP.FreqTable$F.M_RAF[which(chrX.SNP.FreqTable$POS>=88400000 & chrX.SNP.FreqTable$POS<=92000000)]), decreasing=TRUE)][1:2]
pos4 = chrX.SNP.FreqTable$POS[which(!(chrX.SNP.FreqTable$POS>=88400000 & chrX.SNP.FreqTable$POS<=92000000) & !chrX.SNP.FreqTable$POS>=154931044 & !chrX.SNP.FreqTable$POS<=2699520)][order(abs(
       chrX.SNP.FreqTable$F.M_RAF[which(!(chrX.SNP.FreqTable$POS>=88400000 & chrX.SNP.FreqTable$POS<=92000000) & !chrX.SNP.FreqTable$POS>=154931044 & !chrX.SNP.FreqTable$POS<=2699520)]), decreasing=TRUE)][1:2]
chrX.extmSNP = cbind.data.frame(CHR=NA,POS=c(pos1,pos2,pos3,pos4))
chrX.extmSNP$CHR = 'X'          # Need modification

write.table(chrX.extmSNP, file = "chrX.extmSNP.txt", 
            append = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Read .FORMAT file
chrX.extmSNP = read.table('data/chrX/chrX_extmSNP.GT.FORMAT',header=T)

# Read sex info
sex.info = read.table('data/chrX/sex_info3.txt',header=T)

# Read ethnicity info
ethnicity.info = read.csv('data/chrX/sample_ethnicity.csv',header=T)

# SNP table (Xchr)
chrX.extmSNP.table = as.data.frame(matrix(NA,18*8,17))
colnames(chrX.extmSNP.table) = c('SNP','CHR','POS','REGION','A1','A2','Superpopulation',
                                 'SEX','A1A1','A1A2','A2A2','A1','A2','MAF','HWD.delta',
                                 'HWE.p','MAF.p')
sp.group = c('ALL','EAS','EUR','AFR','AMR','SAS')
sex.group = c('F','M','Both')
samples = colnames(chrX.extmSNP)[3:ncol(chrX.extmSNP)]
samples.eth = rep(NA,length(samples))
samples.sex = rep(NA,length(samples))
for(i in 1:length(samples))
{
  samples.eth[i] = ethnicity.info$Super.Population[which(ethnicity.info$Sample==samples[i])]
  samples.sex[i] = sex.info$sexinfo[which(sex.info$IID==samples[i])]
}

ind = 0
for(i in 1:nrow(chrX.extmSNP))
{
  gt = chrX.extmSNP[i,3:ncol(chrX.extmSNP)]
  for(sp in sp.group)
  {
    a1a1.set = a1a2.set = a2a2.set = c() 
    for(sx in sex.group)
    {
      ind = ind + 1
      a1a1 = 0; a1a2 = 0; a2a2 = 0
      for(k in 1:length(gt))
      {
        if((samples.eth[k]==sp|sp=='ALL') & (samples.sex[k]==sx|sx=='Both'))
        {
          if(gt[k]=='0|0'|gt[k]=='0') a1a1 = a1a1+1
          else if(gt[k]=='0|1'|gt[k]=='1|0') a1a2 = a1a2+1
          else if(gt[k]=='1|1'|gt[k]=='1') a2a2 = a2a2+1 
        }
      }
      a1a1.set = c(a1a1.set, a1a1)
      a1a2.set = c(a1a2.set, a1a2)
      a2a2.set = c(a2a2.set, a2a2)
      MAF = (2*a1a1+a1a2)/(2*a1a1+2*a1a2+2*a2a2)
      pos = chrX.extmSNP$POS[i]
      chrX.extmSNP.table$SNP[ind] = as.character(chrX.SNP.FreqTable$ID[which(chrX.SNP.FreqTable$POS==pos)])
      chrX.extmSNP.table$CHR[ind] = 'X'                # Need modification
      chrX.extmSNP.table$POS[ind] = pos
      
      if(pos>=88400000 & pos<=92000000) chrX.extmSNP.table$REGION[ind] = 'PAR3'
      else if(pos>=60001 & pos<=2699250) chrX.extmSNP.table$REGION[ind] = 'PAR1'
      else if(pos>=154931044 & pos<=155260560) chrX.extmSNP.table$REGION[ind] = 'PAR2'
      else chrX.extmSNP.table$REGION[ind] = 'NPR'
      rg = chrX.extmSNP.table$REGION[ind]
      
      chrX.extmSNP.table[ind,5:6] = as.character(chrX.SNP.FreqTable[which(chrX.SNP.FreqTable$POS==chrX.extmSNP$POS[i]),5:6])
      chrX.extmSNP.table[ind,7:8] = c(sp,sx)
      
      if((rg=='NPR'|rg=='PAR3') & sx=='M') 
      {
        chrX.extmSNP.table[ind,9:14] = c(a1a1,NA,a2a2,a1a1,a2a2,MAF)
        chrX.extmSNP.table$HWD.delta[ind] = chrX.extmSNP.table$HWE.p[ind] = chrX.extmSNP.table$`MAF.p`[ind] = NA
      }
      else if((rg=='NPR'|rg=='PAR3') & sx=='Both') 
      {
        chrX.extmSNP.table[ind,9:13] = c(NA,NA,NA,
                                         2*a1a1.set[1]+a1a2.set[1]+a1a1.set[2],2*a2a2.set[1]+a1a2.set[1]+a2a2.set[2])
        chrX.extmSNP.table[ind,14] = chrX.extmSNP.table[ind,12]/(chrX.extmSNP.table[ind,12]+chrX.extmSNP.table[ind,13])
        chrX.extmSNP.table$HWD.delta[ind] = chrX.extmSNP.table$HWE.p[ind] = NA
        chrX.extmSNP.table$`MAF.p`[ind] = 10^(-wald.1df.hwd.xchr(c(a1a1.set[1],a1a2.set[1],a2a2.set[1],
                                                                 a1a1.set[2],a1a2.set[2],a2a2.set[2])))
      }
      else 
      {
        chrX.extmSNP.table[ind,9:14] = c(a1a1,a1a2,a2a2,2*a1a1+a1a2,2*a2a2+a1a2,MAF)
        if(sp!='ALL')
        {
          chrX.extmSNP.table$HWD.delta[ind] = a1a1/(a1a1+a1a2+a2a2) - MAF^2
          chrX.extmSNP.table$HWE.p[ind] = hwe(c(a1a1,a1a2,a2a2))
        }else{  chrX.extmSNP.table$HWD.delta[ind] = chrX.extmSNP.table$HWE.p[ind] = NA  }
        
        if(sx=='Both')  chrX.extmSNP.table$`MAF.p`[ind] = 10^(-wald.1df.hwd.auto(c(a1a1.set[1],a1a2.set[1],a2a2.set[1],
                                                                                      a1a1.set[2],a1a2.set[2],a2a2.set[2])))
        else chrX.extmSNP.table$`MAF.p`[ind] = NA
      }
    }
  }
}


for(i in 1:nrow(chrX.extmSNP.table))
{
  k = ceiling(i/18)
  if(chrX.extmSNP.table$MAF[18*(k-1)+3]>0.5) chrX.extmSNP.table[i,5:6] = rev(chrX.extmSNP.table[i,5:6])
}


chrX.extmSNP.table$A1[which(chrX.extmSNP.table$A1==1)] = 'A'
chrX.extmSNP.table$A1[which(chrX.extmSNP.table$A1==2)] = 'C'
chrX.extmSNP.table$A1[which(chrX.extmSNP.table$A1==3)] = 'G'
chrX.extmSNP.table$A1[which(chrX.extmSNP.table$A1==4)] = 'T'
chrX.extmSNP.table$A2[which(chrX.extmSNP.table$A2==1)] = 'A'
chrX.extmSNP.table$A2[which(chrX.extmSNP.table$A2==2)] = 'C'
chrX.extmSNP.table$A2[which(chrX.extmSNP.table$A2==3)] = 'G'
chrX.extmSNP.table$A2[which(chrX.extmSNP.table$A2==4)] = 'T'










