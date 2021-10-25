####################################################################################
# Form complete SNP table specifying basic SNP-related information (including 
#     positions, stratified allele counts, stratified allele frequency, etc).
# Will need .bim file, .frqx file and .frq file pre-obtained using PLINK commands,
#     the order of SNPs should be the same across files.
# Build GRCh37, autosomal.
# Changing chromosome a to chromosome b by replacing 'chra' with 'chrb'
####################################################################################

library(genio)

chr7phase3.GSNP.bim = read_bim('data/chr7/chr7_SNP_only.recode.bim')
chr7phase3.GSNP.bim$id = as.character(chr7phase3.GSNP.bim$id)

gf.GSNP.male = read.table('data/chr7/chr7.maf.GSNP.male.frqx',sep='\t',header=T)    # Genotype frequency
gf.GSNP.female = read.table('data/chr7/chr7.maf.GSNP.female.frqx',sep='\t',header=T) 
gf.GSNP.all = read.table('data/chr7/chr7.maf.GSNP.all.frqx',sep='\t',header=T) 

raf.GSNP.all = read.table('data/chr7/chr7.maf.GSNP.all.frq',header=T) 
raf.GSNP.female = read.table('data/chr7/chr7.maf.GSNP.female.frq',header=T) 
raf.GSNP.male = read.table('data/chr7/chr7.maf.GSNP.male.frq',header=T) 

A_A1 = 2*gf.GSNP.all$C.HOM.A1. + gf.GSNP.all$C.HET.
A_A2 = 2*gf.GSNP.all$C.HOM.A2. + gf.GSNP.all$C.HET.

F_A1A1 = gf.GSNP.female$C.HOM.A1.
F_A1A2 = gf.GSNP.female$C.HET.
F_A2A2 = gf.GSNP.female$C.HOM.A2.
mis_index = which(gf.GSNP.female$A1!=gf.GSNP.all$A1)
index_val = F_A2A2[mis_index]
F_A2A2[mis_index] = F_A1A1[mis_index]
F_A1A1[mis_index] = index_val

M_A1A1 = gf.GSNP.male$C.HOM.A1.
M_A1A2 = gf.GSNP.male$C.HET.
M_A2A2 = gf.GSNP.male$C.HOM.A2.
mis_index = which(gf.GSNP.male$A1!=gf.GSNP.all$A1)
index_val = M_A2A2[mis_index]
M_A2A2[mis_index] = M_A1A1[mis_index]
M_A1A1[mis_index] = index_val

F_RAF = raf.GSNP.female$MAF; F_RAF[which(gf.GSNP.female$A1!=gf.GSNP.all$A1)] = 1 - F_RAF[which(gf.GSNP.female$A1!=gf.GSNP.all$A1)]
M_RAF = raf.GSNP.male$MAF; M_RAF[which(gf.GSNP.male$A1!=gf.GSNP.all$A1)] = 1 - M_RAF[which(gf.GSNP.male$A1!=gf.GSNP.all$A1)]

chr7.SNP.FreqTable = cbind.data.frame(ID=chr7phase3.GSNP.bim$id, POS=chr7phase3.GSNP.bim$pos,
                                 A1=raf.GSNP.all$A1, A2=raf.GSNP.all$A2, A_A1=A_A1, A_A2=A_A2, F_A1A1=F_A1A1, F_A1A2=F_A1A2,
                                 F_A2A2=F_A2A2, M_A1A1=M_A1A1, M_A1A2=M_A1A2, M_A2A2=M_A2A2, F_RAF=F_RAF,
                                 M_RAF=M_RAF, A_RAF=raf.GSNP.all$MAF, 'F-M_RAF'=F_RAF-M_RAF)

write.csv(chr7.SNP.FreqTable,file='chr7.SNP.FreqTable.csv')
rm(gf.GSNP.male,gf.GSNP.all,gf.GSNP.female,raf.GSNP.male,raf.GSNP.all,raf.GSNP.female)     # Remove temporary files to save space


chr7.SNP.FreqTable = read.csv('chr7.SNP.FreqTable.csv')













