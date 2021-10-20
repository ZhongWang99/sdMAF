####################################################################################
# Form complete SNP table specifying basic SNP-related information (including 
#     positions, stratified allele counts, stratified allele frequency, etc).
# Will need .bim file, .frqx file and .frq file pre-obtained using PLINK commands,
#     the order of SNPs should be the same across files.
# Build GRCh37, chrX.
####################################################################################

require(genio)

xphase3.GSNP.bim = read_bim('data/chrX/xphase3_with_sex.XY.GSNP.bim')
xphase3.GSNP.bim$id = as.character(xphase3.GSNP.bim$id)

gf.GSNP.male = read.table('data/chrX/maf.xphase3.XY.GSNP.male.frqx',sep='\t',header=T)    # Genotype count
gf.GSNP.female = read.table('data/chrX/maf.xphase3.XY.GSNP.female.frqx',sep='\t',header=T) 
gf.GSNP.all = read.table('data/chrX/maf.xphase3.XY.GSNP.all.frqx',sep='\t',header=T) 

raf.GSNP.all = read.table('data/chrX/maf.xphase3.XY.GSNP.all.frq',header=T)  # Genotype frequency
raf.GSNP.female = read.table('data/chrX/maf.xphase3.XY.GSNP.female.frq',header=T) 
raf.GSNP.male = read.table('data/chrX/maf.xphase3.XY.GSNP.male.frq',header=T) 

PAR = xphase3.GSNP.bim$pos
PAR[which(PAR<=2699520)] = 1
PAR[which(PAR>=154931044)] = 2
PAR[which(PAR>=88400000 & PAR<=92000000)] = 3
PAR[which(PAR!=1 & PAR!=2 & PAR!=3)] = -999

A_A1 = 2*gf.GSNP.all$C.HOM.A1. + gf.GSNP.all$C.HET. + gf.GSNP.all$C.HAP.A1.
A_A2 = 2*gf.GSNP.all$C.HOM.A2. + gf.GSNP.all$C.HET. + gf.GSNP.all$C.HAP.A2.

F_A1A1 = gf.GSNP.female$C.HOM.A1.
F_A1A2 = gf.GSNP.female$C.HET.
F_A2A2 = gf.GSNP.female$C.HOM.A2.
mis_index = which(gf.GSNP.female$A1!=gf.GSNP.all$A1)
index_val = F_A2A2[mis_index]
F_A2A2[mis_index] = F_A1A1[mis_index]
F_A1A1[mis_index] = index_val

M_A1A1.A1 = gf.GSNP.male$C.HOM.A1. + gf.GSNP.male$C.HAP.A1.
M_A1A2 = gf.GSNP.male$C.HET.
M_A2A2.A2 = gf.GSNP.male$C.HOM.A2. + gf.GSNP.male$C.HAP.A2.
mis_index = which(gf.GSNP.male$A1!=gf.GSNP.all$A1)
index_val = M_A2A2.A2[mis_index]
M_A2A2.A2[mis_index] = M_A1A1.A1[mis_index]
M_A1A1.A1[mis_index] = index_val

F_RAF = raf.GSNP.female$MAF; F_RAF[which(gf.GSNP.female$A1!=gf.GSNP.all$A1)] = 1 - F_RAF[which(gf.GSNP.female$A1!=gf.GSNP.all$A1)]
M_RAF = raf.GSNP.male$MAF; M_RAF[which(gf.GSNP.male$A1!=gf.GSNP.all$A1)] = 1 - M_RAF[which(gf.GSNP.male$A1!=gf.GSNP.all$A1)]

chrX.SNP.FreqTable = cbind.data.frame(ID=xphase3.GSNP.bim$id, POS=xphase3.GSNP.bim$pos, PAR=PAR,
                                 A1=raf.GSNP.all$A1, A2=raf.GSNP.all$A2, A_A1=A_A1, A_A2=A_A2, F_A1A1=F_A1A1, F_A1A2=F_A1A2,
                                 F_A2A2=F_A2A2, M_A1A1.A1=M_A1A1.A1, M_A1A2=M_A1A2, M_A2A2.A2=M_A2A2.A2, F_RAF=F_RAF,
                                 M_RAF=M_RAF, A_RAF=raf.GSNP.all$MAF, 'F-M_RAF'=F_RAF-M_RAF)

write.csv(chrX.SNP.FreqTable,file='chrX.SNP.FreqTable.csv')















