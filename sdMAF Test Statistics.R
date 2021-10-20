####################################################################################
# 'Wald' Type Test consider whether or nor include Hardy Weinberg Disequilibrium
# Testing sex difference in Minor Allele Frequency between male and female
####################################################################################


# 'Wald' type, 1 d.f. assuming HWD, Xchr 
wald.1df.hwd.xchr <- function(x)
  # 'Wald' type, 1 d.f. assuming HWD, Xchr 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s2; r = r0+r1+r2
  pM = s2/s; pF = (0.5*r1+r2)/r 
  pAA.F = r2/r
  delta.F = pAA.F-pF^2
  stat = (pM-pF)^2/(1/s*(pM*(1-pM))+1/(2*r)*(pF*(1-pF)+delta.F))
  -pchisq(stat,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}


# 'Wald' type, 1 d.f. assuming HWE, Xchr 
wald.1df.hwe.xchr <- function(x)
  # 'Wald' type, 1 d.f. assuming HWE, Xchr 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s2; r = r0+r1+r2
  pM = s2/s; pF = (0.5*r1+r2)/r
  stat = (pM-pF)^2/(1/s*(pM*(1-pM))+1/(2*r)*(pF*(1-pF)))
  -pchisq(stat,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}


# 'Wald' type, 1 d.f. assuming HWD, Autosomal
wald.1df.hwd.auto <- function(x)
  # 'Wald' type, 1 d.f. assuming HWD, Autosomal 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s1+s2; r = r0+r1+r2
  pM = (0.5*s1+s2)/s; pF = (0.5*r1+r2)/r
  pAA.M = s2/s; pAA.F = r2/r
  delta.M = pAA.M-pM^2; delta.F = pAA.F-pF^2
  stat = (pM-pF)^2/(1/(2*s)*(pM*(1-pM)+delta.M)+1/(2*r)*(pF*(1-pF)+delta.F))
  -pchisq(stat,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}


# 'Wald' type, 1 d.f. assuming HWE, Autosomal
wald.1df.hwe.auto <- function(x)
  # 'Wald' type, 1 d.f. assuming HWE, Autosomal
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s1+s2; r = r0+r1+r2
  pM = (0.5*s1+s2)/s; pF = (0.5*r1+r2)/r
  stat = (pM-pF)^2/(1/(2*s)*(pM*(1-pM))+1/(2*r)*(pF*(1-pF)))
  -pchisq(stat,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}












