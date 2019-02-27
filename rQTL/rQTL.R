rm(list=ls())
#install.packages("qtl")
# tutorial
# http://www.rqtl.org/tutorials/
library(qtl)
#?read.cross()
#ls("package:qtl")

# parameters
step_size = 5   # step size (cM)
cofactor_size = 3   # cofactor window (cM)
it_max=1   # iterations
LOD_thr=1.18 # LOD for p-value = 0.001
#LOD_thr=0.59 # LOD for p-value = 0.05

dir=getwd()
# source extra functions
setwd(paste(dir,'src',sep='/'))
files=list.files()
for(file in files){source(file)}
setwd(dir)
####### Make phenotype File for DIPI ##############
phenoREF <- read.csv(paste(dir,'data', 'phenotype.csv', sep='/'))
DIPI <- read.csv(paste(dir,'data', 'DIPI.csv', sep='/'), row.names = 1)

DIPI_geno <- aggregate(DIPI$slope, by= list(DIPI$site, DIPI$RIL), FUN=mean, na.rm=T)
colnames(DIPI_geno) <- c('site', 'RIL', 'DIPI')

DIPI_1 <- DIPI_geno[DIPI_geno$site==1,c('RIL', 'DIPI')]
DIPI_2 <- DIPI_geno[DIPI_geno$site==2, c('RIL','DIPI')]

DIPI_pheno <- merge(phenoREF, DIPI_1, all.x=T, by='RIL')
DIPI_pheno <- merge(DIPI_pheno, DIPI_2, all.x=T, by='RIL')

t.test(DIPI_pheno$DIPI.x, DIPI_pheno$DIPI.y)
#write.csv(DIPI_pheno, file=paste(dir, 'data', 'DIPI_phenotype.csv',sep='/'))

################################
### Reading data & diagnosis ###
################################
data = read.cross(dir = paste(dir,'data',sep='/'),
                  format = "csvs",
                  genfile = "genotype_AB.csv",
                  phefile = 'DIPI_phenotype.csv',
                  F.gen = 11,
                  na.strings=c('NA'))
## Make sure that CSVs are not encoded in UTF-8-BOM, this could produce problems
# convert to RI lines
data=convert2riself(data)
# extract trait and genotpe id
data$pheno=data$pheno[,c('DIPI.y','RIL')]
summary(data)
#plot(data)
head(data$pheno)

# Drop markers with no genotypes
#data=drop.nullmarkers(data)

# estimate pair-wise recombination factor
data = est.rf(data)
#plotRF(data)

# estimate new map
#newmap <- est.map(data, error.prob=0.01)
#plotMap(data, newmap)

# if you were to replace map
#data = replace.map(data,newmap)

# compute genotype error
#data = calc.errorlod(data,
#                     error.prob = .001)
#top.errorlod(data)
#plotGeno(data)

###################
### QTL mapping ###
###################

# compute genotype probability
data = calc.genoprob(data, step = step_size)

# compute imputations
data = sim.geno(data, step = step_size, n.draws = 100)

# Simple Interval Mapping (SIM) with Haley-Knott regression
SIM_hk=scanone(cross = data, method = 'hk', n.cluster = 4)
# Threshold SIM_hk results with LOD_thr
SIM=threshold(SIM_hk, LOD_thr)

#########
# This CIM still needs some work, because the cim from the qtl package is a little weird...
#########
# Composite Interval Maping (CIM) 
# with Haley-Knott regression
# 1 marker covariable
# step size of 5 cM
for(it in 1:it_max){
  if(it == 1){
    CIM_hk = cim(cross = data,
                  n.marcovar=1, 
                  window = cofactor_size, 
                  method = 'hk')
  }else{
    CIM_hk = CIM_hk + cim(cross = data,
                  n.marcovar=1, 
                  window = cofactor_size, 
                  method = 'hk')
  }
}
# compute mean LOD value
CIM_hk$lod=CIM_hk$lod/it_max

SIM=threshold(SIM_hk,LOD_thr)
CIM=threshold(CIM_hk,LOD_thr)

# plot CIM results with LOD = 3 threshold line
plot(SIM_hk,
     CIM_hk,
     SIM,
     col=c("black",'blue','red'),
     type=c('l','l','p'),
     ylim=c(0,ceiling(max(SIM_hk$lod)*1.5)),
     main = 'LOD score from CIM: SIM values in black, CIM values in blue'); abline(h=LOD_thr, col='red')

# multiple QTL model with markers selected form CIM
# chromosome #s from CIM
(chr = as.numeric(SIM$chr))
# positional #s from CIM
(pos = as.numeric(SIM$pos))
# generate QTL object
(QTL = makeqtl(data, chr, pos))

# Define formula as only additive QTL effects
formula=paste('y~',
              paste(
                paste('Q',
                      c(1:nrow(SIM)),
                      sep = ''),
                collapse = '+'),
              sep='')

# Fit multi-QTL model with initial QTLs
QTL_model = fitqtl(data, qtl = QTL,formula = formula, get.ests = T, method = 'hk')
summary(QTL_model)

# make small adjustments to QLT positions to test for better fits
refQTL <- refineqtl(data, qtl=QTL, formula = formula)
# enure refQTL are on actual markers by finding closes real marker
refQTL=find.markerpos(cross = data,
                      marker = find.marker(cross = data,
                                           chr = refQTL$chr, 
                                           pos = refQTL$pos))
# make QTL object
(refQTL=makeqtl(cross = data,chr = refQTL$chr, pos = refQTL$pos, qtl.name = rownames(refQTL)))

# Fit multi-QTL model with  refQTLs
refQTL_model = fitqtl(data, qtl = refQTL,formula = formula, get.ests = T, method = 'hk')
summary(refQTL_model)
