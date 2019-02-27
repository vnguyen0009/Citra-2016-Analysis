######## CODE OBTAINED FROM Christopher Hwang M.S.
### University of Florida Ph.D. Program 
###College of Engineering | Agricultural Biological Engineering


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
geno <- read.csv(paste(dir,'data', 'genotype_AB.csv', sep='/'))
DIPI <- read.csv(paste(dir,'data', 'DIPI.csv', sep='/'), row.names = 1)
genoList <- as.data.frame(geno$RIL[-c(1:2)]) #First two rows are empty
colnames(genoList)[1] <- c('RIL')


DIPI_geno <- aggregate(DIPI$slope, by= list(DIPI$site, DIPI$RIL), FUN=mean, na.rm=T)
colnames(DIPI_geno) <- c('site', 'RIL', 'DIPI')

DIPI_1 <- DIPI_geno[DIPI_geno$site==1,c('RIL', 'DIPI')]
DIPI_2 <- DIPI_geno[DIPI_geno$site==2, c('RIL','DIPI')]

DIPI_pheno <- merge(genoList, DIPI_1, all.x=T, by='RIL')
DIPI_pheno <- merge(DIPI_pheno, DIPI_2, all.x=T, by='RIL')

mean(DIPI_pheno$DIPI.x, na.rm=T)
var(DIPI_pheno$DIPI.x, na.rm=T)

mean(DIPI_pheno$DIPI.y, na.rm=T)
var(DIPI_pheno$DIPI.y, na.rm=T)

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

# compute genotype probability
data = calc.genoprob(data, step = step_size)

# compute imputations
data = sim.geno(data, step = step_size, n.draws = 100)

# Simple Interval Mapping (SIM) with Haley-Knott regression
SIM_hk=scanone(cross = data, method = 'hk', n.cluster = 4)
# Threshold SIM_hk results with LOD_thr
SIM=threshold(SIM_hk, LOD_thr)

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
capture.output(summary(QTL_model), file = paste(dir,'out', 'QTL_model.txt', sep='/'))

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
capture.output(summary(refQTL_model), file = paste(dir,'out', 'refQTL_model.txt',sep='/'))
