rm(list=ls())

library(devtools)
#install_github("vqv/ggbiplot")
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(stats)
library(ggbiplot)

wd <- getwd()

setwd(paste(wd,'data', sep='/'))

DW <- read.csv('DW 4 1 19.csv', row.names=1, na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)
seed <- read.csv('FH seed 9 11 18.csv', row.names=1, na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)

fhDW <- DW[DW$HNO==11,c(1,3:11,17)]

dat <- merge(seed[,c('plot', 'geno', 'X100SD', 'site')], fhDW, all=TRUE, by=c('plot'))

all <- na.omit(dat[,c(2,3,4,7,12)])

dat.pcs <- prcomp(na.omit(all[,c(2:5)]), center = TRUE, scale = TRUE)
summary(dat.pcs)

biplot <- ggbiplot(dat.pcs)
biplotLabels <- ggbiplot(dat.pcs, labels = all$geno)

setwd(paste(wd,'out', sep='/'))
pdf(file='pcsPlot.pdf', onefile = T, width=10, height=10)
print(biplot)
print(biplotLabels)

dev.off()
