# This code is just trying to visialize principal component analysis to find 
# How closely correlated different traits are
# The goal is model a trait against another rather than to days after planting

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
DW$TOTALPOD <- rowSums(DW[,c('DWAP', 'DWGP', 'DWPOD', 'DWMP')], na.rm=T) 

fhDW <- DW[DW$HNO==11,c(1,3:11,17, 19)]

dat <- merge(seed[,c('plot', 'geno', 'X100SD', 'site')], fhDW, all=TRUE, by=c('plot'))

all <- na.omit(dat[,c(2:5,8,13,15)])

meanAll <- aggregate(all[,c(2,4:7)], by=list(all$geno, all$site), na.rm=T, FUN=mean)
colnames(meanAll)[1:2] <- c('geno', 'site')

dat.pcs <- prcomp(meanAll[,c(3:7)], center = TRUE, scale = TRUE)
summary(dat.pcs)

biplot <- ggbiplot(dat.pcs)
biplotLabels <- ggbiplot(dat.pcs, labels = meanAll$geno)

setwd(paste(wd,'out', sep='/'))
pdf(file='pcsPlot.pdf', onefile = T, width=10, height=10)
print(biplot)
print(biplotLabels)

dev.off()
