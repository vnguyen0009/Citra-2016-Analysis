# AUTHOR: VY T. NGUYEN
# LAST EDITED: 8/17/2019
# This code is to clean the dry weight data. There are no dates, only harvest numbers
# data for one ID is in multiple rows and need to be compiled to one

rm(list=ls())

library(dplyr)
library(plyr)

wd <- getwd()

setwd(paste(wd,'data', sep='/'))

source(paste(wd,'src', 'Hdate.R', sep='/')) #function to assign harvest dates
source(paste(wd,'src', 'Pdate.R', sep='/')) #function to assign planting dates

dates <- read.csv('dates.csv',stringsAsFactors = F)
RAW <- read.csv('DW.csv', na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)
genoCode <- read.csv('plotToGeno.csv', stringsAsFactors = F) #id to genotype file
FHdates <- read.csv('FH dates.csv', stringsAsFactors = F) #final harvest dates

#FHdates$H11 <- as.Date(FHdates$H11,tryFormats=c('%m/%d/%Y'))
#dates[,c(2:12)] <- as.Date(dates[,c(2:12)],tryFormats=c('%m/%d/%Y'))

data <- aggregate(RAW[,c(3:11)], by=list(RAW$HNO, RAW$plot), FUN=mean, na.rm=T)
colnames(data)[1:2] <- c('HNO', 'plot')

#data[is.na(data$DWTOTALLEAF)==T,'DWTOTALLEAF'] <- data[is.na(data$DWTOTALLEAF)==T,'DWMSL'] + data[is.na(data$DWTOTALLEAF)==T,'DWBL']

data[is.na(data$DWTOTALLEAF)==T,'DWTOTALLEAF'] <- rowSums(data[is.na(data$DWTOTALLEAF)==T,c('DWMSL', 'DWBL')], na.rm=T)
data[is.na(data$WHOLE)==T,'WHOLE'] <- rowSums(data[is.na(data$WHOLE)==T,
                                                   c('DWTOTALLEAF', 'TOTALSTEM', 'DWAP', 'DWGP', 'DWPOD', 'DWMP')], na.rm=T)

##Categorize D-lines
data$Dlines <- 0
data[grep('D',data$plot),'Dlines'] <- 1 
data$plot <- gsub("[a-zA-Z ]", "", data$plot)
data$plot <- as.numeric(data$plot)

data <- data[is.na(data$plot)==F,]

## CATEGORIZE FOR THE DATES
data[data$plot >= 1 & data$plot < 121 & data$Dlines==0,'category'] <- 1
data[data$plot >= 121 & data$plot < 241 & data$Dlines==0,'category'] <- 2
data[data$plot >= 1 & data$plot < 21 & data$Dlines==1,'category'] <- 3
data[data$plot >= 21 & data$plot < 41 & data$Dlines==1,'category'] <- 4

data[(data$plot >= 241 & data$plot < 361 & data$Dlines==0) |
       (data$plot > 41 & data$plot < 61 & data$Dlines==1),'category'] <- 5

data[(data$plot >= 361 & data$plot <= 480 & data$Dlines==0) |
       (data$plot >= 61 & data$plot <= 80 & data$Dlines==1),'category'] <- 6

data[data$plot==41 & data$Dlines==1,'category'] <- 7

###Setting up site numbers

data[(data$plot >= 1 & data$plot < 241 & data$Dlines==0) |
       (data$plot >= 1 & data$plot < 41 & data$Dlines==1),'site'] <- 1
data[(data$plot >= 241 & data$plot < 481 & data$Dlines==0) |
       (data$plot >= 41 & data$plot < 81 & data$Dlines==1),'site'] <- 2

data <- data[is.na(data$category)==F,]


## SET UP COL FOR PLANTING AND HARVEST DATES
data$pdate <- NA
data$hdate <- NA
data$DAP <- NA

#sample <- data[c(200:210),]
data[data$Dlines==1, 'plot'] <- paste('D',data[data$Dlines==1, 'plot'], sep='')

data$pdate <- apply(data,1,FUN=Pdate)
data$hdate <- apply(data,1,FUN=Hdate)

data$pdate <- as.Date(data$pdate, tryFormats = c("%m/%d/%Y"))
data$hdate <- as.character(data$hdate)
data$hdate <- as.Date(data$hdate, tryFormats = c("%m/%d/%Y"))
data$DAP <- data$hdate - data$pdate

data <- merge(data, genoCode, all.x=T)

write.csv(data,file = paste(wd,'out','DW 4 1 19.csv', sep='/'))


