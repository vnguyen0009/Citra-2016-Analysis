
# This leaf analysis is to fit leaf data to linear or non-linear models for the FSPM
# In this code, I am trying to see if once the data is normalized to leaf age rather than
# days after planting, the leaves follow a similar growth curve

rm(list=ls())
#install.packages('minpack.lm')
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stats)
library(minpack.lm)
library(ggbiplot)

wd <- getwd()
source(paste(wd,'src', 'Hdate.R', sep='/'))
source(paste(wd,'src', 'Pdate.R', sep='/'))
source(paste(wd, 'src', sep='/', 'logCurveFit.R'))

setwd(paste(wd,'data', sep='/'))

data <- read.csv('parentLeafArea.csv', 
                na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)
dates <- read.csv('dates.csv',stringsAsFactors = F)
genoCode <- read.csv('plotToGeno.csv', stringsAsFactors = F)
FHdates <- read.csv('FH dates.csv', stringsAsFactors = F)


data <- data[,c(1:5)]

data$HNO <- gsub("H", "", data$HNO)
data$HNO <- gsub("F", 11, data$HNO)

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
data$DAP <- as.numeric(data$DAP)

line <- data[data$geno=='CAL'& data$site==2,]

leafAge <- function(x){ #This is to calculate the age of a secific leaf based on when it appeared
  plotId <- x['plot'] 
  leafId <- x['Code']
  currentDay <- as.numeric(x['DAP'])
  firstDay <- as.numeric(min(data[data$plot==plotId & data$Code==leafId,'DAP']), 
                         na.rm=T)
  return(currentDay - firstDay)
}

line$LeafAge <- apply(line,1, FUN=leafAge)

MSleaf <- line[grep('V', line$Code),] 
#Only using the first five node because the others seem like outliers


p <- ggplot(MSleaf)

plot1 <- p + geom_smooth(aes(x = MSleaf$LeafAge, 
                             y = MSleaf$cm.2, 
                             color=as.factor(MSleaf$Code)), se=F) + 
  geom_point(aes(x= MSleaf$LeafAge, 
                 y = MSleaf$cm.2, 
                 color=as.factor(MSleaf$Code), alpha=0.5)) #+
  #coord_cartesian(ylim = c()) + 
  #theme(legend.position = "bottom")+
  #labs(title = 'Leaf Area', 
   #    x = 'Days after Planting (DAP)', 
  #     y = 'Leaf Area (cm^2)', 
   #    color = 'Leaf Position' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  #facet_wrap(~parentMeans$geno) +
  #theme_bw() + 
  #scale_alpha(guide = 'none')

print(plot1)

plot(MSleaf$DAP, MSleaf$cm.2)

logCurveFit(x=MSleaf$LeafAge,y=MSleaf$cm.2,
            dat = MSleaf,phi1Start = 100, phi2Start=-3.39,phi3Start = 0.079)


dev.off()
