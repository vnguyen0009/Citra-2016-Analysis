rm(list = ls())

library(plyr)
## Author: Vy T. Nguyen
## Last Edited: 1/4/2019
library(dplyr)
library(ggplot2)


setwd(paste('C:/Users/Vy/Documents'))
wd <- getwd()
setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'DIPI', 'src', sep = '/'))

source('calcDIPI.R')

weather <- read.csv('FAWN_hourly temp.csv', header = T, stringsAsFactors = F, 
                    na.strings = c('', '*', 'NA'))
PI <- read.csv('pIndex_citra 2016.csv', header = T, stringsAsFactors = F, 
               na.strings = c('', '*', 'NA', "#REF!"))
geno <- read.csv('PtoGen.csv', header = T, stringsAsFactors = F, 
                 na.strings = c('', '*', 'NA'))
##head(weather)
##head(PI)
##head(geno)

###### CORRECTING PI DATA ######
### The code in this cell is from PI_citra 2016.pynb
## It is putting the PI file into the correct format that we want
PIT <- PI[c('Plot.ID', 'site', 'Plant.ID', 'PIT1', 'PIT2', 'PIT3', 'T1..DAP.', 'T2.DAP.', 'T3.DAP.')]
#head(PIT)
pit1 <- PI[c('Plot.ID', 'Plant.ID', 'PIT1', 'T1..DAP.', 'site')]
colnames(pit1) <- c('plot', 'rep', 'PI', 'DAP', 'site')

pit2 <- PI[c('Plot.ID', 'Plant.ID', 'PIT2', 'T2.DAP.', 'site')]
colnames(pit2) <- c('plot', 'rep', 'PI', 'DAP', 'site')

pit3 <- PI[c('Plot.ID', 'Plant.ID', 'PIT3', 'T3.DAP.', 'site')]
colnames(pit3) <- c('plot', 'rep', 'PI', 'DAP', 'site')

pit <- rbind(pit1, pit2)
pit <- rbind(pit, pit3)

noNApit <- pit[is.na(pit$PI)==F & is.na(pit$DAP)==F,]

data <- merge(noNApit, geno, all.x = T) 

#head(data)

######## Now to convert dates to DAP ######
weather1 <- weather
weather2 <- weather

s1Pdate <- as.Date('3/23/2016', tryFormats = c('%m/%d/%Y'))
s2Pdate <- as.Date('5/6/2016', tryFormats = c('%m/%d/%Y'))

weather1$DAP <- as.numeric(as.Date(weather$Period, tryFormats = c('%m/%d/%Y')) - s1Pdate)
weather1$site <- 1

weather2$DAP <- as.numeric(as.Date(weather$Period, tryFormats = c('%m/%d/%Y')) - s2Pdate)
weather2$site <- 2

weatherAll <- rbind(weather1, weather2[weather2$DAP > 0 ,])

piTemp <- weatherAll[weatherAll$DAP >=0 & weatherAll$DAP <=30,]
piTemp1 <- aggregate(weatherAll[,c('Tmin', 'Tavg', 'Tmax')], by = list(weatherAll$site), FUN=mean)
#tail(weatherAll)
###### Pull Lines and set parameters #######

selectLines <- data[grep('RIJC016|RIJC208|RIJC217|RIJC223|CAL|JAM',data$geno), ]

selectLines <- selectLines[order(selectLines$geno),]
head(selectLines)
# set up data frame for the parameters Nm, Tbase, Topt1
para <- data.frame('geno' = c('CAL', 'RIJC217', 'RIJC223', 'JAM', 'RIJC016', 'RIJC208'), 
                   'Tbase' = c(7.4, 11.9, 10.3, 7.8, 8.4, 8.4),
                   'Topt1' = c(25.0, 19.8, 21.1, 25, 22.7, 24.0),
                   'Nm' = c(0.38, 0.34, 0.30, 0.38, 0.33, 0.39))

####### Calculate DIPI using a UDF and by() ##########
dp <- by(weatherAll, INDICES = list(weatherAll$Period, weatherAll$site), FUN = calcDIPI)
dpTab <- do.call("rbind", dp)

####### RESULTS means DIPI for each seson of select lines #########
setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'DIPI', 'out', sep = '/'))

#range(data$DAP, na.rm=T)
dpTab1 <- subset(dpTab, DAP >= 15 & DAP <=30)

avgDIPI <- aggregate(dpTab1$DIPI, by = list(dpTab1$site, dpTab1$geno), FUN=mean)
colnames(avgDIPI) <- c('site', 'geno', 'DIPI')
write.csv(avgDIPI, 'selectDIPI.csv')
