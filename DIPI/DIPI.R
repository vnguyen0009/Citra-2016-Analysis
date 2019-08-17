###########################
## Author: Vy T. Nguyen  ##
## Last Edited: 2/8/2019 ##
###########################

# This code to calculate daily increase of plastochron index (DIPI)
# Plastochron index is a way of measuring plant age
# Instead of by days, it is by morphology

rm(list = ls())

library(plyr)
library(dplyr)
library(ggplot2)


wd <- getwd()
setwd(paste(wd, 'src', sep = '/'))

source('calcDIPI.R')

weather <- read.csv('FAWN_report (1).csv', header = T, stringsAsFactors = F, #weather data
                    na.strings = c('', '*', 'NA'))
setwd(paste(wd, 'out', sep = '/'))
###########################################
######## Now to convert dates to DAP ######
###########################################

weather1 <- weather
weather2 <- weather

s1Pdate <- as.Date('3/23/2016', tryFormats = c('%m/%d/%Y')) #planting dates for the seasons
s2Pdate <- as.Date('5/6/2016', tryFormats = c('%m/%d/%Y'))

# season 1 is the control and season 2 is the heat treatment

weather1$DAP <- as.numeric(as.Date(weather$Period, tryFormats = c('%m/%d/%Y')) - s1Pdate)
weather1$site <- 1

weather2$DAP <- as.numeric(as.Date(weather$Period, tryFormats = c('%m/%d/%Y')) - s2Pdate)
weather2$site <- 2

weatherAll <- rbind(weather1, weather2[weather2$DAP > 0 ,])
# piTemp is pulling out the weather data for the days plastochron index (PI) was taken
piTemp <- weatherAll[ifelse(weatherAll$site==1, 
                               weatherAll$DAP >=15 & weatherAll$DAP <=26,
                               weatherAll$DAP >=19 & weatherAll$DAP <=29),]


piTemp1 <- aggregate(piTemp[,c('Tmin', 'Tavg', 'Tmax','srad')], by=list(piTemp$site), FUN=mean, na.rm=T)
colnames(piTemp1)[1] <- c('site')

write.csv(piTemp1,file='weatherTable.csv')
# trying to see if there is a statistically significant difference in the weather of 
# growing seasons
t.test(piTemp$Tavg[piTemp$site==1], piTemp$Tavg[piTemp$site==2])
t.test(piTemp$Tmin[piTemp$site==1], piTemp$Tmin[piTemp$site==2])
t.test(piTemp$Tmax[piTemp$site==1], piTemp$Tmax[piTemp$site==2])
t.test(piTemp$srad[piTemp$site==1], piTemp$srad[piTemp$site==2])

#tail(weatherAll)

#####################################################################
######## set up data frame for the parameters Nm, Tbase, Topt1 ######
#####################################################################

para <- data.frame('geno' = c('CAL', 'RIJC217', 'RIJC223', 'JAM', 'RIJC016', 'RIJC208'), 
                   'Tbase' = c(7.4, 11.9, 10.3, 7.8, 8.4, 8.4),
                   'Topt1' = c(25.0, 19.8, 21.1, 25, 22.7, 24.0),
                   'Nm' = c(0.38, 0.34, 0.30, 0.38, 0.33, 0.39))


######################################################
####### Calculate DIPI using a UDF and by() ##########
######################################################

dp <- by(weatherAll, INDICES = list(weatherAll$Period, weatherAll$site), FUN = calcDIPI)
dpTab <- do.call("rbind", dp)

###################################################################
####### RESULTS means DIPI for each seson of select lines #########
###################################################################

setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'DIPI', 'out', sep = '/'))

#range(data$DAP, na.rm=T)
dpTab1 <- subset(dpTab, DAP >= 0 & DAP <=30)

avgDIPI <- aggregate(dpTab1$DIPI, by = list(dpTab1$site, dpTab1$geno), FUN=mean)
colnames(avgDIPI) <- c('site', 'geno', 'DIPI')
write.csv(avgDIPI, 'selectDIPI.csv')
