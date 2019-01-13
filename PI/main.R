rm(list = ls())

library(plyr)
library(dplyr)
library(ggplot2)

setwd(paste('C:/Users/Vy/Documents'))
wd <- getwd()
setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'PI', 'src', sep = '/'))

source('findSlope.R')
source('plotPI.R')

######LOAD DATA ##########
data <- read.csv('pIndex_citra 2016.csv', header = T, na.strings = c('NA', '', '#NUM!', '#REF!', '#VALUE!'), stringsAsFactors=T)
geno <- read.csv('PtoGen.csv', header = T, stringsAsFactors=T)
flapjack <- read.csv('flapjack.csv', header = T, stringsAsFactors = T) 
#Flapjack was uploaded to find the determinency of the lines


####### CORRECT PI FILE FORMAT #######
PIT <- data[c('Plot.ID', 'site', 'Plant.ID', 'PIT1', 'PIT2', 'PIT3', 'T1..DAP.', 'T2.DAP.', 'T3.DAP.')]
pit1 <- data[c('Plot.ID', 'Plant.ID', 'PIT1', 'T1..DAP.', 'site')]
colnames(pit1) <- c('plot', 'rep', 'PI', 'DAP', 'site')

pit2 <- data[c('Plot.ID', 'Plant.ID', 'PIT2', 'T2.DAP.', 'site')]
colnames(pit2) <- c('plot', 'rep', 'PI', 'DAP', 'site')

pit3 <- data[c('Plot.ID', 'Plant.ID', 'PIT3', 'T3.DAP.', 'site')]
colnames(pit3) <- c('plot', 'rep', 'PI', 'DAP', 'site')

pit <- rbind(pit1, pit2)
pit <- rbind(pit, pit3)

noNApit <- pit[is.na(pit$PI)==F & is.na(pit$DAP)==F,]

###### FIND SLOPES PI VS DAP ########
byPlot = by(noNApit, INDICES=list(noNApit$site, noNApit$plot, noNApit$rep), FUN=findSlope)

results <- do.call("rbind", byPlot)
resultsRIL <- merge(results, geno, by = 'plot', all.x=T)

meanRIL <- aggregate(resultsRIL[,c('slope', 'inter', 'r2')], by = list(resultsRIL$geno, resultsRIL$site), FUN=mean, na.action=na.rm)
colnames(meanRIL)[1:2] <- c('geno', 'site')

seasonResults <- aggregate(results[c('slope', 'inter', 'r2')], by = list(results$site), FUN =mean, na.rm=T)
colnames(seasonResults)[1] = 'site'
s1Slope <- meanRIL[meanRIL$site==1, 'slope']
s2Slope <- meanRIL[meanRIL$site==2, 'slope']

t.test(s1Slope,s2Slope)
p <- ggplot(meanRIL, aes(label = geno))

############ ALL POPULATION PLOTS ###############################
setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'PI', 'out', sep = '/'))
pdf('results.pdf', onefile = T)

plot1 <- p + geom_density(aes(x=meanRIL$slope)) +
  coord_cartesian(ylim = c()) + 
  theme() +
  labs( x = paste('Slope'), 
        title = paste('Histogram of Linear Regression Slopes for Citra 2016', sep=' ') ) + 
  theme_bw() + 
  scale_alpha(guide = 'none')+
  facet_grid(rows = vars(meanRIL$site)) +
  geom_text(aes(label = ifelse(meanRIL$geno %in% meanRIL[grep('CAL|JAM', meanRIL$geno),'geno'],
                               as.character(geno), ''), x = meanRIL$slope, y = 5))

print(plot1)
#################### RIJC PLOTS ###################################

RIJC <- meanRIL[grep('RIJC|CAL|JAM',meanRIL$geno), ]
Fin <- flapjack[, c('RIL', 'Fin')] #set up variable to find determinency

determ <- merge(RIJC, Fin, by.x = 'geno', by.y = 'RIL', all.x=T)

#unique(determ[is.na(determ$Fin),'geno'])

S1RIJC <- RIJC[RIJC$site==1, 'slope']
S2RIJC <- RIJC[RIJC$site==2, 'slope']

t.test(S1RIJC,S2RIJC)

q <- ggplot(RIJC, aes(label = geno))

plot2 <- q + geom_density(aes(x=RIJC$slope)) +
  coord_cartesian(ylim = c()) + 
  theme() +
  labs( x = paste('Slope'), 
        title = paste('RIJC Histogram of Linear Regression Slopes for Citra 2016', sep=' ') ) + 
  theme_bw() + 
  scale_alpha(guide = 'none')+
  facet_grid(rows = vars(RIJC$site)) +
  geom_text(aes(label = ifelse(RIJC$geno %in% RIJC[grep('CAL|JAM', RIJC$geno),'geno'],
                               as.character(geno), ''), x = RIJC$slope, y = 5))

print(plot2)

r <- ggplot(determ[is.na(determ$Fin)==F,], aes())

plot3 <- r + geom_density(aes(x=slope, fill=as.factor(Fin), alpha=0.2)) +
  coord_cartesian(ylim = c()) + 
  theme() +
  labs( x = paste('Slope'), 
        title = paste('PI Citra 2016', sep=' ') ) + 
  theme_bw() + 
  scale_alpha(guide = 'none')+
  facet_grid(rows = vars(site)) + 
  scale_fill_discrete(name="Determinacy",
        labels=c("Indet", "Det"))

print(plot3)
dev.off()

############## SELECT LINES #####################################
selectLines <- meanRIL[grep('RIJC016|RIJC208|RIJC217|RIJC223|CAL|JAM',meanRIL$geno), ]

#selectLines <- selectLines[order(selectLines$geno),]
#write.csv(selectLines, 'selectLines.csv')
