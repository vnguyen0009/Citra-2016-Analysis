## This code is using plastochron index, fitting it to a linear model
## The slope is our DIPI or daily increase of plastochron index

rm(list = ls())

library(plyr)
library(dplyr)
library(ggplot2)

#setwd(paste('C:/Users/Vy/Documents'))
wd <- getwd()
setwd(paste(wd, 'src', sep = '/'))

source('findSlope.R')
source('plotPI.R')

##########################
######LOAD DATA ##########
##########################

PI <- read.csv('pIndex_citra 2016_1-14-19.csv', header = T, na.strings = c('NA', '', '#NUM!', '#REF!', '#VALUE!'), stringsAsFactors=T)
geno <- read.csv('PtoGen.csv', header = T, stringsAsFactors=T)
flapjack <- read.csv('flapjack.csv', header = T, stringsAsFactors = T) 
#Flapjack was uploaded to find the determinency of the lines

setwd(paste(wd, 'out', sep = '/'))
pdf('results.pdf', onefile = T)
######################################
####### CORRECT PI FILE FORMAT #######
######################################

PIT <- PI[c('Plot.ID', 'site', 'Plant.ID', 'PIT1', 'PIT2', 'PIT3', 'T1..DAP.', 'T2.DAP.', 'T3.DAP.')]
#head(PIT)
pit1 <- PI[c('Plot.ID', 'Plant.ID', 'PIT1', 'T1..DAP.', 'site')]
colnames(pit1) <- c('plot', 'rep', 'PI', 'DAP', 'site')
pit1Geno <- merge(pit1, geno, all.x=T)
pit1Mean <- aggregate(pit1Geno[,c('DAP','PI')], by = list(pit1Geno$site, pit1Geno$geno), FUN=mean, na.rm=T)
colnames(pit1Mean)[1:2] <- c('site', 'geno')

pit2 <- PI[c('Plot.ID', 'Plant.ID', 'PIT2', 'T2.DAP.', 'site')]
colnames(pit2) <- c('plot', 'rep', 'PI', 'DAP', 'site')
pit2Geno <- merge(pit2, geno, all.x=T)
pit2Mean <- aggregate(pit2Geno[,c('DAP','PI')], by = list(pit2Geno$site, pit2Geno$geno), FUN=mean, na.rm=T)
colnames(pit2Mean)[1:2] <- c('site', 'geno')

pit3 <- PI[c('Plot.ID', 'Plant.ID', 'PIT3', 'T3.DAP.', 'site')]
colnames(pit3) <- c('plot', 'rep', 'PI', 'DAP', 'site')
pit3Geno <- merge(pit3, geno, all.x=T)
pit3Mean <- aggregate(pit3Geno[,c('DAP','PI')], by = list(pit3Geno$site, pit3Geno$geno), FUN=mean, na.rm=T)
colnames(pit3Mean)[1:2] <- c('site', 'geno')

pit <- rbind(pit1Mean, pit2Mean)
pit <- rbind(pit, pit3Mean)

noNApit <- pit[is.na(pit$PI)==F & is.na(pit$DAP)==F,]
#write.csv(pit3Geno, file='RAW_PI.csv')


#####################################
###### FIND SLOPES PI VS DAP ########
#####################################
byPlot = by(noNApit, INDICES=list(noNApit$site, noNApit$geno), FUN=findSlope)


resultsRIL <- do.call("rbind", byPlot)
#parents <- resultsRIL[grep('CAL|JAM',resultsRIL$geno), ]
#resultsRIL <- merge(results, geno, by = 'plot', all.x=T)

write.csv(resultsRIL, file = 'DIPI.csv')

meanRIL <- aggregate(resultsRIL[,c('slope', 'inter', 'r2')], by = list(resultsRIL$geno, resultsRIL$site), FUN=mean, na.rm=T)
colnames(meanRIL)[1:2] <- c('geno', 'site')

RIJC <- meanRIL[grep('RIJC|CAL|JAM',meanRIL$geno), ]
write.csv(RIJC, file='DIPI_RIJC.csv')
Fin <- flapjack[, c('RIL', 'Fin')] #set up variable to find determinency

determ <- merge(RIJC, Fin, by.x = 'geno', by.y = 'RIL', all.x=T)

unique(determ[is.na(determ$Fin),'geno'])

S1RIJC <- RIJC[RIJC$site==1, c('geno','slope')]
S2RIJC <- RIJC[RIJC$site==2, c('geno','slope')]

t.test(S1RIJC$slope,S2RIJC$slope)

mergedRIJC <- merge(S1RIJC, S2RIJC, by='geno', all=T)
mergedRIJC$difference <- mergedRIJC$slope.x - mergedRIJC$slope.y

hist(mergedRIJC$difference, 
     main=paste('differences between season genotype slopes'),
     xlab=paste('DIPI Difference'))

mean(mergedRIJC$difference, na.rm=T)


############ ALL POPULATION PLOTS ###############################
p <- ggplot(meanRIL, aes(label = geno))

plot1 <- p + geom_density(aes(x=meanRIL$slope)) +
  coord_cartesian(ylim = c()) + 
  theme(text = element_text(size=20)) +
  labs( x = paste('Slope'), 
        title = paste('All Genotypes Daily Increase of Plastechron', sep=' ') ) + 
  theme_bw() + 
  scale_alpha(guide = 'none')+
  facet_grid(rows = vars(meanRIL$site)) +
  geom_text(aes(label = ifelse(meanRIL$geno %in% meanRIL[grep('CAL|JAM', meanRIL$geno),'geno'],
                               as.character(geno), ''), x = meanRIL$slope, y = 5))

print(plot1)

#################### RIJC PLOTS ###################################
q <- ggplot(RIJC, aes(label = geno))

plot2 <- q + geom_density(aes(x=RIJC$slope)) +
  coord_cartesian(ylim = c()) + 
  theme(text = element_text(size=20)) +
  labs( x = paste('DIPI'), 
        title = paste('RIJC Daily Increase of Plastechron Index', sep=' ') ) + 
  theme_bw() + 
  scale_alpha(guide = 'none')+
  facet_grid(rows = vars(RIJC$site)) +
  geom_text(aes(label = ifelse(RIJC$geno %in% RIJC[grep('CAL|JAM', RIJC$geno),'geno'],
                               as.character(geno), ''), x = RIJC$slope, y = 5))

print(plot2)

r <- ggplot(determ[is.na(determ$Fin)==F,], aes())

plot3 <- r + geom_density(aes(x=slope, fill=as.factor(Fin), alpha=0.2)) +
  coord_cartesian(ylim = c()) + 
  theme(text = element_text(size=20)) +
  labs( x = paste('DIPI'), 
        title = paste('RIJC Daily Increase of Plastechron Index', sep=' ') ) + 
  theme_bw() + 
  scale_alpha(guide = 'none')+
  facet_grid(rows = vars(site)) + 
  scale_fill_discrete(name="Determinacy",
        labels=c("Indet", "Det"))

print(plot3)

plot4 <- ggplot() +
  geom_point(data = noNApit[grep('RIJC',noNApit$geno), ], aes(DAP,PI), alpha = 0.1) + 
  geom_smooth(data = noNApit[grep('CAL|JAM',noNApit$geno), ], aes(x=DAP, y=PI, col=geno), 
              method='lm', se=F) +
  labs(x = 'Days After Planting (DAP)', y = paste('Plastechron Index (PI)'), 
       title = paste('Plastechron Index of RIL population over time')) +
  facet_grid(site~.) +
  theme(text = element_text(size=20)) +
  theme_bw()

print(plot4)

dev.off()

############## SELECT LINES #####################################
selectLines <- meanRIL[grep('RIJC016|RIJC208|RIJC217|RIJC223|CAL|JAM',meanRIL$geno), ]

selectLines <- selectLines[order(selectLines$geno),]
write.csv(selectLines, 'selectLines.csv')
