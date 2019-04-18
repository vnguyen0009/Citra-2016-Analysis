rm(list=ls())

library(ggplot2)
library(dplyr)

wd <- getwd()

setwd(paste(wd,'data', sep='/'))

DW <- read.csv('DW 4 1 19.csv', na.strings = c('', 'NA', '*', '.', '?'), 
               stringsAsFactors = F, row.names = 1)
flapjack <- read.csv('flapjack.csv', na.strings = c('', 'NA', '*', '.', '?'), 
                     stringsAsFactors = F)

Fin <- flapjack[,c('RIL', 'Fin')]
Fin[Fin$Fin==1,'GH'] <- 'Indeterminant'
Fin[Fin$Fin==2,'GH'] <- 'Determinant'

DW <- DW[DW$WHOLE < max(DW$WHOLE),]

meanDW <- aggregate(DW[,c('DWTOTALLEAF', 'TOTALSTEM', 'WHOLE', 'DAP')], 
                    by=list(DW$site, DW$HNO, DW$geno), FUN=mean, na.rm=T)
colnames(meanDW)[1:3] <- c('site', 'HNO', 'geno')
meanDW <- merge(meanDW, Fin[,c('GH', 'RIL')], by.x='geno', by.y='RIL', all.x=T)

lines <- c('CAL', 'JAM')

parentMeans <- meanDW[meanDW$geno=='CAL'|meanDW$geno=='JAM',]

dLines <- DW[DW$Dlines==1,]
meanDLines <- aggregate(dLines[,c('DWTOTALLEAF', 'TOTALSTEM', 'WHOLE', 'DAP')], 
                        by=list(dLines$site, dLines$HNO, dLines$geno), FUN=mean, na.rm=T)
colnames(meanDLines)[1:3] <- c('site', 'HNO', 'geno')
################### All geno #####################

p <- ggplot(meanDW)

leafPlot <- p +
  geom_smooth(aes(x = meanDW$DAP, 
                  y = meanDW$DWTOTALLEAF, color=as.factor(meanDW$site)), se=F) + 
  geom_point(aes(x=meanDW$DAP, 
                 y = meanDW$DWTOTALLEAF, color=as.factor(meanDW$site))) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom") +
  labs(title = 'Total Leaf Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Leaf Dry Weight (g)', 
       color = 'Growing Seasons',
       shape = 'Growth Habit') +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meanDW$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

stemPlot <- p +
  geom_smooth(aes(x = meanDW$DAP, 
                  y = meanDW$TOTALSTEM, color=as.factor(meanDW$site)), se=F) + 
  geom_point(aes(x=meanDW$DAP, 
                 y = meanDW$TOTALSTEM, color=as.factor(meanDW$site), alpha = 0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom") +
  labs(title = 'Total Stem Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Stem Dry Weight (g)', 
       color = 'Growing Seasons',
       shape = 'Growth Habit') +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meanDW$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

wholePlot <- p +
  geom_smooth(aes(x = meanDW$DAP, 
                  y = meanDW$WHOLE, color=as.factor(meanDW$site)), se=F) + 
  geom_point(aes(x=meanDW$DAP, 
                 y = meanDW$WHOLE, color=as.factor(meanDW$site), alpha = 0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom") +
  labs(title = 'Total Plant Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Plant Dry Weight (g)', 
       color = 'Growing Seasons',
       shape = 'Growth Habit') +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meanDW$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

pdf(paste(wd, 'out', 'allGenoPlot.pdf', sep='/'), width=40, height=40, onefile = T)
print(leafPlot)
print(stemPlot)
print(wholePlot)
dev.off()

#################Parent Plots#######################################
q <- ggplot(parentMeans)

leafParentPlot <- q + geom_smooth(aes(x = parentMeans$DAP, 
                                  y = parentMeans$DWTOTALLEAF, 
                                  color=as.factor(parentMeans$site)), se=F) + 
  geom_point(aes(x=parentMeans$DAP, y = parentMeans$DWTOTALLEAF, 
                 color=as.factor(parentMeans$site), alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Leaf Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Leaf Dry Weight (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~parentMeans$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

stemParentPlot <- q + geom_smooth(aes(x = parentMeans$DAP, 
                                      y = parentMeans$TOTALSTEM, 
                                      color=as.factor(parentMeans$site)), se=F) + 
  geom_point(aes(x=parentMeans$DAP, y = parentMeans$TOTALSTEM, 
                 color=as.factor(parentMeans$site), alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Stem Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Stem Dry Weight (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~parentMeans$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

wholeParentPlot <- q + geom_smooth(aes(x = parentMeans$DAP, 
                                      y = parentMeans$WHOLE, 
                                      color=as.factor(parentMeans$site)), se=F) + 
  geom_point(aes(x=parentMeans$DAP, y = parentMeans$WHOLE, 
                 color=as.factor(parentMeans$site), alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Plant Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Plant Dry Weight (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~parentMeans$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

pdf(paste(wd, 'out', 'parentPlots.pdf', sep='/'), onefile = T)
print(leafParentPlot)
print(stemParentPlot)
print(wholeParentPlot)
dev.off()

##################D-lines ###########################
r <- ggplot(meanDLines)

leafDPlot <- r + geom_smooth(aes(x = meanDLines$DAP, 
                                      y = meanDLines$DWTOTALLEAF, 
                                      color=as.factor(meanDLines$site)), se=F) + 
  geom_point(aes(x=meanDLines$DAP, y = meanDLines$DWTOTALLEAF, 
                 color=as.factor(meanDLines$site), alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Leaf Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Leaf Dry Weight (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meanDLines$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

stemDPlot <- r + geom_smooth(aes(x = meanDLines$DAP, 
                                      y = meanDLines$TOTALSTEM, 
                                      color=as.factor(meanDLines$site)), se=F) + 
  geom_point(aes(x=meanDLines$DAP, y = meanDLines$TOTALSTEM, 
                 color=as.factor(meanDLines$site), alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Stem Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Stem Dry Weight (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meanDLines$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

wholeDPlot <- r + geom_smooth(aes(x = meanDLines$DAP, 
                                       y = meanDLines$WHOLE, 
                                       color=as.factor(meanDLines$site)), se=F) + 
  geom_point(aes(x=meanDLines$DAP, y = meanDLines$WHOLE, 
                 color=as.factor(meanDLines$site), alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Plant Dry Weight', 
       x = 'Days after Planting (DAP)', 
       y = 'Total Plant Dry Weight (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meanDLines$geno) +
  theme_bw() + 
  scale_alpha(guide = 'none')

pdf(paste(wd, 'out', 'dPlots.pdf', sep='/'), onefile = T)
print(leafDPlot)
print(stemDPlot)
print(wholeDPlot)
dev.off()


