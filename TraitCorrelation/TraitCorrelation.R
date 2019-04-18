rm(list=ls())

library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)

wd <- getwd()

setwd(paste(wd,'data', sep='/'))

DW <- read.csv('DW 4 1 19.csv', row.names=1, na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)
seed <- read.csv('FH seed 9 11 18.csv', row.names=1, na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)
flapjack <- read.csv('flapjack.csv', na.strings = c('', 'NA', '*', '.', '?'), stringsAsFactors = F)

FH <- DW[DW$HNO==11,]
seed100 <- seed[,c('plot', 'X100SD')]
GH <- flapjack[,c('RIL', 'Fin')]

data <- merge(seed100, FH, all.x=T, by='plot')
data <- data[,c(1,2,4:12,15,19)]
meanData <- aggregate(data[,c(2:11)], by = list(data$geno, data$site), na.rm=T, FUN=mean)

colnames(meanData)[1:2] <- c('geno', 'site')

meltedData <- melt(meanData, id = c('site', 'geno', 'X100SD'))
wGH <- merge(meltedData, GH, by.x = 'geno', by.y='RIL', all.x=T)
RIL <- wGH[is.na(wGH$Fin)==F,]

p <- ggplot(meltedData)

allPlot <- p + #geom_smooth(aes(x = meltedData$X100SD, 
             #               y = meltedData$value, 
              #              color=as.factor(meltedData$site)), se=F) + 
  geom_point(aes(x = meltedData$X100SD, 
                 y = meltedData$value, 
                 #color=as.factor(meltedData$site)), 
                alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Population', 
       x = '100 Seed Weight (g)', 
       y = 'Trait Value (g)', 
       color = 'Growing Seasons' ) +
  #scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meltedData$variable, scales='free') +
  theme_bw() + 
  scale_alpha(guide = 'none')

seasonPlot <- p + #geom_smooth(aes(x = meltedData$X100SD, 
  #               y = meltedData$value, 
  #              color=as.factor(meltedData$site)), se=F) + 
  geom_point(aes(x = meltedData$X100SD, 
                 y = meltedData$value, 
                 color=as.factor(meltedData$site), 
                 alpha=0.5)) +
  coord_cartesian(ylim = c()) + 
  theme(legend.position = "bottom")+
  labs(title = 'Total Population', 
       x = '100 Seed Weight (g)', 
       y = 'Trait Value (g)', 
       color = 'Season' ) +
  scale_fill_discrete(name = 'Growing Seasons') + 
  facet_wrap(~meltedData$variable, scales='free') +
  theme_bw() + 
  scale_alpha(guide = 'none')

r <- ggplot(RIL)
determPlot <- r + #geom_smooth(aes(x = meltedData$X100SD, 
    #               y = meltedData$value, 
    #              color=as.factor(meltedData$site)), se=F) + 
    geom_point(aes(x = RIL$X100SD, 
                   y = RIL$value, 
                   color=as.factor(RIL$Fin), 
                   alpha=0.5)) +
    coord_cartesian(ylim = c()) + 
    theme(legend.position = "bottom")+
    labs(title = 'RIL Population', 
         x = '100 Seed Weight (g)', 
         y = 'Trait Value (g)', 
         color = 'Fin' ) +
    #scale_fill_discrete(name = 'Growing Seasons') + 
    facet_wrap(~RIL$variable, scales='free') +
    theme_bw() + 
    scale_alpha(guide = 'none')

pdf(paste(wd, 'out', 'TraitVs100SeedWt.pdf', sep='/'), onefile = T)
print(allPlot)
print(seasonPlot)
print(determPlot)
dev.off()
