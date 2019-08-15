rm(list = ls())

#install.packages("C:/Users/Vy/Downloads/asreml_4.1.0.98.zip", repos = NULL, type ="win.binary")
#install.packages('data.table')
#install.packages('jsonlite')

library(plyr)
library(dplyr)
library(ggplot2)
library(asreml) ##### ASREML LICENSE IS REQUIRED FOR THIS CODE
library(nadiv)
library(asremlPlus)
library(data.table)
library(jsonlite)

setwd(paste('C:/Users/Vy/Documents'))
wd <- getwd()
setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'NormPI', 'src', sep = '/'))

REMLfile<-read.csv("CitraASReml.csv",h=T,na.string='*')
data <- read.csv('pIndex_citra 2016_1-14-19.csv', header=T, na.string = c('NA', '', '#REF!', '#NUM!', '#VALUE!'))

setwd(paste(wd,'work','Data','Citra2016','Citra-2016-Analysis', 'NormPI', 'out', sep = '/'))

PIT <- data[c('Plot.ID', 'site', 'Plant.ID', 'PIT1', 'PIT2', 'PIT3', 'T1..DAP.', 'T2.DAP.', 'T3.DAP.')]
colnames(PIT) <- c('plot', 'site', 'rep', 'PIT1', 'PIT2', 'PIT3', 'T1DAP', 'T2DAP', 'T3DAP')
plotPIT <- aggregate(PIT[,c('PIT1', 'PIT2', 'PIT3', 'T1DAP', 'T2DAP', 'T3DAP')], by=list(PIT$plot), FUN=mean, na.rm=T)
colnames(plotPIT)[1] <- c('plot')

spatial <- merge(REMLfile, plotPIT, by.x=c('IDSORT'), by.y = c('plot'), all.x=T)
y = 'PIT1' #Trait name

spatial$y<-spatial$PIT1 ###Change this to the trait col you are looking for

head(spatial) ###Use head() and tail() to check data
tail(spatial)

spatial$X<-spatial$Row    # X coordinate as a factor
spatial$Y<-spatial$Col    # Y coordinate as a factor
spatial$Row<-as.factor(spatial$Row)
spatial$Col<-as.factor(spatial$Col)
spatial$Rep<-as.factor(spatial$Rep)
spatial$Iblock<-as.factor(spatial$Iblock)
spatial$Genot<-as.factor(spatial$Genot)
head(spatial)
str(spatial)

##### Basic Model without spatial components but design components ######
nospatial<-asreml(fixed=y~1,
                  random=~Rep:Iblock+Genot,
                  data=spatial)
summary(nospatial)$varcomp

plot(nospatial)

(h2b<-nadiv:::pin(nospatial,h2b~V2/(V1+V2+V3)))

###### incorporating Spatial Autocorrelation on X and Y #########
spatial<-spatial[order(spatial$Col,spatial$Row),] ###Necessary even if you don't need variogram
spatial1<-asreml(fixed=y~1,
                 random=~Rep:Iblock+Genot, #rep by block interaction
                 rcov=~ar1(Col):ar1(Row),na.method.Y='include', 
                 data=spatial) #testing whether or not rows and col have effect
summary(spatial1)$varcomp
plot(spatial1)
(h2b<-nadiv:::pin(spatial1,h2b~V2/(V1+V2+V3))) #narrow sense heritability 

####### Using asremlplus ######
reml.lrt.asreml(spatial1,nospatial,positive.zero=FALSE)  # LRT with pvalues
info.crit.asreml(nospatial)                           # AIC and BIC
info.crit.asreml(spatial1)                           # AIC and BIC

###### Detecting Outliers #####
res<-spatial1$residuals
#Look at residuals plot first before you ID the outliers and delete them
ID<-which.max(res)
#ID<-which.min(res)
out <- spatial[ID,]
#spatial$y[spatial$IDSORT==out$IDSORT]<-NA

###### Some predictions with Genotypes as Fixed Effect ######
spatfix<-asreml(fixed=y~Genot,
                random=~Rep:Iblock,
                rcov=~ar1(Col):ar1(Row),na.method.Y='include',
                data=spatial)
plot(spatfix)

pred.spatfix<-predict(spatfix,classify="Genot",sed=TRUE)$predictions$pvals
head(pred.spatfix)
summary(spatfix)$varcomp

write.csv(pred.spatfix,file=paste('normalized',y,'_',m,'.csv'))
