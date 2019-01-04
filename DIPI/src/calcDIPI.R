calcDIPI <- function(x){
  ## x is the subset of an hourly weather data frame
  nGeno <- length(unique(para$geno))
  
  results <- data.frame('geno' = unique(para$geno), ##need to change row number for data frame
                        'site' = numeric(nGeno), 
                        'DAP' = numeric(nGeno), 
                        'DIPI' = numeric(nGeno),
                        'Tavg' = numeric(nGeno),
                        stringsAsFactors=F)
  
  for (jj in unique(para$geno)){
    ## To test function
    ##x <- weatherAll[weatherAll$DAP==30 & weatherAll$site==2,]
    ##jj = 'JAM'
    ##
    Tbase <- as.numeric(para[para$geno==jj,'Tbase'])
    Topt1 <- as.numeric(para[para$geno==jj,'Topt1'])
    Nm <- as.numeric(para[para$geno==jj,'Nm'])
    #geno <- para[para$geno==jj,'geno']
    site <- x[,'site']
    DAP <- x[,'DAP']
    Tday <- x[,'Tavg']
    
    ### You need to use ifelse() for arguments of length > 1
    f <- ifelse(Tday < Tbase, 0, 
                ifelse(Topt1 <= Tday, 1/24,(Tday-Tbase)/(24*(Topt1-Tbase))))
    
    factor = sum(f)
    
    DIP <- Nm*factor
    
    #results$geno[jj] <- jj
    results$site[results$geno== jj] <- site[1]
    results$DAP[results$geno== jj] <- DAP[1]
    results$DIPI[results$geno== jj] <- DIP
    results$Tavg <- mean(Tday)
    
  }
  #results$Tavg <- Tday
  return(results)
}