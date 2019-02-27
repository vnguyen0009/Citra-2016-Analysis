findSlope <- function(x){
  #print(x)
  DAP  <- as.numeric(x[,'DAP'])
  PI <- as.numeric(x[,'PI'])
  
  fit <- lm(PI ~ DAP, data=x, na.action=na.omit)
  
  result <-  data.frame('slope' = summary(fit)$coefficient[2], 
                        'inter' = summary(fit)$coefficient[1],
                        'sd' = sd(x$PI),
                        'geno' = x[1,'geno'], 
                        'site' = x[1,'site'],
                        'r2' = summary(fit)$r.squared)
  #plotPI(DAP,PI,x,fit)
  
  return(result)
}