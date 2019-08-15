logCurveFit <- function(x,y,dat, phi1Start, phi2Start, phi3Start){
  library(minpack.lm)
  MSNODgc <- nlsLM(y ~ phi1/(1+exp(-(phi2+phi3*x))), data=dat,
                 start=list(phi1=phi1Start,phi2=phi2Start,phi3=phi3Start),trace=TRUE)
  summary(MSNODgc)
  
  # Create Predicted Values
  phi1 <- coef(MSNODgc)[1]
  phi2 <- coef(MSNODgc)[2]
  phi3 <- coef(MSNODgc)[3]
  
  xIn <- c(min(x, na.rm=T):max(x, na.rm=T))
  yMSNOD <- phi1/(1+exp(-(phi2+phi3*xIn)))
  predictMSNOD <- data.frame(xIn,yMSNOD)
  
  p <- ggplot(dat)
  
  plotGC <- p + 
    geom_point(aes(x = x, 
                   y = y
                   
    )) +
    geom_line(data=predictMSNOD, aes(x=predictMSNOD$x, y=predictMSNOD$yMSNOD)) +
    coord_cartesian(ylim = c()) + 
    theme(legend.position = "bottom")+
    labs( 
         x = paste('x'), 
         y = paste('y')
    ) +
    theme_bw() + 
    scale_alpha(guide = 'none')
  
  print(plotGC)
  
  return(summary(MSNODgc))
}