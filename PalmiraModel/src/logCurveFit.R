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
    labs(title = paste(selectedLine), 
         x = deparse(substitute(x)), 
         y = deparse(substitute(y))
    ) +
    theme_bw() + 
    scale_alpha(guide = 'none')+
    annotate("text", label = paste("y=", signif(phi1, digits=4),
                                   "/(1+exp(-(", signif(phi2, digits=4),
                                   "+", signif(phi3, digits=4), "*x)))",
                                   "\n sum resid^2 = ", 
                                   signif(sum(resid(MSNODgc)^2), digits=3), 
                                   sep=" "), 
             x = max(x,na.rm=T)/3, y = max(y, na.rm=T)/3, size = 4, color="red")
  
  print(plotGC)
  
  return(summary(MSNODgc))
}