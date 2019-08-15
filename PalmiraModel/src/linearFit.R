linearFit <- function(x,y, dat){
  
  TOTALBRlm <- lm(y ~ x, data=dat, na.action = na.omit)
  summary(TOTALBRlm)
  
  b <- coef(TOTALBRlm)[1]
  m <- coef(TOTALBRlm)[2]
  xMSNOD <- c(min(x, na.rm=T):max(x, na.rm=T))
  yTOTALBR <- m*xMSNOD+b
  predictTOTALBR <- data.frame(xMSNOD,yTOTALBR)
  
  # Plot predicted values
  p <- ggplot(dat)
  
  plotLM <- p +
    geom_point(aes(x = x, 
                   y = y
    )) +
    geom_line(data=predictTOTALBR, aes(x=predictTOTALBR$xMSNOD, y=predictTOTALBR$yTOTALBR)) +
    coord_cartesian(ylim = c()) + 
    theme(legend.position = "bottom")+
    labs(title = paste(selectedLine), 
         x = deparse(substitute(x)), 
         y = deparse(substitute(y))
    ) +
    theme_bw() + 
    scale_alpha(guide = 'none') +
    annotate("text", label = paste("y= ", signif(m, digits=4), "x + ", 
                                   signif(b, digits=4),"\n r^2 = ", 
                                   signif(summary(TOTALBRlm)$r.squared, digits=3), 
                                   sep=" "), 
             x = max(x,na.rm=T)/3, y = max(y, na.rm=T)/3, size = 4, color="red")
  
  print(plotLM)
  return(summary(TOTALBRlm))
  
}