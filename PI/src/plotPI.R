plotPI <- function(x,y, data, fit){
  c <- ggplot(data)
  plot1 <- c + 
    #geom_smooth(aes(x=DAP, y=value), se=F, span=1) +
    geom_point(aes(x=x, y=y))+
    #facet_grid(.~plot) +
    coord_cartesian(ylim = c()) + 
    theme() +
    labs(x = 'Days After Planting (DAP)', y = paste('Plastechron Index (PI)'), 
         title = paste('Plot', data$plot,'Season', data$site, 'Plant', data$rep, sep=' ') ) + 
    theme_bw() + 
    scale_alpha(guide = 'none')  +
    geom_line(aes(x=x, y=predict(fit), color='red'))
  
  print(plot1)
}