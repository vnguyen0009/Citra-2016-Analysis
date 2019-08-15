leafAge <- function(x){
  leafId <- x['variable']
  currentDay <- as.numeric(x['DAP'])
  firstDay <- as.numeric(min(meltLAGeno[meltLAGeno$variable==leafId,'DAP']), 
                         na.rm=T)
  return(currentDay - firstDay)
}