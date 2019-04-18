Hdate <- function(x){
  hno <- as.numeric(x['HNO'])
  cat <- as.numeric(x['category'])
  if (hno == 11){
    plot <- x['plot']
    harvestDate <- FHdates[FHdates$plot==plot,'H11']
  }
  else{
    harvestDate <- dates[cat,paste('H', hno,sep='')]
  }
  return(as.character(harvestDate))
}