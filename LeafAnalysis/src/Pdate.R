Pdate <- function(x){
  cat <- as.numeric(x['category'])
  plantDate <- dates[cat,'Planting.Date']
  return(plantDate)
}