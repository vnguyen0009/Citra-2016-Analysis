threshold=function(interval_mapping_obj, LOD_thr, LOD_range=2){
  # compute differential of lod values (directions of change in lod per cM movement down the genome)
  interval_mapping_obj[2:(nrow(interval_mapping_obj)),'lod_delta']=interval_mapping_obj[-1,'lod']-interval_mapping_obj[-nrow(interval_mapping_obj),'lod']
  
  # markers with lod < thr OR negative movement of lod values are to be removed (rmv=T)
  interval_mapping_obj[interval_mapping_obj$lod<LOD_thr |
                       interval_mapping_obj$lod_delta<0,'rmv']=T
  
  # markers with lod > thr are, for now, labeled as not to be removed
  interval_mapping_obj[interval_mapping_obj$lod>=LOD_thr,'rmv']=F
  
  # for each marker in the genome
  for(i in 1:nrow(interval_mapping_obj)){
    
    # If currently labled as not to be removed
    if(!interval_mapping_obj[i,'rmv']){
      
      # If both the current and next marker have positive los movement
      if(interval_mapping_obj[i,'lod_delta']>0 & interval_mapping_obj[i+1,'lod_delta']>0){
        # label current marker for removal
        interval_mapping_obj[i,'rmv']=T
      
      # If current marker has positive lod movement, but next marker has negative lod movement
      }else if(interval_mapping_obj[i,'lod_delta']>0 & interval_mapping_obj[i+1,'lod_delta']<0){
        # label next marker for removal
        interval_mapping_obj[i+1,'rmv']=T
      
      # If current marker has negative lod movement
      }else if(interval_mapping_obj[i,'lod_delta']<0){
        # label current marker for removal
        interval_mapping_obj[i,'rmv']=T

      }   # end of decision tree for markers that were initially labeled as not to be removed
    }
  }   # end of marker loop
  
  # convert to binary, where rmv=T is 1 rmv=F is 0
  Y=(!interval_mapping_obj[,'rmv'])*1
  n=length(Y)
  # compute special matrix for numeric integration
  Ys=matrix(rep(Y,each=n),nrow=n)
  # step-wise integration of binary value
  sum_Ys=as.matrix(Ys*(lower.tri(matrix(1,n,n),diag = T)*1))%*%matrix(1,nrow=n)
  
  # for each marker that is currently flagged as sig (rmv=F)
  for(sig_mareker_i in 1: (max(sum_Ys)-1)){
    
    # locate index of marker
    (i1=sum(sum_Ys<sig_mareker_i)+1)
    
    # locate index of next sig marker
    (i2=sum(sum_Ys<(sig_mareker_i+1))+1)
    
    # if the neighborhood is spaning the same chromosome
    if(interval_mapping_obj[i1,'chr']==interval_mapping_obj[i2,'chr']){
      # extract df between i1 and i2 as neighbors
      (neighbors = interval_mapping_obj[i1:i2,])
      
      # compute maximumu lod value in this neighborhood
      (max_lod=max(neighbors$lod))
      
      # if all the neighboring lod values are within the (lod - 1) boundary
      if(sum(neighbors$lod>(max_lod-LOD_range))==nrow(neighbors)){
        # temporarely lavel all markers as to be removed
        interval_mapping_obj[i1:i2,'rmv']=T
        # relabel marker with the max_lod score as not to be removed
        interval_mapping_obj[i1:i2 & interval_mapping_obj$lod==max_lod,'rmv']=F
      }
    }
  }   # end of sig marker loop
  
  plot(interval_mapping_obj, interval_mapping_obj[!interval_mapping_obj$rmv,],
       ylim=c(0,ceiling(max(interval_mapping_obj$lod)*1.5)),
       type = c('l','p'),
       col=c('black','red'));abline(h=LOD_thr,col='red')
  
  interval_mapping_obj=interval_mapping_obj[!interval_mapping_obj$rmv,]
  
  return(interval_mapping_obj)
}