pair_rois = function(df, sig_rois){
  sig_rois_paired = matrix(data = NA, nrow = nrow(df), ncol = length(sig_rois) %/% 2)
  
  for(i in 2:length(sig_rois)){
    if(i %% 2 != 0) next # skip odd numbers
    sig_rois_paired[,i/2] =  rowSums(df[,sig_rois[c(i, i+1)]])
  }
  
  sig_rois_df = data.frame(df[,sig_rois[1]], sig_rois_paired)
  colnames(sig_rois_df) = c("Brain Stem", "Cerebellum WM", "Ventral DC", "Fornix", "Ant Limb Int Capsule")
  return(sig_rois_df)
}