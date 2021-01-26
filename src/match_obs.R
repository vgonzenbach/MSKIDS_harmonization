match_obs = function(ref_id, target_id){
  #' Matches two vectors of IDs where vectors may have different formats
  
  row_matches = vector(mode = "integer", length = length(ref_id))
  unmatched = vector(mode = "integer")
  duplicates = list(ref_id = vector(mode = "character"), dup_targets = list())
  
  # This line not generalizable: change input to list of characters
  id_tokens = ref_id %>% sapply(function(x) strsplit(x, "-")) %>% as.data.frame() %>% unname() %>% t()
  
  # Cross-reference tokenized ID
  for(i in 1:nrow(id_tokens)){
    token_match = sapply(id_tokens[i,], grepl, target_id) # a temp array
    row_match = which(rowSums(token_match) == ncol(id_tokens))
    
    if(length(row_match) == 1){
      row_matches[i] = row_match
    }else if(length(row_match) == 0){
      unmatched = c(unmatched, ref_id[i])
      row_matches[i] = NA
    }else if(length(row_match) > 1){
      duplicates$ref_id = c(duplicates$ref_id, ref_id[i])
      duplicates$dup_targets = c(duplicates$dup_targets, list(target_id[row_match])) # OR row_match
      row_matches[i] = NA
    }
    i = i + 1
  }
  
  res = list(row_matches = row_matches, unmatched = unmatched, duplicates = duplicates)
  return(res)
}