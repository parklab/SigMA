#' Assigns a boolean based on a threshold on the 
#' likelihood or mva score for whether the signature
#' is identified 
#'
#' @param df_in input data.frame 
#' @param method 'median_catalog' for likelihood based
#' selection or 'mva' for multivariate analysis score 
#' based selection
#' @param signame name of the signature that user wants 
#' to identify, 'Signature_3' or 'Signature_msi'
#' @param data 'msk', 'seqcap' or 'wgs'
#' @param tumor_type tumor type as listed in 
#' https://github.com/parklab/SigMA/ because the 
#' thresholds are tumor_type specific
#' @param do_strict sets whether a strict 
#' threshold should be applied or a loose one
#'
#' @return a data.frame with a single column which contains 
#' the boolean indicating the presence of the signature

assignment <- function(df_in,
                       method = 'mva',
                       signame = 'Signature_3',
                       data = NULL, 
                       tumor_type = "breast",
                       do_strict = T){
  
  if(method == 'median_catalog') matching = 'ml'
  if(method == 'mva') matching = 'mva'

  if(data == "msk"){
    if(do_strict) cutoffs <- cutoffs_msk_strict
    else cutoffs <- cutoffs_msk 
  }else if(data == seqcap){
    if(do_strict) cutoffs <- cutoffs_exome_strict
    else cutoffs <- cutoffs_exome 
  }

  if(sum(colnames(df_in) == 'total_snvs') == 0)
    df_in$total_snvs <- rowSums(df_in[, 1:96])
 

  if(method == 'median_catalog'){
    if(length(grep(signame, colnames(df_in))) > 1)
      pass <- rowSums(df_in[, grep(signame, colnames(df_in))]) > 0.5
    else
      pass <- df_in[, grep(signame, colnames(df_in))] > 0.5
  }
  if(method == 'mva'){ 
    if(data == "msk" | data == "seqcap") 
      pass <- (df_in[, paste0(signame, '_mva')] >= cutoffs[[tumor_type]]) #0.45/#0.4    
    else if(data == "found") pass <- (df_in[,paste0(signame, '_mva')] >= 0.45) #0.5 
    else pass <- (df_in[, paste0(signame, '_mva')] >= 0.3) #0.05
  }

  df_out <- data.frame(pass = pass)
  colnames(df_out)[[1]] <- paste0('pass_', matching)
 
  return(df_out)
}
