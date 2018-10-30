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
                       do_strict = T,
                       weight_cf){
  
  if(method == 'median_catalog'){
    matching = 'ml'
    if(length(grep(signame, colnames(df_in))) > 1)
      pass <- rowSums(df_in[, grep(signame, colnames(df_in))]) > 0.5
    else
      pass <- df_in[, grep(signame, colnames(df_in))] > 0.5
  }
  else if(method == 'mva'){
    matching = 'mva'

    if(data == "msk"){
      if(weight_cf){
        cutoffs_strict <- cutoffs_msk_strict_cf
        cutoffs <- cutoffs_msk_cf
      }else{
        cutoffs_strict <- cutoffs_msk_strict
        cutoffs <- cutoffs_msk     
      }
    }else if(data == "seqcap"){
      cutoffs_strict <- cutoffs_exome_strict
      cutoffs <- cutoffs_exome 
    }else if(data == "wgs"){ 
      cutoffs_strict <- cutoffs_wgs_strict
      cutoffs <- cutoffs_wgs
    }else{
      stop('invalid data selection')
    }
    pass <- (df_in[, paste0(signame, '_mva')] >= cutoffs[[tumor_type]])     
    pass_strict <- (df_in[, paste0(signame, '_mva')] >= cutoffs_strict[[tumor_type]])     
  }
  else{
    stop('assignments cannot be done for this method')
  }

  if(sum(colnames(df_in) == 'total_snvs') == 0)
    df_in$total_snvs <- rowSums(df_in[, 1:96])
 
  df_out <- data.frame(pass = pass)
  colnames(df_out)[[1]] <- paste0('pass_', matching)

  if(method == 'mva'){
    df_out$pass_mva_strict <- pass_strict
  }

  return(df_out)
}
