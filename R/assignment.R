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
#' @param data 'msk', 'tcga_mc3', 'seqcap' 'wgs' 
#' or 'wgs_pancan'
#' @param tumor_type tumor type as listed in 
#' https://github.com/parklab/SigMA/ because the 
#' thresholds are tumor_type specific
#' @param cutoffs the thresholds to select the samples
#' with the signature
#' @param cut_var 'fpr', 'fdr' or 'sen' the variable
#' that was used to determine cutoffs, only used here
#' to define the column name
#' @param limits the limits that were used on the 
#' variable set by cut_var to determine the cutoffs
#' only used for naming the columns not used in calculation
#'
#' @return a data.frame with a single column which contains 
#' the boolean indicating the presence of the signature


assignment <- function(df_in,
                       method = 'mva',
                       signame = 'Signature_3',
                       data = NULL, 
                       tumor_type = NULL,
                       weight_cf = F,
                       cut_var = NULL,
                       limits = NULL, 
                       cutoffs_custom = NULL, 
                       custom_model = F){

  limits_custom <- F
  if(method == 'median_catalog'){
    matching = 'ml'
    if(length(grep(signame, colnames(df_in))) > 1)
      pass <- rowSums(df_in[, grep(signame, colnames(df_in))]) > 0.5
    else
      pass <- df_in[, grep(signame, colnames(df_in))] > 0.5
    df_out <- data.frame(pass_ml = pass)
  }
  if(method == 'mva'){
    if((data %in% names(cutoffs)) | custom_model){
      if(is.null(cutoffs_custom)){
        matching = 'mva'
        if(data == "msk" & weight_cf){
          cutoffs_strict_this <- cutoffs_msk_strict_cf[[tumor_type]]
          cutoffs_this <- cutoffs_msk_cf[[tumor_type]]
          if(!(tumor_type %in% names(cutoffs_this))) stop('tumor_type not available for the data setting')

        }
        else if(custom_model){
          file_path <- system.file(paste0("extdata/gbm_models/", data, ".rda"),
                               package="SigMA")
          load(file_path)
        
          if(!(tumor_type %in% names(cutoffs_custom))) stop('tumor_type not available for the data setting')
          cutoffs_this <- cutoffs_custom[[tumor_type]]
          cutoffs_strict_this <- cutoffs_strict_custom[[tumor_type]]
        }
        else if(data %in% names(cutoffs)){
          cutoffs_this <- cutoffs[[data]][[tumor_type]]
          cutoffs_strict_this <- cutoffs_strict[[data]][[tumor_type]]
        }
        else{
          stop('invalid data selection')
        }
      }
      else{
        if(sum(is.na(match(limits, c(0.01,0.1)))) == 0){
          cutoffs_strict_this <- cutoffs_custom[limits == 0.01]  
          cutoffs_this <- cutoffs_custom[limits == 0.1]
        }else{
          limits_custom <- T
        }
      }
      if(!limits_custom){
        pass <- (df_in[, paste0(signame, '_mva')] >= cutoffs_this)     
        pass_strict <- (df_in[, paste0(signame, '_mva')] >= cutoffs_strict_this)

        df_out <- data.frame(pass_mva = pass)
        df_out$pass_mva_strict <- pass_strict
      }
    }
  }

  if(limits_custom){
    if(method == "mva"){
      matching <- 'mva'
      variable <- paste0(signame, '_', matching)
    }
    if(method == "median_catalog"){
      matching <- 'ml'
      variable <- paste0(signame, '_', matching)
    }
    for(i in 1:length(limits)){
      if(exists('df_out'))
        df_out$new <- df_in[,variable] > cutoffs_custom[[i]] 
      else
        df_out <- data.frame(new = df_in[,variable] > cutoffs_custom[[i]])
      colnames(df_out)[colnames(df_out) == "new"] <- paste0('pass_', matching, '_', cut_var, '_', round(limits[[i]], digit = 2))
    }
  }

  if(sum(colnames(df_in) == 'total_snvs') == 0)
    df_out$total_snvs <- rowSums(df_out[, 1:96])
 
  return(df_out)
}
