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
                       cutoffs = NULL, 
                       custom = F){

  if(method == 'median_catalog'){
    matching = 'ml'
    if(length(grep(signame, colnames(df_in))) > 1)
      pass <- rowSums(df_in[, grep(signame, colnames(df_in))]) > 0.5
    else
      pass <- df_in[, grep(signame, colnames(df_in))] > 0.5
    df_out <- data.frame(pass_ml = pass)
  }
  if(sum(grepl('tcga_mc3|msk|seqcap|seqcap_probe|wgs|wgs_pancan', data)) > 0){
    if(method == 'mva'){
      matching = 'mva'
      if(data == "msk"){
        if(weight_cf){
          cutoffs_strict <- cutoffs_msk_strict_cf
          cutoffs <- cutoffs_msk_cf
        }else{
          cutoffs_strict <- cutoffs_msk_strict
          cutoffs <- cutoffs_msk     
        }
      }else if(data == "tcga_mc3"){
        cutoffs_strict <- cutoffs_mc3_strict
        cutoffs <- cutoffs_mc3 
      }else if(data == "seqcap"){
        cutoffs_strict <- cutoffs_exome_strict
        cutoffs <- cutoffs_exome 
      }else if(data == "seqcap_probe"){
        cutoffs_strict <- cutoffs_seqcap_probe_strict
        cutoffs <- cutoffs_seqcap_probe
      }else if(data == "wgs_pancan"){ 
        cutoffs_strict <- cutoffs_wgs_strict
        cutoffs <- cutoffs_wgs
      }else if(data == "wgs"){ 
        cutoffs_strict <- cutoffs_wgs_tumor_type_strict
        cutoffs <- cutoffs_wgs_tumor_type
      }else if(custom){
        file_path <- system.file(paste0("extdata/gbm_models/", data, ".rda"),
                             package="SigMA")
        load(file_path)
      }
      else{
        stop('invalid data selection')
      }
      if(is.null(tumor_type)) stop('tumor_type not provided')
      pass <- (df_in[, paste0(signame, '_mva')] >= cutoffs[[tumor_type]])     
      pass_strict <- (df_in[, paste0(signame, '_mva')] >= cutoffs_strict[[tumor_type]])

      df_out <- data.frame(pass_mva = pass)
      df_out$pass_mva_strict <- pass_strict
    }
  }

  if(!is.null(cutoffs) & !is.null(limits) & !is.null(cut_var)){
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
        df_out$new <- df_in[,variable] > cutoffs[[i]] 
      else
        df_out <- data.frame(new = df_in[,variable] > cutoffs[[i]])
      colnames(df_out)[colnames(df_out) == "new"] <- paste0('pass_', matching, '_', cut_var, '_', round(limits[[i]], digit = 2))
    }
  }

  if(sum(colnames(df_in) == 'total_snvs') == 0)
    df_out$total_snvs <- rowSums(df_out[, 1:96])
 
  return(df_out)
}
