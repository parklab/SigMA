#' This code assigns a boolean to the genome based on
#' the tuned cutoff values as a function of number of 
#' SNVs and the likelihood values calculated in advance 
#' with match_to_catalog.R function
#'
#' @param input_files a vector of filenames that was 
#' produced as the output of the run.R
#' @param tune_df a data frame with the tune, i.e. 
#' threshold values as a function of number of snvs 
#' @param method 'weighted_catalog' or 'median_catalog'
#' which ever is used in the tuning
#' @param signame name of the signature as in the catalog
#' which is aimed to be detected and with respect to 
#' which the tuning should be set, with the default 
#' settings the user can use 'Signature_3', or 'Signature_8'
#' options
#'
#' @returns a data vector which contains the pass boolean
#' for the likelihood and cosine similarity thresholds of
#' the given method of choice

assignment <- function(df_in,
                       tune_df = NULL,
                       method = 'median_catalog',
                       signame = 'Signature_3',
                       data = NULL){
  if(exists('df_out')) rm(df_out)  
  
  if(method == 'median_catalog') matching = 'ml'
  if(method == 'gbm') matching = 'gbm'


  if(sum(colnames(df_in) == 'total_snvs') == 0)
    df_in$total_snvs <- rowSums(df_in[, 1:96])
 

  if(method == 'median_catalog'){
    if(length(grep(signame, colnames(df_in))) > 1)
      pass <- rowSums(df_in[, grep(signame, colnames(df_in))]) > 0.5
    else
      pass <- df_in[, grep(signame, colnames(df_in))] > 0.5
  }
  if(method == 'gbm'){ 
    if(data == "msk") pass <- (df_in[, paste0(signame, '_gbm')] >= 0.4) #0.45                     
    else if(data == "found") pass <- (df_in[,paste0(signame, '_gbm')] >= 0.45) #0.5 
    else if(data == "seqcap") pass <- (df_in[, paste0(signame, '_gbm')] >= 0.4) #0.4 
    else pass <- (df_in[, paste0(signame, '_gbm')] >= 0.3) #0.05
  }

  df_out_this <- data.frame(pass = pass)
  colnames(df_out_this)[[1]] <- paste0('pass_', matching)

  if(exists('df_out')) df_out <- rbind(df_out, df_out_this)
  else df_out <- df_out_this
}
