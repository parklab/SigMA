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

assignment <- function(input_files,
                       tune_df = NULL,
                       method = 'weighted_catalog',
                       signame = 'Signature_3'){

  
  if(exists('df_out')) rm(df_out)
  nfiles <- length(input_files)
  
  for(ifile in 1:nfiles){
    df_in <- read.csv(input_files[[1]])
    if(method == 'median_catalog') matching = 'ml'
    if(method == 'weighted_catalog') matching = 'wl'
    if(method == 'cosine_simil') matching = 'c'

    
    indices <- which(colnames(df_in) == paste0(signame, '_', matching))

    if(sum(colnames(df_in) == 'total_snv') == 0)
      df_in$total_snvs <- rowSums(df_in[, 1:96])

    df_in <- df_in[, c('total_snvs', colnames(df_in)[indices])]

    if(method == 'weighted_catalog' | method == 'median_catalog'){ 
      pass <- assign(df_in,                       
                     cuts = tune_df[, 'cutoff_l'],
                     snvs = tune_df$snv_low)
    }
    if(method == 'cosine_simil'){ 
      pass <- assign(df_in,                       
                     cuts = tune_df[, 'cutoff_c'],
                     snvs = tune_df$snv_low)
    }

    df_out_this <- data.frame(pass = pass)
    colnames(df_out_this)[[1]] <- paste0('pass_', matching)

    if(exists('df_out')) df_out <- rbind(df_out, df_out_this)
    else df_out <- df_out_this
  }
  return(df_out)
}


assign <- function(df_m, cuts, snvs){
  results <- rep(FALSE, dim(df_m)[[1]])
  indices <- match(df_m$total_snvs, snvs)
  
  if(sum(is.na(indices)) > 0)
    stop('total snv in sample is out of range of the available tune')

  results <- (df_m[, 2] >= cuts[indices])
  return(results)
}