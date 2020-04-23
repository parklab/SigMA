#' The tune_new_gbm() function is used for tuning a new model 
#' for users who would like to use the algorithm for sequencing 
#' platform without built in models. See SigMA/test/test_tune_new.R
#' for an example. Additional functions for determining the
#' sensitivity, FPR, FDR and thresholds to use for the new model:
#' sen_fpr(), get_threshold(), and for adding and removing gbm
#' models to the algorithm. 
#' 
#' @param input_file the path to the file that contains the
#' data table that will be used in run.R function
#' @param rda_file the path to an rda file to save the tuned 
#' gbm model, if it is NULL function returns the model
#' @param tumor_type see list_tumor_types() for more information
#' @param data sequencing platform see list_data_options() for more
#' information
#' @param norm96 is a 96-dimensional array with weights that adjusts
#' the trinucleotide frequency in the specific sequencing platform
#' to the whole genome frequency. You can use get_trinuc_norm() for
#' calculating the weight from a bed file. If norm96 is NULL the 
#' normalization for the platform determined by data settings is used
#' @param run_SigMA boolean indicating whether the input file
#' is the output of an initial SigMA calculation (to be set to FALSE)
#' or whether it is an unprocessed file (to be set to TRUE).

tune_new_gbm <- function(input_file, 
                         rda_file = NULL, 
                         tumor_type = NULL, 
                         data = NULL, 
                         norm96 = NULL, 
                         run_SigMA = T){

 
  if(is.null(tumor_type) & !run_SigMA){
    if(!grepl('tumor_type', input_file)) stop('tumor_type not provided')
    tumor_type <- unlist(lapply(strsplit(input_file, split = 'tumortype'), 
                         function(x){x[[2]]}))
    tumor_type <- unlist(lapply(strsplit(tumor_type, split = '_'), 
                         function(x){x[[1]]}))
  }
  else if(is.null(tumor_type)){
    stop('tumor_type mising see get_tumor_types()')
  }

  if(run_SigMA){
    output_file <- run(input_file, 
                       tumor_type = tumor_type, 
                       data = data,
                       norm96 = norm96, 
                       do_assign = F,  
                       check_msi = F,
                       do_mva = F)
  }
  else{
    output_file <- input_file
  }

  gbm_model <- tune_gbm_model(output_file) 

  output_table <- predict_prob(output_file, gbm_model)  
  write.table(output_table, 
              file = gsub(input_file, 
                          pattern = '.csv', 
                          replace = '_predictions.csv'), 
              row.names = F, sep = ',', quote = F)

  if(!is.null(rda_file)) save(gbm_model, file = rda_file)
  else return(gbm_model)
}

tune_gbm_model <- function(file){
  df <- read.csv(file)
 
  features_gbm <- c(features_gbm, 'Signature_UV_c1_ml')
  features_gbm <- features_gbm[!is.na(match(features_gbm, colnames(df)))]

  gbm_model <- gbm::gbm(formula = is_sig3 ~ .,
                   distribution = 'bernoulli',
                   data = na.omit(df[df$total_snvs >= 3, c(features_gbm, 'is_sig3')]),
                   n.trees = 5000, 
                   shrinkage = 0.01,
                   bag.fraction = 0.2, 
                   n.minobsinnode= 3,
                   cv.fold = 5)
  bestTreeForPrediction = gbm::gbm.perf(gbm_model)
  gbm_model <- gbm::gbm(formula = is_sig3 ~ .,
                   distribution = 'bernoulli',
                   data = na.omit(df[df$total_snvs >= 3, c(features_gbm, 'is_sig3')]),
                   n.trees = bestTreeForPrediction, 
                   shrinkage = 0.01,
                   bag.fraction = 0.2, 
                   n.minobsinnode= 3)
  rel_infs <- gbm::relative.influence(gbm_model)
  features <- names(rel_infs)[rel_infs/sum(rel_infs) > 0.005]
  gbm_model <- gbm::gbm(formula = is_sig3 ~ .,
                   distribution = 'bernoulli',
                   data = na.omit(df[df$total_snvs >= 3, c(features, 'is_sig3')]),
                   n.trees = bestTreeForPrediction, 
                   shrinkage = 0.01,
                   bag.fraction = 0.2, 
                   n.minobsinnode= 3)
  return(gbm_model)
}

predict_prob <- function(file, gbm_model){
  df <- read.csv(file)
  prediction = predict(object = gbm_model,
                       newdata = df, 
                       n.trees = gbm_model$n.trees, 
                       type = "response")
  prediction[df$total_snvs < 5] <- 0
  df$prob <- prediction
  return(df)
}

# adds the gbm model in the system files of the package
add_gbm_model <- function(name_model, 
                          tumor_type, 
                          gbm_model = NULL, 
                          cutoff = NULL, 
                          cutoff_strict = NULL){

  file_path <- system.file("extdata/gbm_models/",
                           package="SigMA")

  if(file.exists(paste0(file_path, '/', name_model, '.rda')))
    load(paste0(file_path, '/', name_model, '.rda'))
  if(!exists('gbm_models_custom')) gbm_models_custom <- list()

  gbm_models_custom[[tumor_type]] <- gbm_model
  if(!exists('cutoffs_custom')) cutoffs_custom <- list()
  cutoffs_custom[[tumor_type]] <- cutoff
  if(!exists('cutoffs_strict_custom')) cutoffs_strict_custom <- list()
  cutoffs_strict_custom[[tumor_type]] <- cutoff_strict_custom
  save(gbm_models_custom, cutoffs_custom, cutoffs_strict_custom, file = paste0(file_path, '/', name_model, '.rda'))
}

# remove a specific model if the tumor_type is 
# specified the rest of the models are kept and just 
# the one for that tumor type is removed if is null
# all the gbm models with the given name is removed
remove_gbm_model <- function(name_model, 
                             tumor_type = NULL){
  file_path <- system.file(paste0("extdata/gbm_models/", name_model, ".rda"),
                              package="SigMA")
  if(is.null(tumor_type)){
    file.rm(file_path)
  }
  else{
    load(file_path)
    ind <- which(names(gbm_models_custom[[name_model]]) == tumor_type)
    gbm_models_custom <- gbm_models_custom[[-ind]]
    cutoffs_custom <- cutoffs_custom[[-ind]]
    cutoffs_strict_custom <- cutoffs_strict_custom[[-ind]]
    save(gbm_models_custom, cutoffs_custom, cutoffs_strict_custom, file = file_path)
  }
}

# calculates sensitivity FPR and FRD based on the parameter
# provided with the signal setting 
sen_fpr <- function(df, var, signal = 'is_sig3'){

  df <- df[!is.na(df[, var]),]
  df <- df[order(df[,var]),]
  df$signal <- df[, signal]
  #remove
  total_signal <- sum(df$signal)
  total_background <- dim(df)[[1]] - sum(df$signal == 1)

  sen <- rep(0, dim(df)[[1]])
  fpr <- rep(0, dim(df)[[1]])
  fdr <- rep(0, dim(df)[[1]])

  for(ifile in 1:dim(df)[[1]]){
    sen[[ifile]] <- (total_signal - sum(df$signal[1:ifile] == 1))
    sen[[ifile]] <- sen[[ifile]]/total_signal
    fpr[[ifile]] <- sum(df$signal[ifile:dim(df)[[1]]] != 1)
    fpr[[ifile]] <- fpr[[ifile]]/total_background
    fdr[[ifile]] <- sum(df$signal[ifile:dim(df)[[1]]] != 1)
    denom <- sum(df$signal[ifile:dim(df)[[1]]] != 1)
    denom <- denom + (total_signal - sum(df$signal[1:ifile] == 1))
    fdr[[ifile]] <- fdr[[ifile]]/denom
  }
  return(data.frame(sen = sen, fpr = fpr, fdr = fdr, value = df[,var]))
}

# get the threshold to be used to select samples with the signature
# based on the the parameter provided by the var setting.
# the signal defines the true values and cut_var defines whether the limits
# are set on sensitivity, FPR or FDR. For sensitivity the limits are applied
# as a lower bound and for FPR and FDR as an upper bound
get_threshold <- function(df, limits, var = 'prob', 
                          signal = 'is_sig3', cut_var = 'fpr'){

  df_sen_fpr <- sen_fpr(df, var, signal)
  sen_cutoff <- numeric()
  fpr_cutoff <- numeric()
  fdr_cutoff <- numeric()
  cutoff_vec <- numeric() 
  for(limit in limits){
    if(cut_var == "fpr" | cut_var == "fdr"){
      cutoff <- min(df_sen_fpr$value[df_sen_fpr[, cut_var] <= limit], na.rm = T)
      if(length(cutoff) == 0)
        limit <- min(df_sen_fpr[,cut_var], na.rm = T)
    }
    else if(cut_var == "sen"){ 
      cutoff <- max(df_sen_fpr$value[df_sen_fpr[,cut_var] >= limit], na.rm = T)
      if(length(cutoff) == 0)
        limit <- max(df_sen_fpr[,cut_var], na.rm = T)
    }
    else{
      stop('cut var can be fpr, fdr or sen')
    }  

    cutoff_vec <- c(cutoff_vec, cutoff)
    sen_cutoff <- c(sen_cutoff, 
                   max(df_sen_fpr$sen[df_sen_fpr$value > cutoff]))
    fpr_cutoff <- c(fpr_cutoff, 
                    max(df_sen_fpr$fpr[df_sen_fpr$value > cutoff]))
    fdr_cutoff <- c(fdr_cutoff, 
                    max(df_sen_fpr$fdr[df_sen_fpr$value > cutoff]))
  }
  return(list(cutoff = cutoff_vec, 
              sen = sen_cutoff, 
              fpr = fpr_cutoff, 
              fdr = fdr_cutoff))
}