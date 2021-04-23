devtools::load_all()

run_classifier <- function(input_file){

  output_file <- run(input_file, data = 'op', tumor_type = 'other', cosmic_version = 'v3', do_mva = F, check_msi = T)

  # read the input file
  df <- read.csv(output_file)

  # get the columns related to exposures in most likely clusters
  df <- llh_max_characteristics(df, tumor_type = 'other', cosmic_version = 'v3')

  # combine likelihood of different clusters with same signatures  
  lite <- lite_df(df)
  df <- cbind(df, lite[match(df$tumor, lite$tumor), c('Signature_UV_ml', 'Signature_APOBEC_ml', 'Signature_pole_ml', 'Signature_4_ml', 'Signature_msi_ml')])
  df$Signature_msi_pole_ml <- df$Signature_msi_Ma_ml + df$Signature_msi_Mb_ml

  # NNLS exposures
  df_exps <- get_sig_exps(df, 'sigs_all_msi', 'exps_all_msi')
  df <- cbind(df, df_exps[,-na.omit(match(colnames(df), colnames(df_exps)))])

  df$exp_APOBEC <- 0
  if(sum(!is.na(match('exp_sig2', colnames(df))), na.rm = T) > 0)
    df$exp_APOBEC <- df$exp_sig2 
  if(sum(!is.na(match('exp_sig13', colnames(df))), na.rm = T) > 0)
    df$exp_APOBEC <- df$exp_APOBEC + df$exp_sig2 
  df$rat_APOBEC <- df$exp_APOBEC/df$total_snvs

  df$exp_sig7 <- 0
  if(sum(!is.na(match('exp_sig7a', colnames(df))), na.rm = T) > 0)
    df$exp_sig7 <- df$exp_sig7a 
  if(sum(!is.na(match('exp_sig7b', colnames(df))), na.rm = T) > 0)
    df$exp_sig7 <- df$exp_sig7 + df$exp_sig7b
  df$rat_sig7 <- df$exp_sig7/df$total_snvs

  df$rat_sig11 <- 0
  if(sum(!is.na(match('exp_sig11', colnames(df))), na.rm = T) > 0)
    df$rat_sig11 <- df$exp_sig11/df$total_snvs
  else{
    df$exp_sig11 <- 0
  }

  df$exp_pole <- 0
  if(sum(!is.na(match('exp_sig10a', colnames(df))), na.rm = T) > 0)
    df$exp_pole <- df$exp_sig10a 
  if(sum(!is.na(match('exp_sig10b', colnames(df))), na.rm = T) > 0)
    df$exp_pole <- df$exp_pole + df$exp_sig10b
  df$rat_pole <- df$exp_pole/df$total_snvs
    
  df$rat_sig4 <- 0
  if(sum(!is.na(match('exp_sig4', colnames(df))), na.rm = T) > 0)
    df$rat_sig4 <- df$exp_sig4/df$total_snvs
  else{
    df$exp_sig4 <- 0
  }

  # get NNLS exposures of maximum cluster
  df_cluster <- get_sig_exps(df, 'cluster_sigs_all', 'cluster_exps_all')  
  colnames(df_cluster) <- paste0('clust_', colnames(df_cluster))
  df <- cbind(df, df_cluster)

  df$clust_exp_APOBEC <- 0
  if(sum(!is.na(match('clust_exp_sig2', colnames(df))), na.rm = T) > 0)
    df$clust_exp_APOBEC <- df$clust_exp_sig2
  if(sum(!is.na(match('clust_exp_sig13', colnames(df))), na.rm = T) > 0)
    df$clust_exp_APOBEC <- df$clust_exp_APOBEC + df$clust_exp_sig13

  df$clust_exp_UV <- 0
  if(sum(!is.na(match('clust_exp_sig7a', colnames(df))), na.rm = T) > 0)
    df$clust_exp_UV <- df$clust_exp_sig7a
  if(sum(!is.na(match('clust_exp_sig7b', colnames(df))), na.rm = T) > 0)
    df$clust_exp_UV <- df$clust_exp_UV  + df$clust_exp_sig7b
  
  df$clust_exp_msi <- 0
  if(sum(!is.na(match(signames_per_tissue[['msi']], colnames(df)))) > 1)
    df$clust_exp_msi <- rowSums(df[, na.omit(match(gsub(signames_per_tissue[['msi']], pattern = 'Signature_', replace ='clust_exp_sig'), colnames(df)))])
  else if(sum(!is.na(match(signames_per_tissue[['msi']], colnames(df)))) == 1)
    df$clust_exp_msi <- df[, na.omit(match(gsub(signames_per_tissue[['msi']], pattern = 'Signature_', replace ='clust_exp_sig'), colnames(df)))]
  else 
    df$clust_exp_msi <- 0 

  df$exp_msi <- 0
  if(sum(!is.na(match(signames_per_tissue[['msi']], colnames(df)))) > 1)
    df$exp_msi <- rowSums(df[, na.omit(match(gsub(signames_per_tissue[['msi']], pattern = 'Signature_', replace ='exp_sig'), colnames(df)))])
  else if(sum(!is.na(match(signames_per_tissue[['msi']], colnames(df)))) == 1)
    df$exp_msi <- df[, na.omit(match(gsub(signames_per_tissue[['msi']], pattern = 'Signature_', replace ='exp_sig'), colnames(df)))]
  df$rat_msi <- df$exp_msi/df$total_snvs

  df$clust_exp_pole <- 0
  if(sum(!is.na(match(signames_per_tissue[['pole']], colnames(df)))) > 1)
    df$clust_exp_pole <- rowSums(df[, na.omit(match(gsub(signames_per_tissue[['pole']], pattern = 'Signature_', replace ='clust_exp_sig'), colnames(df)))])
  else if(sum(!is.na(match(signames_per_tissue[['pole']], colnames(df)))) == 1)
    df$clust_exp_pole <- df[, na.omit(match(gsub(signames_per_tissue[['pole']], pattern = 'Signature_', replace ='clust_exp_sig'), colnames(df)))]
  else 
    df$clust_exp_pole <- 0
  
  prediction = predict(object = gbm_model_OP,
                       newdata = df,
                       n.trees = gbm_model_OP$ntrees,
                       type = "response")

  df_pred <- data.frame(prediction)
  num_tree_str <- strsplit(colnames(df_pred)[grepl('is_uv', colnames(df_pred))], split = '\\.')[[1]][[2]]
  colnames(df_pred) <- paste0('prob_', gsub(colnames(df_pred), pattern = paste0('\\.', num_tree_str), replace = ''))
 
  df <- cbind(df, df_pred)
  df <- get_classes(df, 
                    composite_classes = c('is_tobacco.is_apobec', 'is_temozolomide.is_msi','is_pole.is_msi'), 
                    single_classes = c('is_tobacco', 'is_apobec', 'is_temozolomide','is_msi', 'is_pole', 'is_uva', 'none'),
                    thresholds_single = c(0.2, 0.8), composite_thresholds = c(0.2, 0.5),
                    exception_thresh_high = c(0.2, 0.2), exception_class_high = c('is_msi', 'is_pole')) 

  return(df)
}

get_classes <- function(df, 
                        composite_classes, 
                        single_classes,
                        thresholds_single,
                        composite_thresholds,
                        exception_thresh_high = NULL,
                        exception_class_high = NULL,
                        exception_thresh_low = NULL, 
                        exception_class_low = NULL
                        ){

  for(class_this in single_classes){
    df[,paste0(class_this,'_predict')] <- as.integer(df[,paste0('prob_', class_this)] > thresholds_single[[1]]) + as.integer(df[,paste0('prob_', class_this)] > thresholds_single[[2]])
  }
  for(class_this in composite_classes){
    class1 <- strsplit(class_this, split = '\\.')[[1]][[1]]
    class2 <- strsplit(class_this, split = '\\.')[[1]][[2]]
    inds <- which(df[,paste0('prob_', class1)] + df[,paste0('prob_', class2)] + df[,paste0('prob_', class_this)] > 0.9 & df[,paste0('prob_', class_this)] > composite_thresholds[[1]])
    inds_zero <- which(df[, paste0(class1, '_predict')] == 0)
    df[inds[inds %in% inds_zero], paste0(class1, '_predict')] <- 1
    inds_zero <- which(df[, paste0(class2, '_predict')] == 0)
    df[inds[inds %in% inds_zero], paste0(class2, '_predict')] <- 1
    inds <- which(df[,paste0('prob_', class_this)] > composite_thresholds[[2]])
    df[inds, paste0(class1, '_predict')] <- 2
    df[inds, paste0(class2, '_predict')] <- 2
    inds <- which(df[,paste0('prob_', class1)]  < 0.005 & df[,paste0(class1,'_predict')] != 0)
    df[inds, paste0(class1, '_predict')] <- 0
    inds <- which(df[,paste0('prob_', class2)]  < 0.005 & df[,paste0(class2,'_predict')] != 0)
    df[inds, paste0(class2, '_predict')] <- 0
  }
  if(!is.null(exception_class_high)){
    for(i in 1:length(exception_class_high)){
      class_this <- exception_class_high[[i]]
      thresh <- exception_thresh_high[[i]]
      df[df[,paste0('prob_', class_this)] > thresh, paste0(class_this, '_predict')] <- 2
    }
  }
  if(!is.null(exception_class_low)){
    for(i in 1:length(exception_class_low)){
      class_this <- exception_class_low[[i]]
      thresh <- exception_thresh_low[[i]]
      df[df[,paste0('prob_', class_this)] > thresh & df[,paste0(class_this,'_predict')] == 0, paste0(class_this, '_predict')] <- 1
    }
  }
  return(df)
}

get_info <- function(df, sample, thresh_unique = 0.0, msi_col = 'is_msi_predict', pole_col = 'is_pole_predict'){
  classes_full_names <- list()
  classes_full_names[['is_pole']] <- 'POLE-exo domain mutations'
  classes_full_names[['is_msi']] <- 'Mismatch repair deficiency'
  classes_full_names[['is_apobec']] <- 'Activity of APOBEC cytidine deaminases'
  classes_full_names[['is_temozolomide']] <- 'Temozolomide'
  classes_full_names[['is_tobacco']] <- 'Tobacco smoke'
  classes_full_names[['is_uv']] <- 'UV radiation'
  classes_full_names[['none']] <- 'No POLE/MMRD/APOBEC/Temozol/Tobacco/UV'

  this <- df[df$tumor == sample,]
  classes <- gsub(colnames(this)[grep('predict', colnames(this))], pattern = '_predict', replace ='')
  class_names <- unlist(classes_full_names)[match(classes, names(classes_full_names))]

  vals <- this[, paste0(classes, '_predict')]
  if(sum(vals == 2, na.rm = T) > 0)
    message(paste0('Signatures that are present:', paste0(unlist(class_names)[vals == 2], collapse = ', '),collapse = ''))
  if(sum(vals == 1, na.rm = T) > 0)
    message(paste0('Signatures that are present (L.C.):', paste0(unlist(class_names)[vals == 1], collapse = ', '), collapse = ''))
 
  signames_clust <- gsub(colnames(this)[grep('clust_exp_sig', colnames(this))], pattern = 'clust_exp_sig', replace = 'SBS')
  vals <- this[, grep('clust_exp_sig', colnames(this))]
  message(paste0('Likelihood of most likely cluster:', round(this$max_ml_msi, digit = 3)))
  message(paste0('Most likely cluster has:', paste0(paste0(signames_clust[vals > 0], ' ', round(vals[vals > 0]*100, digit = 0), '%'), collapse = ',')))

  # get columns with signature names, exposures and likelihood ratios
  signames <- gsub(colnames(this)[grepl('exp_sig', colnames(this)) & !grepl('clust', colnames(this))], pattern = 'exp_sig', replace = 'SBS')
  vals <- this[, gsub(signames, pattern = 'SBS', replace = 'exp_sig')]/this$total_snvs
  cols_l_rat <- paste0(gsub(signames, pattern = 'SBS', replace = 'Signature_'), '_l_rat_msi')
  inds <- which(is.na(match(cols_l_rat, colnames(this))))
  if(length(inds) > 0){
    signames <- signames[-inds]
    vals <- vals[-inds] 
    cols_l_rat <- cols_l_rat[-inds]
  }
  lrats = this[,cols_l_rat]
  
  # get signatures with non-zero exposures
  lrats <- lrats[vals > 0]
  signames <- signames[vals > 0]
  vals <- vals[vals > 0]

  signames <- signames[lrats > thresh_unique] 
  vals <- vals[lrats > thresh_unique]
  lrats <- lrats[lrats > thresh_unique]  

  # order with the largest contribution reported first
  order_inds <- order(-vals)
  lrats <- lrats[order_inds]
  signames <- signames[order_inds]
  vals <- vals[order_inds]
  
  message(paste0('NNLS result is:', paste0(paste0(signames, ' ', round(vals*100, digit = 0), '% - uniqueness score:', round(lrats, digit = 2), collapse = ',')))) 

  if(this[,msi_col] == 0 & this[,pole_col] == 0){
    df_exps <- get_sig_exps(this, 'sigs_all', 'exps_all')

    signames <- gsub(colnames(df_exps)[grepl('exp_sig', colnames(df_exps)) & !grepl('clust', colnames(df_exps))], pattern = 'exp_sig', replace = 'SBS')
    vals <- df_exps[, gsub(signames, pattern = 'SBS', replace = 'exp_sig')]/this$total_snvs
    cols_l_rat <- paste0(gsub(signames, pattern = 'SBS', replace = 'Signature_'), '_l_rat')
    inds <- which(is.na(match(cols_l_rat, colnames(this))))
    if(length(inds) > 0){
      signames <- signames[-inds]
      vals <- vals[-inds]
      cols_l_rat <- cols_l_rat[-inds]
    }
    lrats = this[,cols_l_rat]

    # get signatures with non-zero exposures
    lrats <- lrats[vals > 0]
    cols_l_rat <- cols_l_rat[vals > 0]
    signames <- signames[vals > 0]
    vals <- vals[vals > 0]

    cols_l_rat <- cols_l_rat[lrats > thresh_unique]
    signames <- signames[lrats > thresh_unique]
    vals <- vals[lrats > thresh_unique]
    lrats <- lrats[lrats > thresh_unique]

    order_inds <- order(-vals)
    lrats <- lrats[order_inds]
    signames <- signames[order_inds]
    vals <- vals[order_inds]

    message(paste0('NNLS result without MMRD and POLE-exo signatures is:', paste0(paste0(signames, ' ', round(vals*100, digit = 0), '% - uniqueness score:', round(lrats, digit = 2), collapse = ','))))
  }
}