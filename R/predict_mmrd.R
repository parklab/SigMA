#' Predicts the MVA scores for three classes (MMRP, MMRD and POLE)
#' and categorizes samples. The new info is added as new columns
#' to the input data table, (MMRP_prob, MMRD_prob, POLE_prob, MMRD_categ,
#' respectively). Also the most likely cluster and subtype are determined.
#' Should be used after the run() function with check_msi = TRUE settings.
#' An example use case can be found in examples/run_SigMA_MMRD.R
#' 
#' @param input_file file path containing indel features and the 
#' signature calculation results by SigMA.
#' @param data determines the type of sequencing platform
#' used for dataset see list_data_options() and find_data_setting()
#' @param df alternatively a data.frame containing the same info as 
#' input_file can be given as input. Set to NULL by default. If df is 
#' not null and the output file is not saved and instead the new data
#' frame is returned


get_exps_mmrd <- function(df){
  sigs_mmrd <- c(signames_per_tissue$msi, signames_per_tissue$pole)
  exp_list <- apply(df, 1, function(x, sigs_mmrd){
    exps <- c(as.numeric(unlist(strsplit(as.character(x['exps_all_msi']), split = '\\_'))))
    sigs <- c(unlist(strsplit(as.character(x['sigs_all_msi']), split = '\\.')))
    
    exps_mmrd <- rep(0, length(sigs_mmrd))
    names(exps_mmrd) <- sigs_mmrd
   
    exps_mmrd[!is.na(match(sigs_mmrd, sigs))] <- exps[na.omit(match(sigs_mmrd, sigs))]
    return(exps_mmrd)
  }, sigs_mmrd = sigs_mmrd)
 
  df_exp <- data.frame(t(exp_list))
  colnames(df_exp) <- gsub(colnames(df_exp), pattern = 'Signature_', replace = 'exp_sig')

  df_exp$exp_sig10 <- df_exp$exp_sig10a + df_exp$exp_sig10b

  return(cbind(df, df_exp))
}

predict_mmrd <- function(data, input_file = NULL, df = NULL){

  if(is.null(data)){
    stop('data parameter needs to be set')
  }
  if(!is.null(input_file))
    df <- read.csv(input_file)

  if(sum(colnames(df) == "msisensor") > 0) tag <- 'with_msisensor'
  else tag <- 'without_msisensor'

  if(sum(grepl('msi', colnames(df)) & grepl('Signature_', colnames(df))) == 0)  stop('first run SigMA with check_msi = T parameters')

  model <- gbm_models_mmrd[[tag]][[data]]

  if(sum(colnames(df) == "Signature_msi_ml") == 0)
    df$Signature_msi_ml <- rowSums(df[,grepl('Signature_msi_', colnames(df)) & grepl('_ml', colnames(df))])
  if(sum(colnames(df) == "Signature_pole_ml") == 0)
    df$Signature_pole_ml <- rowSums(df[,grepl('Signature_pole_', colnames(df)) & grepl('_ml', colnames(df))])
  if(sum(colnames(df) == "Signature_msi_pole_ml") == 0)
    df$Signature_msi_pole_ml <- rowSums(df[,grepl('Signature_msi_La|Signature_msi_Lb|Signature_msi_Ma|Signature_msi_Mb', colnames(df)) & grepl('_ml', colnames(df))])

  if(sum(colnames(df) %in% c('exp_sig6','exp_sig15','exp_sig20','exp_sig21','exp_sig26')) == 0)
    df <- get_exps_mmrd(df)
    
  terms <- attributes(gbm_models_mmrd$without_msisensor$msk$Terms)$term.labels

  missing_terms <- terms[!(terms %in% colnames(df))]
  for(col in missing_terms){
    df$new <- 0
    colnames(df)[colnames(df) == "new"] <- col
  }
  
  predictions = predict(object = model,
                        newdata = df[, names(gbm::relative.influence(model))],
                        n.trees = model$n.trees,
                        type = "response")
  df$MMRP_prob <- predictions[1:(length(predictions)/3)]
  df$POLE_prob <- predictions[(length(predictions)/3+1):(2*length(predictions)/3)]
  df$MMRD_prob <- predictions[(2*length(predictions)/3 + 1):length(predictions)]

  df$MMRD_category <- 'Unclassified'
  df$MMRD_category[df$MMRP_prob > cutoffs_mmrd[[tag]][[data]][['MMRP']]] <- 'MMRP'
  df$MMRD_category[df$POLE_prob > cutoffs_mmrd[[tag]][[data]][['POLE']]] <- 'POLE'
  df$MMRD_category[df$MMRD_prob > cutoffs_mmrd[[tag]][[data]][['MMRD']]] <- 'MMRD'

  t <- df[,grepl('Signature_pole|Signature_msi', colnames(df))]
  t <- t[,!(colnames(t) %in% c('Signature_msi_ml', 'Signature_pole_ml', 'Signature_msi_pole_ml'))]

  df$max_clust <- unlist(lapply(apply(t, 1, function(x, cols){cols[which(x == max(x))]}, cols = colnames(t)), function(x){x[[1]]}))
  df$max_clust <- gsub(gsub(df$max_clust, pattern = 'Signature_msi_', replace = ''), pattern = '_ml_msi', replace = '')
  df$max_clust[df$max_clust == "Signature_pole_c1"] <- "pole"
  df$max_clust[df$max_clust == "Signature_pole_c2"] <- "pole"
  df$max_clust[df$MMRD_category == "MMRP"] <- 'MMRP'

  df$subtype_4 <- df$Signature_pole_c1_ml_msi + df$Signature_pole_c2_ml_msi
  df$subtype_3 <- df$Signature_msi_La_ml_msi + df$Signature_msi_Lb_ml_msi
                   + df$Signature_msi_Ma_ml_msi + df$Signature_msi_Mb_ml_msi
  df$subtype_2 <- rowSums(df[ ,na.omit(match(paste0('Signature_msi_', c('F', 'G', 'H', 'I', 'J'), '_ml_msi')
                                           , colnames(df)))])
  df$subtype_1 <- rowSums(df[ ,na.omit(match(paste0('Signature_msi_', c('A', 'B', 'C', 'D', 'E'), '_ml_msi')
                                           , colnames(df)))])

  t <- df[, paste0('subtype_', 1:4)]
  df$max_subtype <- unlist(lapply(apply(t, 1, function(x, cols){cols[x == max(x)]}, cols = colnames(t)), function(x){x[[1]]}))
  df <- df[,!(colnames(df) %in% paste0('subtype_', 1:4))]
  df$max_subtype[df$MMRD_category == "MMRP"] <- 'MMRP'

  output_file <- gsub(input_file, pattern = '.csv', replace = paste0('_mmrd_', data, '.csv'))

  message(paste0('The MMRD predictions are in: ', output_file))
   
  if(!is.null(input_file)){
    write.table(df, 
      file = output_file,
      row.names = F, 
      sep = ',', 
      quote = F)
  }
  else
    return(df)
}

