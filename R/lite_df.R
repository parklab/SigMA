#' produces a data.frame with fewer columns for easier use
#'
#' @param merged_output is the input data.frame

lite_df <- function(merged_output){
  lite <- merged_output[, c('tumor', 'total_snvs')]
  merged_output$exps_all <- as.character(merged_output$exps_all)
  merged_output$sigs_all <- as.character(merged_output$sigs_all)
  
  if(sum(colnames(merged_output) == "Signature_3_ml") > 0) 
     merged_output <- merged_output[, -which(colnames(merged_output) ==  "Signature_3_ml")]

  all_cols <- colnames(merged_output)

  # likelihood columns to add up for different categories

  cols_ml <- all_cols[grepl('Signature_', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols)]
  cols_ml_msi <- all_cols[grepl('Signature_', all_cols) & grepl('_ml', all_cols) & grepl('_msi', all_cols)]

  groups <- unique(unlist(lapply(strsplit(cols_ml, split = '_'), function(x){x[[2]]})))
  groups_msi <- unique(unlist(lapply(strsplit(cols_ml_msi, split = '_'), function(x){x[[2]]})))
  groups_msi <- groups_msi[is.na(match(groups_msi, groups))]

  for(group in groups){
    cols_this <- cols_ml[grepl(paste0('Signature_', group), cols_ml)] 
    if(length(cols_this) > 1)
      lite$new <- rowSums(merged_output[,cols_this]) 
    else 
      lite$new <- merged_output[,cols_this]
    colnames(lite)[colnames(lite) == "new"] <- paste0('Signature_', group, '_ml')
  }
  for(group in groups_msi){
    cols_this <- cols_ml_msi[grepl(paste0('Signature_', group), cols_ml_msi)] 
    lite$new <- rowSums(merged_output[,cols_this]) 
    colnames(lite)[colnames(lite) == "new"] <- paste0('Signature_', group, '_ml')
  }

  if(sum(colnames(merged_output) == "MMRD_category") > 0)
    lite$MMRD_category <- merged_output$MMRD_category
  
  # signature 3 related values
  col_sig3_c <- na.omit(match('Signature_3_c', all_cols))
  col_exp_sig3 <- na.omit(match('exp_sig3', all_cols))
  col_sig3_l_rat <- na.omit(match('Signature_3_l_rat', all_cols))
  col_sig3_mva <- na.omit(match('Signature_3_mva', all_cols))
  col_pass_mva <- na.omit(match('pass_mva', all_cols))
  col_pass_mva_strict <- na.omit(match('pass_mva_strict', all_cols))
  col_pass_ml <- na.omit(match('pass_ml', all_cols))

  # nnls results
  col_sigs_all <- na.omit(match('sigs_all', all_cols))
  col_exps_all <- na.omit(match('exps_all', all_cols))

  inds_keep <- c(col_sig3_c,
                 col_exp_sig3,
                 col_sig3_l_rat,
                 col_sig3_mva,
                 col_pass_mva,
                 col_pass_mva_strict,
                 col_sigs_all,
                 col_exps_all)
  if(sum(colnames(lite) == "pass_mva") == 0) inds_keep <- c(inds_keep, col_pass_ml)

  sum_cols <- function(col_indices){
    return(unlist(apply(merged_output, 1,
                        function(x, inds){
                          sum(as.numeric(x[inds]))
                        },
                        inds = col_indices)))
  }

  if(length(inds_keep) > 0)
    lite <- cbind(lite, merged_output[, inds_keep])

  cols_ml <- colnames(lite)[grepl('_ml', colnames(lite)) & !grepl('pass', colnames(lite))]
  
  categs <- rep('', dim(lite)[[1]])
  for(i in seq_len(dim(lite)[[1]])){
    if(sum(grepl('msi_ml', colnames(lite))) > 0){
      if(sum(colnames(lite) == "MMRD_category") > 0){
        if(lite$MMRD_category[[i]] %in% c('MMRD', 'POLE'))
          categs[[i]] <- lite$MMRD_category[[i]]
      }
      else{
        if(lite$Signature_msi_ml[[i]] >= 0.99 & lite$total_snvs[[i]] > 10)
          categs[[i]] <- 'Signature_msi'
        else if(lite$Signature_pole_ml[[i]] > 0.9999)
          categs[[i]] <- 'Signature_pole'
      }
    }
    if(!(categs[[i]] %in% c('MMRD', 'POLE', 'Signature_msi', 'Signature_pole'))){
      if(sum(colnames(lite) %in% c('pass_mva')) > 0){ 
        if(lite$pass_mva_strict[[i]])
          categs[[i]] <- 'Signature_3_hc'
        else if(lite$pass_mva[[i]])
          categs[[i]] <- 'Signature_3_lc'
        else{
          cols_ml_this <- cols_ml[!(cols_ml %in% c('Signature_3_ml', 'Signature_msi_ml', 'Signature_pole_ml'))]
          ind_max <- which(max(lite[i, cols_ml_this]) == lite[i,])
          categs[[i]] <- gsub(paste0(colnames(lite)[ind_max], collapse = ':'),
                              pattern = '_ml',
                              replace = '')
        }
      }
      else{
        cols_ml_this <- cols_ml[!(cols_ml %in% c('Signature_msi_ml', 'Signature_pole_ml'))]
        ind_max <- which(max(lite[i, cols_ml_this]) == lite[i,])
        categs[[i]] <- gsub(paste0(colnames(lite)[ind_max], collapse = ':'),
                            pattern = '_ml',
                            replace = '')
      }
    }
  }
  lite$categ <- categs
  
  if(length(grep('msi_ml', colnames(lite))) > 0){
    inds <- which(lite$categ %in% c("Signature_msi", "Signature_pole", "MMRD", "POLE"))
    merged_output$exps_all_msi <- as.character(merged_output$exps_all_msi)
    merged_output$sigs_all_msi <- as.character(merged_output$sigs_all_msi)
    lite$exps_all[inds] <- merged_output$exps_all_msi[inds]
    lite$sigs_all[inds] <- merged_output$sigs_all_msi[inds]
  }
  
  return(lite)
}