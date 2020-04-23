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
  cols_apobec_ml <- which(grepl('APOBEC', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sig3_ml <- which(grepl('Signature_3', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sig8_ml <- which(grepl('Signature_8', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sig9_ml <- which(grepl('Signature_9', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sig17_ml <- which(grepl('Signature_17', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sig18_ml <- which(grepl('Signature_18', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sig28_ml <- which(grepl('Signature_28', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_sigN1_ml <- which(grepl('Signature_N1', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_smoke_ml <- which(grepl('Signature_smoke', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_clock_ml <- which(grepl('Signature_clock', all_cols) & grepl('_ml', all_cols) & !grepl('_msi', all_cols))
  cols_pole_ml <- which(grepl('Signature_pole', all_cols) & grepl('_ml', all_cols))
  cols_msi_ml <- which(grepl('Signature_msi', all_cols) & grepl('_ml', all_cols))

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

 
  ncols <- dim(lite)[[2]]
  if(length(cols_apobec_ml) > 0)
    lite$Signature_APOBEC_ml <- sum_cols(cols_apobec_ml)
  if(length(cols_sig3_ml) > 0)
    lite$Signature_3_ml <- sum_cols(cols_sig3_ml)
  if(length(cols_sig8_ml) > 0)
    lite$Signature_8_ml <- sum_cols(cols_sig8_ml)
  if(length(cols_sig9_ml) > 0)
    lite$Signature_9_ml <- sum_cols(cols_sig9_ml)
  if(length(cols_sig17_ml) > 0)
    lite$Signature_17_ml <- sum_cols(cols_sig17_ml)
  if(length(cols_sig18_ml) > 0)
    lite$Signature_18_ml <- sum_cols(cols_sig18_ml)
  if(length(cols_sig28_ml) > 0)
    lite$Signature_28_ml <- sum_cols(cols_sig28_ml)
  if(length(cols_sigN1_ml) > 0)
    lite$Signature_N1_ml <- sum_cols(cols_sigN1_ml)
  if(length(cols_smoke_ml) > 0)
    lite$Signature_smoke_ml <- sum_cols(cols_smoke_ml)
  if(length(cols_pole_ml) > 0)
    lite$Signature_pole_ml <- sum_cols(cols_pole_ml)
  if(length(cols_msi_ml) > 0)
    lite$Signature_msi_ml <- sum_cols(cols_msi_ml)
  if(length(cols_clock_ml) > 0)
    lite$Signature_clock_ml <- sum_cols(cols_clock_ml)


  if(ncols != dim(lite)[[2]]){
    cols_ml <- colnames(lite)[(ncols + 1):(dim(lite)[[2]])]
  }

  if(length(inds_keep) > 0)
    lite <- cbind(lite, merged_output[, inds_keep])

  if(sum(colnames(lite) == "pass_mva") > 0){
    categs <- rep('', dim(lite)[[1]])
    if(exists('cols_ml')){
      if(sum(grepl('msi_ml', colnames(lite))) > 0){
        for(i in seq_len(dim(lite)[[1]])){
          if(lite$Signature_msi_ml[[i]] >= 0.95 & lite$total_snvs[[i]] > 10)
            categs[[i]] <- 'Signature_msi'
          else if(lite$Signature_pole_ml[[i]] > 0.9999)
            categs[[i]] <- 'Signature_pole'
          else if(lite$pass_mva_strict[[i]])
            categs[[i]] <- 'Signature_3_hc'
          else if(lite$pass_mva[[i]])
            categs[[i]] <- 'Signature_3_lc'
          else{
            cols_ml_this <- cols_ml[cols_ml != 'Signature_3_ml' & cols_ml != 'Signature_msi_ml' & cols_ml != 'Signature_pole_ml']
            ind_max <- which(max(lite[i, cols_ml_this]) == lite[i,])
            categs[[i]] <- gsub(paste0(colnames(lite)[ind_max], collapse = ':'),
                              pattern = '_ml',
                              replace = '')
        
          }
        }
     }
     else{
       for(i in seq_len(dim(lite)[[1]])){
          if(lite$pass_mva_strict[[i]])
            categs[[i]] <- 'Signature_3_hc'
          else if(lite$pass_mva[[i]])
            categs[[i]] <- 'Signature_3_lc'
          else{
            cols_ml_this <- cols_ml[cols_ml != 'Signature_3_ml' & cols_ml != 'Signature_msi_ml' & cols_ml != 'Signature_pole_ml']
            ind_max <- which(max(lite[i, cols_ml_this]) == lite[i,])
            categs[[i]] <- gsub(paste0(colnames(lite)[ind_max], collapse = ':'),
                              pattern = '_ml',
                              replace = '')
        
          }
        }
      }
      lite$categ <- categs
    }
  }
  else if(sum(colnames(lite) == "pass_ml") > 0){
    categs <- rep('', dim(lite)[[1]])
    if(exists('cols_ml')){
      if(sum(grepl('msi_ml', colnames(lite))) > 0){
        for(i in seq_len(dim(lite)[[1]])){
          if(lite$Signature_msi_ml[[i]] >= 0.95 & lite$total_snvs[[i]] > 10)
            categs[[i]] <- 'Signature_msi'
          else if(lite$Signature_pole_ml[[i]] > 0.9999)
            categs[[i]] <- 'Signature_pole'
          else if(lite$pass_ml[[i]]) 
            categs[[i]] <- 'Signature_3 (No MVA)'
          else{
            cols_ml_this <- cols_ml[cols_ml != 'Signature_3_ml' & cols_ml != 'Signature_msi_ml' & cols_ml != 'Signature_pole_ml']
            ind_max <- which(max(lite[i, cols_ml_this]) == lite[i,])
            categs[[i]] <- gsub(paste0(colnames(lite)[ind_max], collapse = ':'),
                              pattern = '_ml',
                              replace = '')
        
          }
        }
     }
     else{
       for(i in seq_len(dim(lite)[[1]])){
          if(lite$pass_ml[[i]])
            categs[[i]] <- 'Signature_3 (No MVA)'
          else{
            cols_ml_this <- cols_ml[cols_ml != 'Signature_3_ml' & cols_ml != 'Signature_msi_ml' & cols_ml != 'Signature_pole_ml']
            ind_max <- which(max(lite[i, cols_ml_this]) == lite[i,])
            categs[[i]] <- gsub(paste0(colnames(lite)[ind_max], collapse = ':'),
                              pattern = '_ml',
                              replace = '')
        
          }
        }
      }
      lite$categ <- categs
    }
  }
  
  if(length(grep('msi_ml', colnames(lite))) > 0){
    inds <- which(lite$categ == "Signature_msi" | lite$categ == "Signature_pole")
    merged_output$exps_all_msi <- as.character(merged_output$exps_all_msi)
    merged_output$sigs_all_msi <- as.character(merged_output$sigs_all_msi)
    lite$exps_all[inds] <- merged_output$exps_all_msi[inds]
    lite$sigs_all[inds] <- merged_output$sigs_all_msi[inds]
  }
  
  return(lite)
}