#' Generates a detailed plot per sample 
#'
#' @param file the csv file produced by SigMA
#' @param sample name to be plotted

plot_detailed <- function(file = NULL, sample = NULL){
  text_size = 10

  df <- read.csv(file)

  df <- df[, -grep('Signature_N', colnames(df))]
  df <- df[, -grep('Signature_17v2', colnames(df))]

  tumor_type <- unlist(strsplit(unlist(strsplit(file,
                       split = 'tumortype_'))[[2]],
                       split = '_platform_'))[[1]]

  platform <- unlist(strsplit(unlist(strsplit(file,
                     split = 'platform_'))[[2]],
                     split = '.csv'))[[1]]
  
  signames <- signames_per_tissue[[tumor_type]]

  this <- df[df$tumor == sample,]
###cosine simil################################################################################################## 
  # strip cosine similarity values
  inds <- grep('_c', colnames(this))
  inds_rm <- grep('_c_diff|_ml', colnames(this))
  cos_vals <- this[, inds[-na.omit(match(inds_rm, inds))]]
  cos_vals <- cos_vals[, grep('Signature', colnames(cos_vals))]

  # tumor type specific signatures vs others
  cos_vals_tt <- cos_vals[, paste0(signames, '_c')]
  cos_vals_tt <- cos_vals_tt[, order(-cos_vals_tt[1,])]
  nsig_tt <- dim(cos_vals_tt)[[2]]

  cos_vals_ntt <- cos_vals[, -match(paste0(signames, '_c'), colnames(cos_vals))]
  cos_vals_ntt <- cos_vals_ntt[, order(-cos_vals_ntt[1,])]

  cos_vals <- cbind(cos_vals_tt, cos_vals_ntt)  
  signame_cos <- gsub(gsub(colnames(cos_vals),
                           pattern = '_c',
                           replace = ''),
                      pattern = '_', replace = '')
  # plot
  df_cos <- data.frame(cosine = unlist(cos_vals[1,]), 
                       signame = signame_cos,
                       color = c(rep('tumor-type specific', nsig_tt),
                                 rep('others', dim(cos_vals)[[2]] - nsig_tt)))
  df_cos$color <- as.character(df_cos$color)
  df_cos$color[df_cos$signame == "Signature3"] <- 'sig3'

  df_cos <- transform(df_cos, color = factor(color, levels = c('sig3', 'others', 'tumor-type specific')))
  df_cos <- transform(df_cos, signame = factor(signame, levels = signame_cos))

  plot_cos <- ggplot2::ggplot(df_cos, ggplot2::aes(x = signame, y = cosine)) 
  plot_cos <- plot_cos + ggplot2::geom_bar(stat = 'identity', ggplot2::aes(fill = color))
  plot_cos <- plot_cos + ggplot2::scale_fill_manual(values = as.character(c('#f9bb68', col_pos_neg)))
  plot_cos <- plot_cos + theme_def + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                           axis.title.y = ggplot2::element_text(size = text_size),
                                           axis.text.y = ggplot2::element_text(size = text_size),
                                           axis.text.x = ggplot2::element_text(size = text_size, 
                                                                      hjust = 1, 
                                                                      angle = 90))
  plot_cos <- plot_cos + ggplot2::ylab('Cosine similarity') 

###likelihood##################################################################################################
  # columns related to likelihood
  if(length(grep('ml_msi', colnames(this))) > 0){
    inds_msi <- grep('_msi_c1_ml_msi|_msi_c2_ml_msi|_msi_c3_ml_msi|_msi_c4_ml_msi', colnames(this))
    inds_pole <- grep('_pole_c1_ml_msi', colnames(this))
    inds_rm <- grep('_ml_msi', colnames(this))
    inds <- grep('_ml', colnames(this)[-inds_rm])
    ml_vals <- this[, c(inds, inds_msi, inds_pole)]
  }else{
    inds <- grep('_ml', colnames(this))
    ml_vals <- this[, inds]
  }

  ml_vals <- ml_vals[, grep('Signature', colnames(ml_vals))]

  # determine the groups
  groups <-  unique(c(unlist(lapply(strsplit(colnames(ml_vals),
                             split = paste0('_c', 1:10, '_ml', collapse = '|')), 
                      function(x){ x[[1]] }))))
  
  is_sig3 <- (groups == "Signature_3")

  # sum the likelihood of clusters within each group
  vals <- numeric(length(groups))
  for(i in 1:length(groups)){
    vals[[i]] <- sum(ml_vals[, grep(groups[[i]], colnames(ml_vals))])
  }

  groups <- gsub(groups, pattern = 'Signature_msi', replace = 'MSI')
  groups <- gsub(groups, pattern = 'Signature_pole', replace = 'POLE')
  groups <- gsub(groups, pattern = 'Signature_clock', replace = 'Clock')
  groups <- gsub(groups, pattern = 'Signature_APOBEC', replace = 'APOBEC')
  groups <- gsub(groups, pattern = 'Signature_', replace = 'Signature')

  df_ml <- data.frame(likelihood = vals,
                      group = groups, 
                      is_sig3)
  df_ml <- df_ml[order(-df_ml$likelihood), ]
  df_ml <- transform(df_ml, is_sig3 = factor(is_sig3, levels = c(T, F)))
  df_ml <- transform(df_ml, group = factor(group, levels = as.character(df_ml$group)))

  # plot likelihood
  plot_ml <- ggplot2::ggplot(df_ml, ggplot2::aes(x = group, y = likelihood))
  plot_ml <- plot_ml + ggplot2::geom_bar(stat = 'identity', ggplot2::aes(fill = is_sig3)) 
  plot_ml <- plot_ml + theme_def + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                         axis.title.y = ggplot2::element_text(size = text_size),
                                         axis.text.y = ggplot2::element_text(size = text_size),
                                         axis.text.x = ggplot2::element_text(size = text_size,
                                                                      hjust = 1,
                                                                      angle = 90))
  plot_ml <- plot_ml + ggplot2::scale_fill_manual(values = c('#f9bb68', col_pos_neg[[2]]))
  plot_ml <- plot_ml + ggplot2::ylab('Likelihood') 

###exposures################################################################################################## 
  exps <- as.numeric(unlist(strsplit(as.character(this$exps_all), split = '_')))
  sigs <- unlist(strsplit(as.character(this$sigs_all), split = '\\.'))

  if(length(grep('ml_msi', colnames(this))) > 0){
    this$Signature_msi_ml <- rowSums(this[, paste0('Signature_msi_c', 1:4, '_ml_msi')])
    if(this$Signature_msi_ml > 0.99 | this$Signature_pole_c1_ml_msi > 0.99){
      exps <- as.numeric(unlist(strsplit(as.character(this$exps_all_msi), split = '_')))
      sigs <- unlist(strsplit(as.character(this$sigs_all_msi), split = '\\.'))
    }
  }
  sigs <- gsub(sigs, pattern = '_', replace = '')

  df_exp <- data.frame(exps = exps, sigs = sigs, is_sig3 = (sigs == "Signature3"))
  df_exp <- df_exp[order(-df_exp$exps), ]

  df_exp <- transform(df_exp, is_sig3 = factor(is_sig3, levels = c(T, F)))
  df_exp <- transform(df_exp, sigs = factor(sigs, levels = as.character(sigs)))

  plot_exp <- ggplot2::ggplot(df_exp, ggplot2::aes(x = sigs, y = exps)) 
  plot_exp <- plot_exp + ggplot2::geom_bar(stat = 'identity', ggplot2::aes(fill = is_sig3))
  plot_exp <- plot_exp + theme_def + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                           axis.title.y = ggplot2::element_text(size = text_size),
                                           axis.text.y = ggplot2::element_text(size = text_size),
                                           axis.text.x = ggplot2::element_text(size = text_size,
                                                                      hjust = 1,
                                                                      angle = 90))
  if(sum(sigs == "Signature3") > 0) 
    plot_exp <- plot_exp + ggplot2::scale_fill_manual(values = c('#f9bb68', col_pos_neg[[2]]))
  else
    plot_exp <- plot_exp + ggplot2::scale_fill_manual(values = c(col_pos_neg[[2]]))
  
  plot_exp <- plot_exp + ggplot2::ylab('NNLS exposure') 


  tribase <- plot_tribase_dist(as.data.frame(t(this[1, 1:96])))

  lay <- rbind(c(1, 1, 1, 1, 1),
               c(3,4, 2, 2, 2))

  plot_arr <- gridExtra::grid.arrange(tribase, plot_cos, plot_ml, plot_exp,
                           layout_matrix = lay)
  return(plot_arr)
  
}