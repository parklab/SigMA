#' Generates summary plot 
#'
#' @param file the csv file produced by SigMA

plot_summary <- function(file = NULL){
  df <- read.csv(file)
  
  # read SigMA settings from file name
  msk_string <- gsub(platform_names[['msk']],
                     pattern = ' ', replace = '')

  exome_string <- gsub(platform_names[['seqcap']],
                     pattern = ' ', replace = '')

  wgs_string <- gsub(platform_names[['wgs']],
                     pattern = ' ', replace = '')

  tumor_type <- unlist(strsplit(unlist(strsplit(file, 
                                                split = 'tumortype_'))[[2]],
                       split = '_platform_'))[[1]]

  platform <- unlist(strsplit(unlist(strsplit(file, 
                                              split = 'platform_'))[[2]],
                     split = '_cf'))[[1]]

  if(platform == msk_string) data = 'msk'
  if(platform == exome_string) data = 'seqcap'
  if(platform == wgs_string) data = 'wgs'

  if(data == 'msk'){
    cutoff_strict <- cutoffs_msk_strict[[tumor_type]]
    cutoff <- cutoffs_msk[[tumor_type]]
  }else if(data == 'seqcap'){
    cutoff_strict <- cutoffs_exome_strict[[tumor_type]]
    cutoff <- cutoffs_exome[[tumor_type]]
  }else if(data == 'wgs'){
    cutoff_strict <- cutoffs_wgs_strict[[tumor_type]]
    cutoff <- cutoffs_wgs[[tumor_type]]
  }

  text_size = 10

  # 3-base distributions
  if(sum(df$pass_mva) > 0){ 
    inds_pos <- which(df$pass_mva)
    pos_tribase <- plot_tribase_dist(as.data.frame(colSums(df[inds_pos, 1:96])),
                                     signame = paste0('  Aggregated mutations from ', 
                                                      sum(df$pass_mva), ' Sig3+ sample(s)'))
  }
  if(sum(!df$pass_mva) > 0){
    inds_neg <- which(!df$pass_mva)
    neg_tribase <- plot_tribase_dist(as.data.frame(colSums(df[inds_neg, 1:96])),
                                     signame = paste0('  Aggregated mutations from ', 
                                                      sum(df$pass_mva == F), ' Sig3- sample(s)'))
    neg_tribase <- neg_tribase + ggplot2::theme(legend.position = "none")
  }

  # cosine density distribution
  colors_cos <- c(col_pos_neg[[2]], col_pos_neg[[1]])

  if(length(unique(df$pass_mva)) == 1){
    if(unique(df$pass_mva))
      colors_cos <- col_pos_neg[[2]]
    else
      colors_cos <- col_pos_neg[[1]]
  }

  
  cos_bin_width = 0.1
  cos_bins <- seq(0, 1, by = cos_bin_width)
  counts_cos_pos <- rep(0, length(cos_bins) - 1)
  counts_cos_neg <- rep(0, length(cos_bins) - 1)
  bin_center <- rep(0, length(cos_bins) - 1)
  
  for(ibin in seq_len(length(cos_bins) - 1)){
    cos_vals_pos <- df$Signature_3_c[df$pass_mva]
    cos_vals_neg <- df$Signature_3_c[!df$pass_mva]
    counts_cos_pos[[ibin]] <- sum(cos_vals_pos >= cos_bins[[ibin]] & cos_vals_pos < cos_bins[[ibin + 1]]) 
    counts_cos_neg[[ibin]] <- sum(cos_vals_neg >= cos_bins[[ibin]] & cos_vals_neg < cos_bins[[ibin + 1]]) 

    bin_center[[ibin]] <- (cos_bins[[ibin]] + cos_bins[[ibin + 1]])/2
  }

  
  if(sum(counts_cos_pos) > 0){
    df_cos <- data.frame(cos = bin_center, count = counts_cos_pos, group = 'Sig3+')
    if(sum(counts_cos_neg) > 0){
      df_cos <- rbind(df_cos, data.frame( cos = bin_center, count = counts_cos_neg, group = 'Sig3-' ))
    }
  }
  else if(sum(counts_cos_neg) > 0){
    df_cos <- data.frame( cos = bin_center, count = counts_cos_neg, group = 'Sig3-' )
  }
  else{
    stop('no samples')
  }

  plot_cos <- ggplot2::ggplot(df_cos, ggplot2::aes(x = cos, y = count))
  plot_cos <- plot_cos + ggplot2::geom_bar(stat = 'identity', position = 'dodge', 
                                           ggplot2::aes(fill = group))
  plot_cos <- plot_cos + ggplot2::scale_fill_manual(values = colors_cos)
  plot_cos <- plot_cos + theme_def + ggplot2::theme(axis.title.x = ggplot2::element_text(size = text_size),
                                                    axis.title.y = ggplot2::element_text(size = text_size),
                                                    axis.text.x = ggplot2::element_text(size = text_size),
                                                    axis.text.y = ggplot2::element_text(size = text_size))

  plot_cos <- plot_cos + ggplot2::xlab('Cos simil of Sig3')
  plot_cos <- plot_cos + ggplot2::xlim(0, 1)

  # likelihood density distribution
  inds <- grep('_ml', colnames(df))
  if(length(grep('_ml_msi', colnames(df))) > 0){
    inds_rm <- grep('_ml_msi', colnames(df))
    inds <- inds[-na.omit(match(inds_rm, inds))]
  }

  ind_tumor <- grep('tumor', colnames(df))
  inds <- c(inds, ind_tumor)
  rm(ind_tumor)

  df_ml <- df[, inds]
  df_ml <- df_ml[, grep('Signature_3', colnames(df_ml))]
  df_ml$Signature_3_ml <- rowSums(df_ml)
  df_ml$tumor <- as.character(df$tumor)
  df_ml$pass_mva <- df$pass_mva
  df_ml$Signature_3_mva <- df$Signature_3_mva 

  ml_bin_width = 0.1
  ml_bins <- seq(0, 1.001, by = ml_bin_width)
  counts_ml_pos <- rep(0, length(ml_bins) - 1)
  counts_ml_neg <- rep(0, length(ml_bins) - 1)
  bin_center <- rep(0, length(ml_bins) - 1)
  
  ml_vals_pos <- df_ml$Signature_3_ml[df_ml$pass_mva]
  ml_vals_neg <- df_ml$Signature_3_ml[!df_ml$pass_mva]

  
  for(ibin in seq_len(length(ml_bins) - 1)){
    counts_ml_pos[[ibin]] <- sum(ml_vals_pos >= ml_bins[[ibin]] & ml_vals_pos < ml_bins[[ibin + 1]]) 
    counts_ml_neg[[ibin]] <- sum(ml_vals_neg >= ml_bins[[ibin]] & ml_vals_neg < ml_bins[[ibin + 1]]) 

    bin_center[[ibin]] <- (ml_bins[[ibin]] + ml_bins[[ibin + 1]])/2
  }

  if(sum(counts_ml_pos) > 0){
    df_ml_binned <- data.frame(ml = bin_center, count = counts_ml_pos, group = 'Sig3+')
    if(sum(counts_ml_neg) > 0){
      df_ml_binned <- rbind(df_ml_binned, data.frame( ml = bin_center, count = counts_ml_neg, group = 'Sig3-' ))
    }
  }
  else if(sum(counts_ml_neg) > 0){
    df_ml_binned <- data.frame( ml = bin_center, count = counts_ml_neg, group = 'Sig3-' )
  }
  else{
    stop('no samples')
  }

  plot_ml <- ggplot2::ggplot(df_ml_binned, ggplot2::aes(x = ml, y = count))
  plot_ml <- plot_ml + ggplot2::geom_bar(stat = 'identity', position = 'dodge', 
                                           ggplot2::aes(fill = group))
  plot_ml <- plot_ml + ggplot2::scale_fill_manual(values = colors_cos) 
  # here add a line to strip away the legend
  plot_ml <- plot_ml + theme_def + ggplot2::theme(axis.title.x = ggplot2::element_text(size = text_size),
                                                  axis.title.y = ggplot2::element_text(size = text_size),
                                                  axis.text.x = ggplot2::element_text(size = text_size),
                                                  axis.text.y = ggplot2::element_text(size = text_size))


  plot_ml <- plot_ml + ggplot2::xlab('Likelihood of Sig3')
  plot_ml <- plot_ml + ggplot2::xlim(0, 1)  

  # score bar plot
  df_ml$Signature_3_mva <- df_ml$Signature_3_mva 
  df_ml$group <- 'Sig3 -'
  df_ml$group[df_ml$Signature_3_mva > cutoff] <- 'lc Sig3 +'
  df_ml$group[df_ml$Signature_3_mva > cutoff_strict] <- 'Sig3 +'

  df_color <- data.frame(color = c(col_pos_neg[[1]], "#dee9fc", col_pos_neg[[2]]),
                         group = c('Sig3 -', 'lc Sig3 +', 'Sig3 +'))
  size_nchar <- mean(unlist(apply(df_ml, 1, function(x){ nchar(x['tumor']) })))
  df_ml$tumor_label <- df_ml$tumor
  if(size_nchar > 20){
    df_ml$tumor_label <- paste0('tumor', 1:dim(df_ml)[[1]])
  }
  
  
  df_ml$group_fac <- df_ml$group
  df_ml <- transform(df_ml, group_fac = factor(group_fac, 
                                               levels = unique(df_ml$group)))
  df_ml <- transform(df_ml, tumor_label = factor(tumor_label, levels = df_ml$tumor_label[order(df_ml$Signature_3_mva)]))

  plot_score <- ggplot2::ggplot(df_ml, 
                                ggplot2::aes(x = tumor_label, 
                                             y = Signature_3_mva))
  plot_score <- plot_score + ggplot2::geom_bar(ggplot2::aes(fill = group_fac),
                                               stat = 'identity')
  plot_score <- plot_score + ggplot2::scale_fill_manual(values = as.character(df_color$color[match(unique(df_ml$group), 
                                                                                             df_color$group)]))
  plot_score <- plot_score + theme_def 
  plot_score <- plot_score + ggplot2::theme(panel.border = ggplot2::element_blank(),
                                            axis.line.y = ggplot2::element_line(size = 0.3, 
                                                                                linetype = "solid", 
                                                                                colour = "black"),
                                            axis.line.x = ggplot2::element_line(size = 0.3, 
                                                                                linetype = "solid", 
                                                                                colour = "black"),
                                            axis.text.x = ggplot2::element_text(size = text_size,
                                                                                hjust = 1,
                                                                                angle = 45), 
                                            axis.title.y = ggplot2::element_text(size = text_size),
                                            axis.text.y = ggplot2::element_text(size = text_size),
                                            legend.position = 'top')
  plot_score <- plot_score + ggplot2::ylab('MVA score') + ggplot2::xlab('')
  plot_score <- plot_score + ggplot2::geom_hline(ggplot2::aes(yintercept = cutoff))
  plot_score <- plot_score + ggplot2::geom_hline(ggplot2::aes(yintercept = cutoff_strict), 
                                                 linetype = 'dashed')

  # Signature 3 exposures 
  colors_exp <- colors_cos
  if(length(unique(df$pass_mva[df$exp_sig3 > 0])) == 1){
    if(unique(df$pass_mva[df$exp_sig3 > 0]))
      colors_exp <- c(colors_cos[[2]], colors_cos[[1]])
  }

  exp_bin_width <- (max(df$exp_sig3) + 1)/10
  exp_bins <- seq(0.001, max(df$exp_sig3) + 1, by = exp_bin_width)
  counts_exp_pos <- rep(0, length(exp_bins) - 1)
  counts_exp_neg <- rep(0, length(exp_bins) - 1)
  bin_center <- rep(0, length(exp_bins) - 1)
  
  for(ibin in seq_len(length(exp_bins) - 1)){
    exp_vals_pos <- df$exp_sig3[df$pass_mva]
    exp_vals_neg <- df$exp_sig3[!df$pass_mva]
    counts_exp_pos[[ibin]] <- sum(exp_vals_pos >= exp_bins[[ibin]] & exp_vals_pos < exp_bins[[ibin + 1]]) 
    counts_exp_neg[[ibin]] <- sum(exp_vals_neg >= exp_bins[[ibin]] & exp_vals_neg < exp_bins[[ibin + 1]]) 

    bin_center[[ibin]] <- (exp_bins[[ibin]] + exp_bins[[ibin + 1]])/2
  }

  if(sum(counts_exp_pos) > 0){
    df_exp <- data.frame(exp = bin_center, count = counts_exp_pos, group = 'Sig3+')
    if(sum(counts_exp_neg) > 0){
      df_exp <- rbind(df_exp, data.frame( exp = bin_center, count = counts_exp_neg, group = 'Sig3-' ))
    }
  }
  else if(sum(counts_exp_neg) > 0){
    df_exp <- data.frame( exp = bin_center, count = counts_exp_neg, group = 'Sig3-' )
  }

  if(sum(counts_exp_pos) > 0 | sum(counts_exp_neg) > 0){
    plot_exp <- ggplot2::ggplot(df_exp, ggplot2::aes(x = exp, y = count))
    plot_exp <- plot_exp + ggplot2::geom_bar(stat = 'identity', position = 'dodge',
                                             ggplot2::aes(fill = group))
    plot_exp <- plot_exp + ggplot2::scale_fill_manual(values = colors_exp)
    plot_exp <- plot_exp + theme_def + ggplot2::theme(axis.title.x = ggplot2::element_text(size = text_size),
                                                      axis.title.y = ggplot2::element_text(size = text_size),
                                                      axis.text.x = ggplot2::element_text(size = text_size),
                                                      axis.text.y = ggplot2::element_text(size = text_size))
    plot_exp <- plot_exp + ggplot2::xlab('Sig3 exposure')
    plot_exp <- plot_exp + ggplot2::xlim(0, max(df$exp_sig3))
  }

  # combined plot
  pdf('summary_short_sig3.pdf', width = 1.2*9.8/2.4, 7)
  if(sum(df$pass_mva) > 0 & sum(!df$pass_mva > 0)){
    if(sum(counts_exp_pos) > 0 | sum(counts_exp_neg) > 0){
      lay <- rbind(c(6, 6, 6),
                   c(3, 4, 5),
                   c(1, 1, 1),
                   c(2, 2, 2))
   
      plot_arr <- gridExtra::grid.arrange(pos_tribase,
                                     neg_tribase,
                                     plot_cos,
                                     plot_ml,
                                     plot_exp,
                                     plot_score,
                                     layout_matrix = lay, 
                                     heights = c(1.5, 0.8, 1.1, 0.85))
    }else{
      lay <- rbind(c(5, 5),
                   c(3, 4),
                   c(1, 1),
                   c(2, 2))
   
      plot_arr <- gridExtra::grid.arrange(pos_tribase,
                                     neg_tribase,
                                     plot_cos,
                                     plot_ml,
                                     plot_score,
                                     layout_matrix = lay, 
                                     heights = c(1.5, 0.8, 1.1, 0.85))
    }
  }else if(sum(df$pass_mva) > 0 & sum(!df$pass_mva == 0)){
    if(sum(counts_exp_pos) > 0 | sum(counts_exp_neg) > 0){
      lay <- rbind(c(5, 5, 5),
                   c(2, 3, 4),
                   c(1, 1, 1))
      plot_arr <- gridExtra::grid.arrange(pos_tribase,
                                     plot_cos,
                                     plot_ml,
                                     plot_exp,
                                     plot_score,
                                     layout_matrix = lay, 
                                     heights = c(1.5, 0.8, 1.1, 0.85))
    }else{
      lay <- rbind(c(4, 4),
                   c(2, 3),
                   c(1, 1))
      plot_arr <- gridExtra::grid.arrange(pos_tribase,
                                     plot_cos,
                                     plot_ml,
                                     plot_score,
                                     layout_matrix = lay, 
                                     heights = c(1.5, 0.8, 1.1, 0.85))
    }
  }else if(sum(df$pass_mva) == 0 & sum(!df$pass_mva > 0)){
    if(sum(counts_exp_pos) > 0 | sum(counts_exp_neg) > 0){
      lay <- rbind(c(5, 5, 5),
                   c(2, 3, 4),
                   c(1, 1, 1))
      plot_arr <- gridExtra::grid.arrange(neg_tribase,
                                     plot_cos,
                                     plot_ml,
                                     plot_exp,
                                     plot_score,
                                     layout_matrix = lay, 
                                     heights = c(1.5, 0.8, 1.1, 0.85))
    }
    else{
      lay <- rbind(c(4, 4),
                   c(2, 3),
                   c(1, 1))
      plot_arr <- gridExtra::grid.arrange(neg_tribase,
                                     plot_cos,
                                     plot_ml,
                                     plot_score,
                                     layout_matrix = lay, 
                                     heights = c(1.5, 0.8, 1.1, 0.85))
    }
  }

  dev.off()                      
  return(plot_arr)
}