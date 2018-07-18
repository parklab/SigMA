#' Generates summary plot 
#'
#' @param file the csv file produced by SigMA
#' @param detailed boolean determines whether a detailed 
#' summary plot or a short summary should be produced, if
#' set to true creates a plot for each sample, if false 
#' makes a single plot telling what fraction of samples had
#' were determined to be Signature 3-positive
#'

plot_summary <- function(file = NULL, detailed = F){
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
                       split = '_'))[[1]]

  platform <- unlist(strsplit(unlist(strsplit(file, 
                                              split = 'platform_'))[[2]],
                     split = '.csv'))[[1]]


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

  text_size = 7

  # 3-base distributions 
  pos_tribase <- plot_tribase_dist(as.data.frame(t(df[which(df$pass_mva), 1:96])))
  neg_tribase <- plot_tribase_dist(as.data.frame(t(df[-which(df$pass_mva), 1:96])))
  neg_tribase <- neg_tribase + ggplot2::theme(legend.position = "none")

  # cosine density distribution
  plot_cos <- ggplot2::ggplot(df, ggplot2::aes(x = Signature_3_c))
  plot_cos <- plot_cos + ggplot2::geom_density(ggplot2::aes(color = pass_mva))
  plot_cos <- plot_cos + ggplot2::scale_color_manual(values = col_pos_neg)
  plot_cos <- plot_cos + theme_def
  plot_cos <- plot_cos + ggplot2::xlab('cosine Sig3')
  plot_cos <- plot_cos + ggplot2::xlim(0, 1)

  # likelihood density distribution
  inds <- grep('_ml', colnames(df))
  ind_tumor <- grep('tumor', colnames(df))
  inds <- c(inds, ind_tumor)
  rm(ind_tumor)

  df_ml <- df[, inds]
  df_ml <- df_ml[, grep('Signature_3', colnames(df_ml))]
  df_ml$Signature_3_ml <- rowSums(df_ml)
  df_ml$tumor <- as.character(df$tumor)
  df_ml$pass_mva <- df$pass_mva
  df_ml$Signature_3_mva <- df$Signature_3_mva 

  plot_ml <- ggplot2::ggplot(df_ml, ggplot2::aes(x = Signature_3_ml))
  plot_ml <- plot_ml + ggplot2::geom_density(ggplot2::aes(color = pass_mva))
  plot_ml <- plot_ml + ggplot2::scale_color_manual(values = col_pos_neg) 
  # here add a line to strip away the legend
  plot_ml <- plot_ml + theme_def
  plot_ml <- plot_ml + ggplot2::xlab('likelihood Sig3')
  plot_ml <- plot_ml + ggplot2::xlim(0, 1)  

  # score bar plot
  df_ml$Signature_3_mva <- df_ml$Signature_3_mva 
  df_ml$group <- 'Sig3 -'
  df_ml$group[df_ml$Signature_3_mva > cutoff] <- 'lc Sig3 +'
  df_ml$group[df_ml$Signature_3_mva > cutoff_strict] <- 'Sig3 +'
  df_color <- data.frame(color = c(col_pos_neg[[1]], "#ffdf12", col_pos_neg[[2]]),
                         group = c('Sig3 -', 'lc Sig3 +', 'Sig3 +'))
  size_nchar <- mean(unlist(apply(df_ml, 1, function(x){ nchar(x['tumor']) })))
  df_ml$tumor_label <- df_ml$tumor
  if(size_nchar > 20){
    df_ml$tumor_label <- paste0('tumor', 1:dim(df_ml)[[1]])
  }
  df_ml$group_fac <- df_ml$group
  df_ml <- transform(df_ml, group_fac = factor(group_fac, 
                                               levels = unique(df_ml$group)))


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
                                            legend.position = 'top')
  plot_score <- plot_score + ggplot2::ylab('MVA score') + ggplot2::xlab('')
  plot_score <- plot_score + ggplot2::geom_hline(ggplot2::aes(yintercept = cutoff))
  plot_score <- plot_score + ggplot2::geom_hline(ggplot2::aes(yintercept = cutoff_strict), 
                                                 linetype = 'dashed')

  # Signature 3 exposures 
  plot_exp <- ggplot2::ggplot(df[df$exp_sig3 > 0,], ggplot2::aes(x = exp_sig3))
  plot_exp <- plot_exp + ggplot2::geom_density(ggplot2::aes(color = pass_mva))
  plot_exp <- plot_exp + ggplot2::scale_color_manual(values = col_pos_neg)
  plot_exp <- plot_exp + theme_def
  plot_exp <- plot_exp + ggplot2::xlab('exp Sig3')
  plot_exp <- plot_exp + ggplot2::xlim(0, max(df$exp_sig3))

  # combined plot
  lay <- rbind(c(6, 6, 6),
               c(3, 4, 5),
               c(1, 1, 1),
               c(2, 2, 2))

  pdf('summary_short_sig3.pdf', width = 1.2*9.8/2.4, 7)
  plot_arr <- gridExtra::grid.arrange(pos_tribase,
                                 neg_tribase,
                                 plot_cos,
                                 plot_ml,
                                 plot_exp,
                                 plot_score,
                                 layout_matrix = lay, 
                                 heights = c(1.3, 0.7, 1.4, 1))

  grid::grid.text(paste0(sum(df$pass_mva), " Sig3+ sample(s)"),
            x = grid::unit(0.1, "npc"), y = grid::unit(0.47, "npc"),
            gp = grid::gpar(fontsize = text_size))

  grid::grid.text(paste0(sum(df$pass_mva == F), " Sig3- sample(s)"),
            x = grid::unit(0.1, "npc"), y = grid::unit(0.24, "npc"),
            gp = grid::gpar(fontsize = text_size))

  dev.off()                      

  if(detailed){
    # strip cosine similarity values
    inds <- grep('_c', colnames(df))
    inds_rm <- grep('_c_diff|_ml', colnames(df))
    df_cos <- df[, inds[-na.omit(match(inds_rm, inds))]]
    df_cos <- df_cos[, grep('Signature', colnames(df_cos))]

    # strip likelihood values
    inds <- grep('_ml', colnames(df))  
    df_ml <- df[, inds]
    df_ml <- df_ml[, grep('Signature', colnames(df_ml))]
  }
}