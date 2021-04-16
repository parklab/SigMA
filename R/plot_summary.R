#' Generates summary plot 
#'
#' @param file the csv file produced by SigMA

plot_summary <- function(file = NULL, do_mva = T, data = NULL, tumor_type = NULL){
  df <- read.csv(file)

  if(do_mva) df$pass <- df$pass_mva_strict
  else df$pass <- df$pass_ml
   

  # read SigMA settings from file name
  
  platform_names_collapsed <- lapply(platform_names, function(x){gsub(x, pattern  = " ", replace = "")})
  if(is.null(tumor_type)){
    tumor_type <- unlist(strsplit(unlist(strsplit(file, 
                                                 split = 'tumortype_'))[[2]],
                          split = '_platform_'))[[1]]
  }

  if(is.null(data)){
 
    platform <- unlist(strsplit(unlist(strsplit(file, 
                                                split = 'platform_'))[[2]],
                       split = '_cf'))[[1]]


    if(platform %in% platform_names_collapsed){
      data <- names(platform_names)[which(platform_names_collapsed == platform)]
    }
    else if(grepl('Custom', platform)){
      data <- gsub('Custom', platform)
      custom <- T
    }
    else{
      stop('data setting cannot be extracted from file provide data and tumor_type arguments to plot_summary function')
    }
  }

  if(platform %in% platform_names_collapsed){
    cutoff_strict <- cutoffs_strict[[data]][[tumor_type]]
    cutoff <- cutoffs[[data]][[tumor_type]]
  }
  else if(custom){
    file_path <- system.file(paste0("extdata/gbm_models/", data, ".rda"),
                             package="SigMA")
    if(!file.exists(file_path)) stop('the provided custom data setting does not exist')
    if(!(tumor_type %in% names(cutoffs_custom))) stop('tumor_type not available for the data setting')
    cutoffs_this <- cutoffs_custom[[tumor_type]]
    cutoffs_strict_this <- cutoffs_strict_custom[[tumor_type]]
  }else{
    stop('the data parameter provided does not exist')
  }

   

  text_size = 10

  # 3-base distributions
  inds_pos <- which(df$pass)
  pos_tribase <- plot_tribase_dist(as.data.frame(colSums(df[inds_pos, 1:96])),
                                     signame = paste0('  Aggregated mutations from ', 
                                                      sum(df$pass), ' Sig3+ sample(s)'))
  inds_neg <- which(!df$pass)
  neg_tribase <- plot_tribase_dist(as.data.frame(colSums(df[inds_neg, 1:96])),
                                   signame = paste0('  Aggregated mutations from ', 
                                                    sum(df$pass == F), ' Sig3- sample(s)'))
  neg_tribase <- neg_tribase + ggplot2::theme(legend.position = "none")

  # cosine density distribution
  colors_cos <- c(col_pos_neg[[2]], col_pos_neg[[1]])

  if(length(unique(df$pass)) == 1){
    if(unique(df$pass))
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
    cos_vals_pos <- df$Signature_3_c[df$pass]
    cos_vals_neg <- df$Signature_3_c[!df$pass]
    counts_cos_pos[[ibin]] <- sum(cos_vals_pos >= cos_bins[[ibin]] & cos_vals_pos < cos_bins[[ibin + 1]]) 
    counts_cos_neg[[ibin]] <- sum(cos_vals_neg >= cos_bins[[ibin]] & cos_vals_neg < cos_bins[[ibin + 1]]) 

    bin_center[[ibin]] <- (cos_bins[[ibin]] + cos_bins[[ibin + 1]])/2
  }

  df_cos <- data.frame(cos = numeric(), count = integer(), group = character())
  if(sum(counts_cos_pos) > 0){  
    df_cos <- rbind(df_cos, data.frame(cos = bin_center, count = counts_cos_pos, group = 'Sig3+'))
  }
  if(sum(counts_cos_neg) > 0){
    df_cos <- rbind(df_cos, data.frame( cos = bin_center, count = counts_cos_neg, group = 'Sig3-' ))
  }

  plot_cos <- ggplot2::ggplot(df_cos, ggplot2::aes(x = cos, y = count))
  plot_cos <- plot_cos + ggplot2::geom_bar(stat = 'identity', position = 'dodge', 
                                           ggplot2::aes(fill = group))
  plot_cos <- plot_cos + ggplot2::scale_fill_manual(values = colors_cos)
  plot_cos <- plot_cos + ggplot2::theme_bw()
  plot_cos <- plot_cos + ggplot2::theme(
                                legend.title = ggplot2::element_blank(),
                                legend.justification = c(1,1),
                                legend.text = ggplot2::element_text(size = text_size, face = "bold"),
                                legend.position="top",
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank())

  plot_cos <- plot_cos + ggplot2::theme(axis.title.x = ggplot2::element_text(size = text_size),
                                                    axis.title.y = ggplot2::element_text(size = text_size),
                                                    axis.text.x = ggplot2::element_text(size = text_size),
                                                    axis.text.y = ggplot2::element_text(size = text_size))

  plot_cos <- plot_cos + ggplot2::xlab('Cos simil of Sig3')
#  plot_cos <- plot_cos + ggplot2::xlim(0, 1)
 
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

  if(!('Signature_3_ml' %in% colnames(df_ml))){
    inds_3ml <- grep('Signature_3', colnames(df_ml))
    df_ml <- df_ml[, c(colnames(df_ml), 'tumor')]

    if(length(inds_3ml) > 1) df_ml$Signature_3_ml <- rowSums(df_ml[,inds_3ml])
    else if(length(inds_3ml) == 1) colnames(df_ml)[[inds_3ml]] <- 'Signature_3_ml'
  }

  df_ml$tumor <- as.character(df$tumor)
  df_ml$pass <- df$pass
  if(do_mva){
    df_ml$Signature_3_mva <- df$Signature_3_mva 
  }
  ml_bin_width = 0.1
  ml_bins <- seq(0, 1.001, by = ml_bin_width)
  counts_ml_pos <- rep(0, length(ml_bins) - 1)
  counts_ml_neg <- rep(0, length(ml_bins) - 1)
  bin_center <- rep(0, length(ml_bins) - 1)
  
  ml_vals_pos <- df_ml$Signature_3_ml[df_ml$pass]
  ml_vals_neg <- df_ml$Signature_3_ml[!df_ml$pass]

  for(ibin in seq_len(length(ml_bins) - 1)){
    ml_bins[[ibin]]
    counts_ml_pos[[ibin]] <- sum(ml_vals_pos >= ml_bins[[ibin]] & ml_vals_pos < ml_bins[[ibin + 1]]) 
    counts_ml_neg[[ibin]] <- sum(ml_vals_neg >= ml_bins[[ibin]] & ml_vals_neg < ml_bins[[ibin + 1]]) 
    bin_center[[ibin]] <- (ml_bins[[ibin]] + ml_bins[[ibin + 1]])/2
  }

  df_ml_binned <- data.frame(ml = numeric(), count = integer(), group = character())
  if(sum(counts_ml_pos) > 0){
    df_ml_binned <- rbind(df_ml_binned, data.frame(ml = bin_center, count = counts_ml_pos, group = 'Sig3+'))
  }
  if(sum(counts_ml_neg) > 0){
    df_ml_binned <- rbind(df_ml_binned, data.frame( ml = bin_center, count = counts_ml_neg, group = 'Sig3-' ))
  }


  plot_ml <- ggplot2::ggplot(df_ml_binned, ggplot2::aes(x = ml, y = count))
  plot_ml <- plot_ml + ggplot2::geom_bar(stat = 'identity', position = 'dodge', 
                                           ggplot2::aes(fill = group))
  plot_ml <- plot_ml + ggplot2::scale_fill_manual(values = colors_cos) 
  # here add a line to strip away the legend
  plot_ml <- plot_ml + ggplot2::theme_bw()
  plot_ml <- plot_ml + ggplot2::theme(
                                legend.title = ggplot2::element_blank(),
                                legend.justification = c(1,1),
                                legend.text = ggplot2::element_text(size = text_size, face = "bold"),
                                legend.position="top",
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank())
  plot_ml <- plot_ml + ggplot2::theme(axis.title.x = ggplot2::element_text(size = text_size),
                                                  axis.title.y = ggplot2::element_text(size = text_size),
                                                  axis.text.x = ggplot2::element_blank(),
                                                  axis.text.y = ggplot2::element_text(size = text_size))


  plot_ml <- plot_ml + ggplot2::xlab('Likelihood of Sig3')
  plot_ml <- plot_ml + ggplot2::xlim(0, 1)  


  if(do_mva){
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


    if(do_mva) df_ml <- transform(df_ml, tumor_label = factor(tumor_label, levels = df_ml$tumor_label[order(df_ml$Signature_3_mva)]))
    else  df_ml <- transform(df_ml, tumor_label = factor(tumor_label, levels = df_ml$tumor_label[order(df_ml$Signature_3_ml)]))

    plot_score <- ggplot2::ggplot(df_ml, 
                                  ggplot2::aes(x = tumor_label, 
                                               y = Signature_3_mva))
    plot_score <- plot_score + ggplot2::geom_bar(ggplot2::aes(fill = group_fac),
                                                 stat = 'identity')
    plot_score <- plot_score + ggplot2::scale_fill_manual(values = as.character(df_color$color[match(unique(df_ml$group), 
                                                                                               df_color$group)]))
    plot_score <- plot_score + ggplot2::theme_bw()
    plot_score <- plot_score + ggplot2::theme(
                                legend.title = ggplot2::element_blank(),
                                legend.justification = c(1,1),
                                legend.text = ggplot2::element_text(size = text_size, face = "bold"),
                                legend.position="top",
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank())
 
    plot_score <- plot_score + ggplot2::theme(panel.border = ggplot2::element_blank(),
                                              axis.line.y = ggplot2::element_line(size = 0.3, 
                                                                                  linetype = "solid", 
                                                                                  colour = "black"),
                                              axis.line.x = ggplot2::element_line(size = 0.3, 
                                                                                  linetype = "solid", 
                                                                                  colour = "black"),
                                              axis.text.x = ggplot2::element_text(size = text_size), 
                                              axis.title.y = ggplot2::element_text(size = text_size),
                                              axis.text.y = ggplot2::element_text(size = text_size),
                                              legend.position = 'top')
    plot_score <- plot_score + ggplot2::ylab('MVA score') + ggplot2::xlab('')
    plot_score <- plot_score + ggplot2::geom_hline(ggplot2::aes(yintercept = cutoff))
    plot_score <- plot_score + ggplot2::geom_hline(ggplot2::aes(yintercept = cutoff_strict), 
                                                   linetype = 'dashed')
  }
  # Signature 3 exposures 
  colors_exp <- colors_cos
  if(length(unique(df$pass[df$exp_sig3 > 0])) == 1){
    if(unique(df$pass[df$exp_sig3 > 0]))
      colors_exp <- colors_cos[seq(from = length(colors_cos), to = 1, by = -1)]
  }

  exp_bin_width <- (max(df$exp_sig3) + 1)/10
  exp_bins <- seq(0.001, max(df$exp_sig3) + 1, by = exp_bin_width)
  counts_exp_pos <- rep(0, length(exp_bins) - 1)
  counts_exp_neg <- rep(0, length(exp_bins) - 1)
  bin_center <- rep(0, length(exp_bins) - 1)
  
  for(ibin in seq_len(length(exp_bins) - 1)){
    exp_vals_pos <- df$exp_sig3[df$pass]
    exp_vals_neg <- df$exp_sig3[!df$pass]
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
  return(list(trinucleotide_pos = pos_tribase, trinucleotide_neg = neg_tribase, plot_cos = plot_cos, plot_ml = plot_ml, plot_score = plot_score))
}