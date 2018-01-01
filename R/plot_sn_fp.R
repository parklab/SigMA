#' plot sensitivity vs false positive rate for different settings:
#' * by changing number of total SNVs with two options (1) with
#' a cut off on the likelihood/cosine similarity (2) without
#' * by changing the cutoff for two options (1) by requiring 
#' matching to the signature of interest first (2) only using 
#' the likelihood/cosine similarity values
#'
#' @param file1 input file with a single signature whose name
#' is given by signame1
#' @param file2 input file with a single signature whose name 
#' is given by signame2
#' @param signame1 main signature of interest
#' @param signame2 signature to test signame1 against
#' @param dependence can be "total_snvs" or "cut_off", it is 
#' the parameter to vary to span the false positive rate and 
#' sensitivity
#' @param snv_ranges bins on total SNVs used if dependence 
#' is set "to total_snvs"
#' @param cutoff_low the values lower threshold on cutoff, 
#' i.e. likelihood or cosine similarity, used if dependence
#' is set to "cutoff"
#' @param with_matching determines whether cutoff dependent plot 
#' will have prior matching
#' @param do_cutoff for the "total_snv" plot TRUE adds additional
#' cutoff on likelihood or cosine similarity on top of matching
#'
#' @examples
#' plot_sn_fp('genomes_sig3_mc_output.csv 
#'           'genomes_sig5_mc_output.csv', 
#'           'Signature_3', 
#'           'Signature_5')


plot_sn_fp <- function(file1, 
                       file2, 
                       signame1, 
                       signame2, 
                       dependence = "total_snvs", 
                       snv_ranges = c(5, 25, 45, 65, 85, 100),
                       cutoff_low = c(0, 0.01, 0.02, 0.05, 0.1, 
                                      0.15, 0.2, 0.25, 0.3, 
                                      0.35, 0.4, 0.45, 0.5, 
                                      0.55, 0.6, 0.65, 0.7, 
                                      0.75, 0.8, 0.85, 0.875, 
                                      0.9, 0.92, 0.94, 0.96, 
                                      0.98, 0.99, 0.995, 0.997),
                       with_matching = FALSE,
                       do_cutoff = FALSE,
                       max_allowed_fp = 0.5, 
                       plot_dir = '.',
                       write_output = F,
                       output_file = 'test_sn_fp_output.csv',
                       output_dir = '.')
{
  print('snv ranges')
  print(snv_ranges)

  library(ggplot2)
  color_l_c <- c('#76ACF1', '#0B148B')

  input1 <- read.csv(file1)
  input2 <- read.csv(file2)
  
  if(sum(colnames(input1) == 'total_snvs') == 0) input1$total_snvs <- rowSums(input1[, 1:96])
  if(sum(colnames(input2) == 'total_snvs') == 0) input2$total_snvs <- rowSums(input2[, 1:96])
 
  input1$truth <- signame1
  input2$truth <- signame2

  merged <- rbind(input1, input2)
 
  if(!do_cutoff){ # if no cut off should be applied set it to 0
    cutoff_c <- rep(0, length(snv_ranges) - 1)
    cutoff_l <- rep(0, length(snv_ranges) - 1)
  }else{   
    list_cutoff <- tune_cutoff_vs_nsnv(input1 = input1, 
                                       input2 = input2, 
                                       signame1 = signame1, 
                                       signame2 = signame2, 
                                       snv_ranges = snv_ranges, 
                                       cutoff_low = cutoff_low, 
                                       with_matching = with_matching,  
                                       max_allowed_fp = max_allowed_fp, 
                                       output_dir = output_dir, 
                                       output_file = output_file,
                                       plot_dir = plot_dir)
    cutoff_c <- list_cutoff$cutoff_c
    cutoff_l <- list_cutoff$cutoff_l
  }


  if(dependence == "total_snvs"){
    # calculate sensitivity and false positive rate as a function of NSNV
    
    # for the likelihood method    
    df_l <- sn_fp_vs_nsnv(merged, signame1, signame2, snv_ranges, 'sig_max_l', cutoff_l, with_matching)
    # for the cosine similarity method
    df_c <- sn_fp_vs_nsnv(merged, signame1, signame2, snv_ranges, 'sig_max_c', cutoff_c, with_matching)   
    df <- rbind(df_l, df_c) 
    
    # set method to char for plotting
    df$method <- as.character(df$method)
    
    if(write_output) write.table(df, sprintf('%s/%s_%s', output_dir, dependence, output_file), row.names = F, quote = F, sep = ',')

    # plot sensitivity vs FPR 
    plot <- ggplot2::ggplot(df[df$truth == signame1,], aes(x = fp, y = sn)) 
    plot <- plot + ggplot2::geom_line(aes(color = method))
    plot <- plot + ggplot2::scale_color_manual(values = color_l_c)
    plot <- plot + ggplot2::theme_bw()
    plot <- plot + ggplot2::xlab('FPR') + ggplot2::ylab('Sensitivity')
    plot <- plot + ggplot2::ylim(0, 1)

    ggplot2::ggsave(plot, 
                    file = sprintf('%s/sn_fp_%s_%s_%s_%s_cutoff%d.jpg', 
                                   plot_dir, 
                                   dependence,
                                   output_file,  
                                   signame1, 
                                   signame2, 
                                   do_cutoff),
                    width = 5,
                    height = 4)
    
    # plot sensitivity as a function of NSNV
    plot <- ggplot2::ggplot(df[df$truth == signame1,], aes(x = (nsnv_low + nsnv_high)/2., y = sn))
    plot <- plot + ggplot2::geom_line(aes(color = method))
    plot <- plot + ggplot2::scale_color_manual(values = color_l_c)
    plot <- plot + ggplot2::theme_bw()
    plot <- plot + ggplot2::xlab('# SNV') + ggplot2::ylab('Sensitivity')
    ggplot2::ggsave(plot, 
                   file = sprintf('%s/sn_nsnv_%s_%s_%s_%s_cutoff%d.jpg',
                                 plot_dir, 
                                 dependence,
                                 output_file, 
                                 signame1, 
                                 signame2, 
                                 do_cutoff),
                   width = 5,
                   height = 4)


    # plot specificity as a function of NSNV
    plot <- ggplot2::ggplot(df[df$truth == signame1,], aes(x = (nsnv_low + nsnv_high)/2., y = fp)) 
    plot <- plot + ggplot2::geom_line(aes(color = method))
    plot <- plot + ggplot2::scale_color_manual(values = color_l_c)
    plot <- plot + ggplot2::theme_bw()
    plot <- plot + ggplot2::xlab('# SNV') + ggplot2::ylab('FPR')
    ggplot2::ggsave(plot, 
                    file = sprintf('%s/fp_nsnv_%s_%s_%s_%s_cutoff%d.jpg', 
                                   plot_dir,
                                   dependence, 
                                   output_file,
                                   signame1, 
                                   signame2, 
                                   do_cutoff),
                    width = 5,
                    height = 4)

  }

  # no selection is applied on number of SNVs cut off is varied
  if(dependence == "cutoff"){
    if(with_matching){  # first matching is done then cutoff is varied
      # first calculate the effect of matching on sensitivity and false 
      # positive rate, these values are applied as overall scales later
      # when the scan over cutoff values is made

      # for the likelihood method
      scale_df_l <- sn_fp_vs_nsnv(merged, 
                                  signame1, 
                                  signame2, 
                                  c(min(merged$total_snvs), max(merged$total_snvs)), 
                                  'sig_max_l',
                                  cutoff = 0,
                                  with_matching)

      scale_df_l$truth <- as.character(scale_df_l$truth)

      scale_fp_l <- scale_df_l$fp[scale_df_l$truth == signame1]
      scale_sn_l <- scale_df_l$sn[scale_df_l$truth == signame1]

      #for the cosine similarity method
      scale_df_c <- sn_fp_vs_nsnv(merged, 
                                  signame1, 
                                  signame2, 
                                  c(min(merged$total_snvs), max(merged$total_snvs)), 
                                  'sig_max_c',
                                  cutoff = 0, 
                                  with_matching)

      scale_df_c$truth <- as.character(scale_df_c$truth)
      scale_fp_c <- scale_df_c$fp[scale_df_c$truth == signame1]
      scale_sn_c <- scale_df_c$sn[scale_df_c$truth == signame1]

      # calculate relative sensitivity and false positive rate by
      # changing the cut off 
   
      # for the likelihood method
      df_l <- calc_sn_fp_cut(input1$max_l[input1$sig_max_l == signame1], 
                             input2$max_l[input2$sig_max_l == signame1], 
                             cutoff_low)
      # for the cosine similarity method
      df_c <- calc_sn_fp_cut(input1$max_c[input1$sig_max_c == signame1], 
                             input2$max_c[input2$sig_max_c == signame1], 
                             cutoff_low)

      #apply the scales obtained from just matching
      df_l$sn <- df_l$sn*scale_sn_l
      df_l$fp <- df_l$fp*scale_fp_l

      df_c$sn <- df_c$sn*scale_sn_c
      df_c$fp <- df_c$fp*scale_fp_c

      df_l$method <- 'sig_max_l'
      df_c$method <- 'sig_max_c'

    }else{  # no prior matching is applied just cutoff is varied
      # calculate sensitivity and false positive rate by varying 
      # the cutoff values

      # for likelihood method
      df_l <- calc_sn_fp_cut(input1[, paste0(signame1, '_l')], 
                             input2[, paste0(signame1, '_l')], 
                             cutoff_low)
      # for cosine similarity method
      df_c <- calc_sn_fp_cut(input1[, paste0(signame1, '_c')], 
                             input2[, paste0(signame1, '_c')], 
                             cutoff_low)

      df_l$method <- 'sig_max_l'
      df_c$method <- 'sig_max_c'
    }

    df <- rbind(df_l, df_c)

    # set method to be char instead of factor for consistency of colors
    df$method <- as.character(df$method)

    if(write_output) write.table(df, sprintf('%s/%s_%s', output_dir, dependence, output_file), row.names = F, quote = F, sep = ',')

    #plot the roc curve
    plot <- ggplot2::ggplot(df, aes(x = fp, y = sn, color = method)) + ggplot2::geom_line()
    plot <- plot + ggplot2::theme_bw() + ggplot2::xlab('FPR') + ggplot2::ylab('Sensitivity')
    plot <- plot + ggplot2::scale_color_manual(values = color_l_c)
    if(with_matching) plot <- plot + ggplot2::labs(title = sprintf('Samples matched to %s', signame1)) 
    else plot <- plot + ggplot2::labs(title = 'No prior matching')

    ggplot2::ggsave(plot, 
                    file = sprintf('%s/sensitivity_falsepos_vs_cutoff_%s_%s_%s_with_matching%d.jpg', 
                                   plot_dir,
                                   output_file, 
                                   signame1, 
                                   signame2, 
                                   with_matching), 
                    height = 4, 
                    width = 5)
  }
}

calc_sn_fp_matching <- function(df, signame1, signame2, matching, cutoff){
  df[, matching] <- as.character(df[, matching])
  df$truth <- as.character(df$truth)

  total1 <- sum(df$truth == signame1)
  total2 <- sum(df$truth == signame2)

  tp1 <- sum(df$truth == signame1 &
             df[, matching] == signame1 &
             df[, paste0('max', strsplit(matching, 'sig_max')[[1]][[2]])] >= cutoff)
  tp2 <- sum(df$truth == signame2 &
             df[, matching] == signame2 &
             df[, paste0('max', strsplit(matching, 'sig_max')[[1]][[2]])] >= cutoff)
  fp1 <- sum(df$truth == signame2 &
             df[, matching] == signame1 &
             df[, paste0('max', strsplit(matching, 'sig_max')[[1]][[2]])] >= cutoff)
  fp2 <- sum(df$truth == signame1 &
             df[, matching] == signame2 &
             df[, paste0('max', strsplit(matching, 'sig_max')[[1]][[2]])] >= cutoff)
   

  sn1 <- tp1/total1
  fp1 <- fp1/total2
  sn2 <- tp2/total2
  fp2 <- fp2/total1
  return(list(fp1 = fp1, sn1 = sn1, fp2 = fp2, sn2 = sn2))
}

calc_sn_fp_cut <- function(vals_true, 
                           vals_false, 
                           cutoff_low){
  falsepos <- rep(0, length(cutoff_low))
  sensitivity <- rep(0, length(cutoff_low))
  
  number_neg <- length(vals_false)
  number_pos <- length(vals_true)
  
  
  for(icut in 1:length(cutoff_low)){
    if(number_neg != 0) falsepos[[icut]] <- sum(vals_false > cutoff_low[[icut]])/number_neg
    if(number_pos != 0)sensitivity[[icut]] <- sum(vals_true > cutoff_low[[icut]])/number_pos
  }
  
  return(data.frame(cutoff_low = cutoff_low,
                    sn = sensitivity,
                    fp = falsepos))
}



sn_fp_vs_nsnv <- function(df, 
                          signame1, 
                          signame2, 
                          snv_ranges, 
                          matching, 
                          cutoff, 
                          with_matching){

  fp_vec1 <- rep(0, length(snv_ranges) - 1)
  fp_vec2 <- rep(0, length(snv_ranges) - 1)
  sn_vec1 <- rep(0, length(snv_ranges) - 1)
  sn_vec2 <- rep(0, length(snv_ranges) - 1)
  
  for(isnv in 1:(length(snv_ranges) - 1)){
    df_this <- df[df$total_snvs >= snv_ranges[[isnv]] &
                  df$total_snvs < snv_ranges[[isnv + 1]],]
    if(with_matching){
      sn_fp <- calc_sn_fp_matching(df_this, 
                                   signame1, 
                                   signame2, 
                                   matching, 
                                   cutoff[[isnv]])
      fp_vec1[[isnv]] <- sn_fp$fp1
      fp_vec2[[isnv]] <- sn_fp$fp2
      sn_vec1[[isnv]] <- sn_fp$sn1
      sn_vec2[[isnv]] <- sn_fp$sn2
    }
    else{
      df_sig1 <- calc_sn_fp_cut(df[df$truth == signame1, paste0(signame1, '_l')],
                                df[df$truth == signame2, paste0(signame1, '_l')],
                                cutoff[[isnv]])

      df_sig2 <- calc_sn_fp_cut(df[df$truth == signame2, paste0(signame2, '_l')],
                                df[df$truth == signame1, paste0(signame2, '_l')],
                                cutoff[[isnv]])

      fp_vec1[[isnv]] <- df_sig1$fp
      fp_vec2[[isnv]] <- df_sig2$fp
      sn_vec1[[isnv]] <- df_sig1$sn
      sn_vec2[[isnv]] <- df_sig2$sn
    }
     
  }

  df1 <- data.frame(fp = fp_vec1, 
                    sn = sn_vec1, 
                    nsnv_low = snv_ranges[1:(length(snv_ranges) - 1)],
                    nsnv_high = snv_ranges[2:length(snv_ranges)],
                    truth = signame1,
                    method = matching)
  df2 <- data.frame(fp = fp_vec2, 
                    sn = sn_vec2, 
                    nsnv_low = snv_ranges[1:(length(snv_ranges) - 1)],
                    nsnv_high = snv_ranges[2:length(snv_ranges)],
                    truth = signame2,
                    method = matching)
  df <- rbind(df1, df2)
  return(df)
}

tune_cutoff <- function(sn_vec, fp_vec, cutoff_low, max_allowed_fp){
  indices_keep <- which(fp_vec <= max_allowed_fp)

  sn_vec <- sn_vec[indices_keep]
  fp_vec <- fp_vec[indices_keep]

  cutoff_low <- cutoff_low[indices_keep]    
  ind_max <- which(max(sn_vec - 2 * fp_vec) == (sn_vec - 2 * fp_vec))[[1]]

  return(cutoff_low[[ind_max]])
}

tune_cutoff_vs_nsnv <- function(input1, 
                                input2, 
                                signame1, 
                                signame2, 
                                snv_ranges,  
                                cutoff_low, 
                                with_matching, 
                                max_allowed_fp, 
                                output_dir, 
                                output_file, 
                                plot_dir){

  cutoff_l <- rep(0, length(snv_ranges) - 1)
  cutoff_c <- rep(0, length(snv_ranges) - 1)
  
  for(isnv in 1:(length(snv_ranges) - 1)){
    input1_this <- input1[input1$total_snvs >= snv_ranges[[isnv]] & input1$total_snvs < snv_ranges[[isnv + 1]], ]
    input2_this <- input2[input2$total_snvs >= snv_ranges[[isnv]] & input2$total_snvs < snv_ranges[[isnv + 1]], ]

   
    if(with_matching){
      vals_pos_l <- input1_this$max_l[input1_this$sig_max_l == signame1]
      vals_neg_l <- input2_this$max_l[input2_this$sig_max_l == signame1]

      vals_pos_c <- input1_this$max_c[input1_this$sig_max_c == signame1]
      vals_neg_c <- input2_this$max_c[input2_this$sig_max_c == signame1]

      scale_sn_l <- (sum(input1_this$sig_max_l == signame1)/dim(input1_this)[[1]])
      scale_fp_l <- (sum(input2_this$sig_max_l == signame1)/dim(input2_this)[[1]])

      scale_sn_c <- (sum(input1_this$sig_max_c == signame1)/dim(input1_this)[[1]])
      scale_fp_c <- (sum(input2_this$sig_max_c == signame1)/dim(input2_this)[[1]])

    }else{

      scale_sn_l <- 1
      scale_fp_l <- 1
      scale_sn_c <- 1
      scale_fp_c <- 1
      
      vals_pos_l <- input1_this[, paste0(signame1, '_l')]
      vals_neg_l <- input2_this[, paste0(signame1, '_l')]
          
      vals_pos_c <- input1_this[, paste0(signame1, '_c')]
      vals_neg_c <- input2_this[, paste0(signame1, '_c')]  
      
    }

    df_l_1 <- calc_sn_fp_cut(vals_pos_l, 
                             vals_neg_l, 
                             cutoff_low)
    
    cutoff_l_1 <- tune_cutoff(scale_sn_l*df_l_1$sn, 
                              scale_fp_l*df_l_1$fp, 
                              cutoff_low,
                              max_allowed_fp = max_allowed_fp)
    df_c_1 <- calc_sn_fp_cut(vals_pos_c, 
                             vals_neg_c,
                             cutoff_low)

    df_l_1$method <- 'l'
    df_l_1$sn <- df_l_1$sn * scale_sn_l
    df_l_1$fp <- df_l_1$fp * scale_fp_l

    df_c_1$method <- 'c'
    df_c_1$sn <- df_c_1$sn * scale_sn_c
    df_c_1$fp <- df_c_1$fp * scale_fp_c
   
     
    df_1 <- rbind(df_l_1, df_c_1)    
    plot <- ggplot2::ggplot(df_1, aes(x = fp, y = sn)) 
    plot <- plot + ggplot2::geom_line(aes(color = method))
    plot <- plot + ggplot2::theme_bw()
    ggsave(plot, file = sprintf('%s/tune_%d_%d.pdf', 
                                 plot_dir, 
                                 min(snv_ranges), 
                                 max(snv_ranges)))

    cutoff_c_1 <- ggplot2::tune_cutoff(scale_sn_c*df_c_1$sn, 
                                       scale_fp_c*df_c_1$fp, 
                                       cutoff_low, 
                                       max_allowed_fp = max_allowed_fp)

    cutoff_l[[isnv]] <- cutoff_l_1
    cutoff_c[[isnv]] <- cutoff_c_1
  }

  write.table(data.frame(cutoff_l = cutoff_l, 
                         cutoff_c = cutoff_c, 
                         snv_low = snv_ranges[1:(length(snv_ranges) - 1)]),
              sprintf('%s/tuned_cuts_vs_nsnvs_%s', output_dir, output_file),
              row.names = F,
              col.names = T, 
              quote = F, 
              sep = ',')
  return(list(cutoff_l = cutoff_l, cutoff_c = cutoff_c))
  
}