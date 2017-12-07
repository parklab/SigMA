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
#' @param cutoff_l if do_cutoff is set to TRUE this parameter is 
#' used for likelihood method
#' @param cutoff_c if do_cutoff is set to TRUE this parameter is 
#' used for cosine method
#' @param snv_ranges ranges used for plotting  
#'
#' @examples
#' plot_sn_fp('genomes_sig3_mc_output.csv', 
#'           'genomes_sig5_mc_output.csv', 
#'           'Signature_3', 
#'           'Signature_5')

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

calc_sn_fp_cut <- function(vals_true, vals_false, cutoff_low){
  falsepos <- rep(0, length(cutoff_low))
  sensitivity <- rep(0, length(cutoff_low))
  
  number_neg <- length(vals_false)
  number_pos <- length(vals_true)
  
  for(icut in 1:length(cutoff_low)){
    falsepos[[icut]] <- sum(vals_false > cutoff_low[[icut]])/number_neg
    sensitivity[[icut]] <- sum(vals_true > cutoff_low[[icut]])/number_pos
  }
  
  return(data.frame(cutoff_low = cutoff_low,
                    sn = sensitivity,
                    fp = falsepos))
}

sn_fp_vs_nsnv <- function(df, signame1, signame2, snv_ranges, matching, cutoff){

  fp_vec1 <- rep(0, length(snv_ranges) - 1)
  fp_vec2 <- rep(0, length(snv_ranges) - 1)
  sn_vec1 <- rep(0, length(snv_ranges) - 1)
  sn_vec2 <- rep(0, length(snv_ranges) - 1)
  
  for(isnv in 1:(length(snv_ranges) - 1)){
    df_this <- df[df$total_snvs >= snv_ranges[[isnv]] &
                  df$total_snvs < snv_ranges[[isnv + 1]],]
    sn_fp <- calc_sn_fp_matching(df_this, signame1, signame2, matching, cutoff)
    fp_vec1[[isnv]] <- sn_fp$fp1
    fp_vec2[[isnv]] <- sn_fp$fp2
    sn_vec1[[isnv]] <- sn_fp$sn1
    sn_vec2[[isnv]] <- sn_fp$sn2
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

plot_sn_fp <- function(file1, 
                       file2, 
                       signame1, 
                       signame2, 
                       dependence = "total_snvs", 
                       snv_ranges = c(5, 15, 25, 35, 45, 55, 
                                      65, 75, 85, 100),
                       cutoff_low = c(0, 0.002, 0.005, 0.007, 
                                      0.01, 0.02, 0.05, 0.1, 
                                      0.15, 0.2, 0.25, 0.3, 
                                      0.35, 0.4, 0.45, 0.5, 
                                      0.55, 0.6, 0.65, 0.7, 
                                      0.75, 0.8, 0.85, 0.875, 
                                      0.9, 0.92, 0.94, 0.96, 
                                      0.98, 0.99, 0.995, 0.997),
                       with_matching = FALSE,
                       do_cutoff = FALSE,
                       cutoff_c = 0.15, 
                       cutoff_l = 0.25)
{
  library(ggplot2)
  color_l_c <- c('#76ACF1', '#0B148B')

  input1 <- read.csv(file1)
  input2 <- read.csv(file2)
  
  input1$truth <- signame1
  input2$truth <- signame2

  merged <- rbind(input1, input2)
 
  if(!do_cutoff){
    cutoff_c <- 0
    cutoff_l <- 0
  }
  if(dependence == "total_snvs"){
    
    df_l <- sn_fp_vs_nsnv(merged, signame1, signame2, snv_ranges, 'sig_max_l', cutoff_l)
    df_c <- sn_fp_vs_nsnv(merged, signame1, signame2, snv_ranges, 'sig_max_c', cutoff_c)   
    df <- rbind(df_l, df_c) 
    
    df$method <- as.character(df$method)

    plot <- ggplot(df, aes(x = fp, y = sn)) 
    plot <- plot + geom_line(aes(color = method))
    plot <- plot + scale_color_manual(values = color_l_c)
    plot <- plot + facet_wrap( ~ truth)
    plot <- plot + theme_bw()
    plot <- plot + xlab('FPR') + ylab('Sensitivity')
    plot <- plot + xlim(0, 0.15)
    plot <- plot + ylim(0, 1)
    ggsave(plot, 
           file = sprintf('sn_fp_%s_%s_%s_cutoff%d.jpg', 
                           dependence, 
                           signame1, 
                           signame2, 
                           do_cutoff),
           width = 8,
           height = 4)
  }

  # no selection is applied on number of SNVs cut off is varied
  if(dependence == "cutoff"){
    if(with_matching){  # first matching is done then cutoff is varied
      scale_df_l <- sn_fp_vs_nsnv(merged, 
                                  signame1, 
                                  signame2, 
                                  c(min(merged$total_snvs), max(merged$total_snvs)), 
                                  'sig_max_l',
                                  cutoff = 0)

      scale_df_l$truth <- as.character(scale_df_l$truth)

      scale_fp_l <- scale_df_l$fp[scale_df_l$truth == signame1]
      scale_sn_l <- scale_df_l$sn[scale_df_l$truth == signame1]

      scale_df_c <- sn_fp_vs_nsnv(merged, 
                                  signame1, 
                                  signame2, 
                                  c(min(merged$total_snvs), max(merged$total_snvs)), 
                                  'sig_max_c',
                                  cutoff = 0)

      scale_df_c$truth <- as.character(scale_df_c$truth)
      scale_fp_c <- scale_df_c$fp[scale_df_c$truth == signame1]
      scale_sn_c <- scale_df_c$sn[scale_df_c$truth == signame1]

      df_l <- calc_sn_fp_cut(input1$max_l[input1$sig_max_l == signame1], 
                             input2$max_l[input2$sig_max_l == signame1], 
                             cutoff_low)
      df_c <- calc_sn_fp_cut(input1$max_c[input1$sig_max_c == signame1], 
                             input2$max_c[input2$sig_max_c == signame1], 
                             cutoff_low)

      df_l$sn <- df_l$sn*scale_sn_l
      df_l$fp <- df_l$fp*scale_fp_l

      df_c$sn <- df_c$sn*scale_sn_c
      df_c$fp <- df_c$fp*scale_fp_c

      df_l$method <- 'sig_max_l'
      df_c$method <- 'sig_max_c'
    }else{  # no prior matching is applied just cutoff is varied
      scale_fp_l <- 1
      scale_sn_l <- 1
      scale_fp_c <- 1
      scale_sn_c <- 1
      df_l <- calc_sn_fp_cut(input1[, paste0(signame1, '_l')], 
                             input2[, paste0(signame1, '_l')], 
                             cutoff_low)
      df_c <- calc_sn_fp_cut(input1[, paste0(signame1, '_c')], 
                             input2[, paste0(signame1, '_c')], 
                             cutoff_low)

      df_l$method <- 'sig_max_l'
      df_c$method <- 'sig_max_c'
    }
    df <- rbind(df_l, df_c)
    
    df$method <- as.character(df$method)

    plot <- ggplot(df, aes(x = fp, y = sn, color = method)) + geom_line()
    plot <- plot + theme_bw() + xlab('FPR') + ylab('Sensitivity')
    plot <- plot + scale_color_manual(values = color_l_c)

    if(with_matching) plot <- plot + labs(title = sprintf('Samples matched to %s', signame1)) 
    else plot <- plot + labs(title = 'No prior matching')

    ggsave(plot, 
           file = sprintf('sensitivity_falsepos_vs_cutoff_%s_%s_with_matching%d.jpg', 
                          signame1, 
                          signame2, 
                          with_matching), 
           height = 4, 
           width = 5)
  }
}
