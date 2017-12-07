library(ggplot2)

# write a function that takes two input csv files
# and two signature names and calculates the sp an sn
# makes a plot

calc_sn_sp_matching <- function(df, signame1, signame2, matching){
  df[, matching] <- as.character(df[, matching])
  df$truth <- as.character(df$truth)

  total1 <- sum(df$truth == signame1)
  total2 <- sum(df$truth == signame2)
  tp1 <- sum(df$truth == signame1 &
                df[,matching] == signame1)
  tp2 <- sum(df$truth == signame2 &
                df[, matching] == signame2)
  fp1 <- sum(df$truth == signame2 &
                df[, matching] == signame1)
  fp2 <- sum(df$truth == signame1 &
                df[, matching] == signame2)

  sn1 <- tp1/total1
  sp1 <- (total2 - fp1)/total2
  sn2 <- tp2/total2
  sp2 <- (total1 - fp2)/total1
  return(list(sp1 = sp1, sn1 = sn1, sp2 = sp2, sn2 = sn2))
}

calc_sn_sp_cut <- function(vals_true, vals_false, cutoff_low){
  specificity <- rep(0, length(cutoff_low))
  sensitivity <- rep(0, length(cutoff_low))

  number_neg <- length(vals_false)
  number_pos <- length(vals_true)

  for(icut in 1:length(cutoff_low)){
    specificity[[icut]] <- sum(qual_vals_false > cutoff_low[[icut]])/number_neg
    sensitivity[[icut]] <- sum(qual_vals_true > cutoff_low[[icut]])/number_pos
  }
  return(data.frame(cutoff_low = cutoff_low,
                    sn = sensitivity,
                    sp = specificity))
}

sn_sp_vs_nsnv <- function(df, signame1, signame2, snv_ranges, matching){
  sp_vec1 <- rep(0, length(snv_ranges) - 1)
  sp_vec2 <- rep(0, length(snv_ranges) - 1)
  sn_vec1 <- rep(0, length(snv_ranges) - 1)
  sn_vec2 <- rep(0, length(snv_ranges) - 1)

  for(isnv in 1:(length(snv_ranges) - 1)){
    df_this <- df[df$total_snvs >= snv_ranges[[isnv]] &
                  df$total_snvs < snv_ranges[[isnv + 1]],]
    sn_sp <- calc_sn_sp_matching(df_this, signame1, signame2, matching)
    sp_vec1[[isnv]] <- sn_sp$sp1
    sp_vec2[[isnv]] <- sn_sp$sp2
    sn_vec1[[isnv]] <- sn_sp$sn1
    sn_vec2[[isnv]] <- sn_sp$sn2
  }

  df1 <- data.frame(sp = sp_vec1, 
                    sn = sn_vec1, 
                    nsnv_low = snv_ranges[1:(length(snv_ranges) - 1)],
                    nsnv_high = snv_ranges[2:length(snv_ranges)],
                    truth = signame1,
                    method = matching)
  df2 <- data.frame(sp = sp_vec2, 
                    sn = sn_vec2, 
                    nsnv_low = snv_ranges[1:(length(snv_ranges) - 1)],
                    nsnv_high = snv_ranges[2:length(snv_ranges)],
                    truth = signame2,
                    method = matching)
  df <- rbind(df1, df2)
  return(df)
}

plot_sn_sp <- function(file1, 
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

  input1 <- read.csv(file1)
  input2 <- read.csv(file2)
  
  input1$truth <- signame1
  input2$truth <- signame2
  if(do_cutoff){
    input1_l <- input1[input1$max_l > cutoff_l, ]
    input2_l <- input2[input2$max_l > cutoff_l, ]

    input1_c <- input1[input1$max_c > cutoff_c, ]
    input2_c <- input2[input2$max_c > cutoff_c, ]
  }else{
    input1_l <- input1
    input2_l <- input2
    input1_c <- input1
    input2_c <- input2
  }

  
  merged <- rbind(input1, input2)
  merged_l <- rbind(input1_l, input2_l)
  merged_c <- rbind(input1_c, input2_c)
 
  if(dependence == "total_snvs"){
    df_l <- sn_sp_vs_nsnv(merged_l, signame1, signame2, snv_ranges, 'sig_max_l')
    df_c <- sn_sp_vs_nsnv(merged_c, signame1, signame2, snv_ranges, 'sig_max_c')   
    df <- rbind(df_l, df_c) 
   
    plot <- ggplot(df, aes(x = 1 - sp, y = sn)) 
    plot <- plot + geom_line(aes(color = method))
    plot <- plot + facet_wrap( ~ truth)
    plot <- plot + theme_bw()
    plot <- plot + xlab('FPR') + ylab('Sensitivity')
    ggsave(plot, file = sprintf('sn_sp_%s_%s_%s.jpg', dependence, signame1, signame2))

    write.table(df, 'test.txt', row.names = F, quote = F)
  }

  if(dependence == "cut_off"){
    if(with_matching){
      scale_df_l <- sn_sp_vs_nsnv(merged, 
                                  signame1, 
                                  signame2, 
                                  c(min(merged$total_snvs, max(merged$total_snvs))), 
                                  'sig_max_l')
      scale_sp_l <- scale_df_l$sp1
      scale_sn_l <- scale_df_l$sn1
      scale_df_c <- sn_sp_vs_nsnv(merged, 
                                  signame1, 
                                  signame2, 
                                  c(min(merged$total_snvs, max(merged$total_snvs))), 
                                  'sig_max_c')
      scale_sp_c <- scale_df_c$sp1
      scale_sn_c <- scale_df_c$sn1

      df_l <- calc_sn_sp_cut(input1$max_l[input1$sig_max_l == signame1], input2$max_l[input2$sig_max_l == signame1], cutoff_low)
      df_c <- calc_sn_sp_cut(input1$max_c[input1$sig_max_c == signame1], input2$max_c[input2$sig_max_c == signame1], cutoff_low)
      df_l$sn <- df_l$sn*scale_sn_l
      df_l$sp <- df_l$sn*scale_sp_l
      df_c$sn <- df_c$sn*scale_sn_c
      df_c$sp <- df_c$sn*scale_sp_c

      df_l$method <- 'sig_max_l'
      df_c$method <- 'sig_max_c'
    }else{
      scale_sp_l <- 1
      scale_sn_l <- 1
      scale_sp_c <- 1
      scale_sn_c <- 1
      df_l <- calc_sn_sp_cut(input1[, paste0(signame1, '_l')], input2[, paste0(signame1, '_l')], cutoff_low)
      df_c <- calc_sn_sp_cut(input1[, paste0(signame1, '_c')], input2[, paste0(signame1, '_c')], cutoff_low)
    }
    df <- rbind(df_l, df_c)

    plot <- ggplot(df, aes(x = sp, y = sn, color = method)) + geom_line()
    plot <- plot + theme_bw() + xlab('relative FPR') + ylab('relative sensitivity')
    if(with_matching) plot <- plot + labs(title = sprintf('Samples matched to %s', signame1)) 
    else plot <- plot + labs('No prior matching')
    ggsave(plot, file = sprintf('sensitivity_specificity_vs_cutoff_%s_%s_with_matching%d.jpg', signame1, signame2, with_matching))
  }
}
