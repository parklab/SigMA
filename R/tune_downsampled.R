#' After running downsample_random(), run this to 
#' calculate the likelihoods, sensitivity and
#' false positive rate as a function of total number
#' of snvs. Tuning as a function of number of mutations
#' is done and a file with tune cutoffs and another
#' file with sensitivity and false positive rate for
#' this cutoff is produced
#' 
#' @param directory is the directory where the outputs of
#' downsample_random() are
#' @param info_file is the file with the tumor column
#' that matches the file names in the 'directory', this 
#' is the file specified as 'genome_file' in 
#' downsample_random()
#' @param add_info set to F by default this boolean should
#' be set to T if the info_file is desired to be combined
#' @param use_weight boolean for whether weights are going
#' to be used for the likelihood calculation
#' @param signame_info the signature that is to be determined
#' with likelihoods the name as in the info file should be 
#' used
#' @param snv_ranges snv ranges to be used for the tuning
#' @param max_allowed_fp max value of false positive rate
#' to be allowed
#' @param method 'weighted_catalog' or 'median_catalog' first

tune_downsampled <- function(directory, 
                             info_file, 
                             add_info = T, 
                             signame_info = "Signature_3", 
                             snv_ranges = c(3, 100, 200, 300, 400, 500),
                             max_allowed_fp = 0.32,
                             method = 'all', 
                             exome = F
                            ){
  source('SOddIUM/R/run.R')
  source('plot_sn_fp.R')
  source('SOddIUM/R/match_to_catalog.R')

  nsnv_ranges <- length(snv_ranges) - 1
  
  for(isnv in 1:nsnv_ranges){
    print(sprintf('%s - %s', snv_ranges[[isnv]], snv_ranges[[isnv + 1]]))

    #define outputs directories and files for plot_sn_fp function
    directory_plot = sprintf('%s/%d_%d_exome%d/plots/', 
                             directory,
                             snv_ranges[[isnv]], 
                             snv_ranges[[isnv + 1]],
                             as.integer(exome))
    directory_output = sprintf('%s/%d_%d_exome%d/output_sn_fp/', 
                               directory,
                               snv_ranges[[isnv]], 
                               snv_ranges[[isnv + 1]],
                               as.integer(exome))
    output_file = sprintf('output_plot_sn_fp_maxfp_%d_method_%s.csv', 
                          as.integer(round(max_allowed_fp * 100, digit = 0)),
                          method)
    #create directories for the primary and secondary outputs
    dir.create(sprintf('%s/%d_%d_exome%d/', 
                       directory, 
                       snv_ranges[[isnv]], 
                       snv_ranges[[isnv + 1]],
                       as.integer(exome)))
    dir.create(directory_plot)
    dir.create(directory_output)
    combine_csv_files(directory = directory, 
                      info_file = info_file, 
                      add_info = add_info, 
                      signame_info = signame_info, 
                      snv_min = snv_ranges[[isnv]], 
                      snv_max = snv_ranges[[isnv + 1]],
                      directory_plot = directory_plot, 
                      directory_output = directory_output,
                      output_file = output_file, 
                      max_allowed_fp = max_allowed_fp,
                      method = method)
  
   
    tune_output <- read.csv(sprintf('%s/tuned_cuts_vs_nsnvs_%s', directory_output, output_file))
    sn_fp_output <- read.csv(sprintf('%s/%s_%s', directory_output, 'total_snvs', output_file))
  
    if(exists('merged_tune')) merged_tune <- rbind(merged_tune, tune_output)
    else merged_tune <- tune_output

    if(exists('merged_sn_fp')) merged_sn_fp <- rbind(merged_sn_fp, sn_fp_output)
    else merged_sn_fp <- sn_fp_output
  }
  
  write.table(merged_tune, 'merged_560bc_tune.csv', row.names = F, quote = F, sep = ",")
  write.table(merged_sn_fp, 'merged_560bc_sn_fp.csv', row.names = F, quote = F, sep = ",")
}

combine_csv_files <- function(directory, info_file, add_info = F, signame_info = 'Sig3', snv_min = 0, snv_max = 500, directory_plot, directory_output, output_file, max_allowed_fp, method){
  load('SOddIUM/R/sysdata.rda')
   
  if(add_info){
     info <- read.csv(info_file)
     info <- info[, grep('Sig|tumor|tumour', colnames(info))]
  }

  if(exists('merged_df_pos')) rm(merged_df_pos)
  if(exists('merged_df_neg')) rm(merged_df_neg)

  files <- list.files(path = directory, pattern = "*.csv")

  for(file in files){
    tumor <- unlist(strsplit(file, split = '.csv'))
    
    df <- read.csv(sprintf('%s/%s', directory, file))
    total_snvs <- rowSums(df[ ,1:96])

    df <- df[total_snvs >= snv_min & total_snvs < snv_max, ]

    if(add_info) info_this <- info[as.character(info$tumor) == tumor,]
    df$code <- paste0(tumor, 1:dim(df)[[1]])

    if(info_this[signame_info] > 0) sig_exists <- 'pos'
    else sig_exists <- 'neg'
 
    df <- cbind(df, info_this)
  
    write.table(df, 
                sprintf('%s/%d_%d_exome%d/annotated_method_%s_%s_%s', 
                        directory,
                        snv_min,
                        snv_max,  
                        as.integer(exome),
                        method,
                        sig_exists,
                        file), 
                sep = ',', 
                row.names = F, 
                quote = F)

    df <- run(sprintf('%s/%d_%d_exome%d/annotated_method_%s_%s_%s', 
                      directory, 
                      snv_min, 
                      snv_max, 
                      as.integer(exome), 
                      method, 
                      sig_exists, 
                      file),
             sprintf('%s/%d_%d_exome%d/output_method_%s_%s_%s', 
                     directory, 
                     snv_min, 
                     snv_max,
                     as.integer(exome),
                     method,
                     sig_exists, 
                     file),
             method = method,
             exome = exome)
 
    if(sig_exists == "pos"){
      if(exists('merged_df_pos')) merged_df_pos <- rbind(df, merged_df_pos)
      else merged_df_pos <- df    
    }
    if(sig_exists == "neg"){
      if(exists('merged_df_neg')) merged_df_neg <- rbind(df, merged_df_neg)
      else merged_df_neg <- df    
    }
  }

  dir.create(sprintf('%s/%d_%d_exome%d/output_pos/', 
                     directory, 
                     snv_min, 
                     snv_max,
                     as.integer(exome)))

  dir.create(sprintf('%s/%d_%d_exome%d/output_neg/', 
                     directory, 
                     snv_min, 
                     snv_max,
                     as.integer(exome)))
  file_pos = sprintf('%s/%d_%d_exome%d/output_pos/merged_method_%s_pos.csv', 
                     directory, 
                     snv_min, 
                     snv_max,
                     as.integer(exome),
                     method)
  file_neg = sprintf('%s/%d_%d_exome%d/output_neg/merged_method_%s_neg.csv', 
                     directory, 
                     snv_min, 
                     snv_max,
                     as.integer(exome),
                     method)

  print(sprintf('writing %s', file_pos))
  write.table(merged_df_pos, 
              file = file_pos, 
              sep = ',', 
              row.names = F,
              quote = F)
  print(sprintf('writing %s', file_neg))
  write.table(merged_df_neg, 
              file = file_neg, 
              sep = ',', 
              row.names = F,
              quote = F)

  if(method == "weighted_catalog") signame1 <- 'Signature_3'
  if(method == "median_catalog") signame1 <- 'Signature_3_c1'

  plot_sn_fp(file_pos,
             file_neg,
             'Signature_3',
             'no_Signature_3',
             dependence = 'total_snvs',
             snv_ranges = c(snv_min, snv_max),
             with_matching = F,
             do_cutoff = T,
             max_allowed_fp = max_allowed_fp,
             plot_dir = directory_plot,
             write_output = T,
             output_file = output_file,
             output_dir = directory_output)

}



