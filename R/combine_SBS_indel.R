#' combines SBS matrix with indel counts that overlap and
#' that do not overlap with repeat regions
#'
#' @param file_SBS file path containing data frame with 
#' an SBS matrix that can be created with make_matrix() function
#' @param file_overlap_repeat file path containing data frame
#' with indels overlaping with repeat regions
#' @param mut_bed_file mutation bed file that can be created
#' with examples/convert_muts_to_bed.R script
#' @param return_df if set to TRUE merged output data frame
#' is returned. If set to FALSE the merged data table is saved 
#' in the input file, file_SBS, replacing the old file

combine_SBS_indel <- function(file_SBS, file_overlap_repeat, mut_bed_file, return_df = FALSE, msisensor_file = NULL){
  df <- read.csv(file_SBS)
  df <- df[, !(colnames(df) %in% c('nins', 'ndel', 'nmsi_ins', 'nmsi_del'))]
  df_repeat <- read.csv(file_overlap_repeat)
  df_muts <- read.delim(mut_bed_file, header = F)
  df_muts$V4 <- as.factor(df_muts$V4)

  t_ins <- table(df_muts$V4[df_muts$V5 == "INS"])
  t_del <- table(df_muts$V4[df_muts$V5 == "DEL"])

  df_indel <- data.frame(tumor = unlist(names(t_ins)),
                         nins = c(unlist(t_ins)), 
                         ndel = c(unlist(t_del)))

  df <- cbind(df, 
              df_repeat[match(df$tumor, df_repeat$tumor), colnames(df_repeat) != "tumor"],
              df_indel[match(df$tumor, df_indel$tumor), colnames(df_indel) != "tumor"])

  df$nmsi_ins[is.na(df$nmsi_ins)] <- 0
  df$nmsi_del[is.na(df$nmsi_del)] <- 0
  df$nins[is.na(df$nins)] <- 0
  df$ndel[is.na(df$ndel)] <- 0

  if(!is.null(msisensor_file)){
    df_msisensor <- read.csv(msisensor_file)
    df$msisensor <- df_msisensor$msisensor[match(df$tumor, df_msisensor$tumor)]
  }

  if(return_df)  return(df)
  else{
    write.table(df, file = file_SBS, row.names = F, sep = ',', quote = F)
  }
}