#' After running downsample_random(), this code allowes to
#' merge the outputs, adding the initial snv distribution
#' before downsampling and any other info if needed
#' 
#' @param directory is the directory where the outputs of
#' downsample_random() are
#' @param info_file is the file with the tumor column
#' that matches the file names in the 'directory', this 
#' is the file specified as 'genome_file' in 
#' downsample_random()
#' @param add_info set to F by default this boolean should
#' be set to T if the info_file is desired to be combined

combine_csv_files <- function(directory, info_file, add_info = F){
  files <- list.files(directory, pattern = '*.csv')

  if(add_info) info <- read.csv(info_file)

  if(exists('merged_df')) rm(merged_df)

  for(file in files){
    df <- read.csv(sprintf('%s/%s', directory, file))

    tumor <- unlist(strsplit(file, split = '.csv'))
    if(add_info) info_this <- info[as.character(info$tumor) == tumor,]

    df$code <- paste0(tumor, 1:dim(df)[[1]])
    df <- cbind(df, info_this)
   
    if(exists('merged_df')) merged_df <- rbind(df, merged_df)
    else merged_df <- df    
  }

  write.table(merged_df, 
              sprintf('merged_%s.csv', directory), 
              sep = ',', 
              row.names = F,
              quote = F)
}