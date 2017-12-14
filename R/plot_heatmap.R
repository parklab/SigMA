#' plots heat map of likelihood (or cosine similarity) for each genome
#' the clustering is made based on the likelihood of each genome to
#' all the signatures
#' 
#' @param inputfile csv file which was obtained by running run() function
#' @param matching 'l' for likelihood and 'c' for cosine similarity
#' @param snv_ranges bins on total number of SNVs used for annotating 
#  the heat map
#' @param file_name the name of the file where the plot will be saved 

plot_heatmap <- function(inputfile, 
                         matching = 'l',
                         snv_ranges = c(5, 25, 45, 65, 
                                        85, 105, 500, 1000),
                         file_name = "heatmap_test.jpg"){

  df <- read.csv(inputfile)
  df$total_snvs <- rowSums(df[,1:96])
  df <- df[df$total_snvs >= snv_ranges[[1]] & df$total_snvs < snv_ranges[[length(snv_ranges)]],]
  df <- na.omit(df)

  #get the columns with likelihood (or cosine similarity) info
  signames <- colnames(df)[grep('Signature_', colnames(df))]
  # choose likelihood or cosine similarity
  signames_m <- signames[grep(paste0('_', matching), signames)]

  # group based on NSNV for the annotation
  group_snv <- apply(df, 1, function(x){
                                 for(isnv in 1:(length(snv_ranges) - 1)){
                                   if(as.numeric(x['total_snvs']) >= snv_ranges[[isnv]] &
                                      as.numeric(x['total_snvs']) < snv_ranges[[isnv + 1]])
                                     return(paste0(snv_ranges[[isnv]], '-',
                                       snv_ranges[[isnv + 1]]))
                                 }
                               })

  
  groups <- data.frame(NSNV = unlist(group_snv))
  rownames(groups) <- 1:dim(groups)[[1]]

  # restrict to the likelihood or cosine similarity values 
  df_for_map <- df[, c(signames_m)]
  df_for_map <- df_for_map * 100
  
  for(icol in 1:dim(df_for_map)[[2]]){
    df_for_map[df_for_map[, icol] < 0.001, icol] <- 0
  }

  #make the gray scale color map based on number of SNV bins
  step_col <- 100/(length(snv_ranges) - 1)  
  groups_col <- rep('', length(snv_ranges) - 1)
  names_groups <- rep('', length(snv_ranges) - 1)
  for(isnv in 1:(length(snv_ranges) - 1)){
    names_groups[[isnv]] <- paste0(snv_ranges[[isnv]], '-', snv_ranges[[isnv + 1]])
    groups_col[[isnv]] <- paste0('gray', round(step_col * (isnv - 1) + 1, digit = 0))
  }
  names(groups_col) <- names_groups
  anno_colors <- list(NSNV = groups_col)

  #construct the transposed df for heat map
  mat <- as.data.frame(t(as.matrix(df_for_map)))
  colnames(mat) <- 1:dim(groups)[[1]]

  #plot the heat map
  jpeg(file_name, width = 1000, height = 1000)
  pheatmap::pheatmap(mat, 
                    annotation = groups, 
                    annotation_colors = anno_colors, 
                    show_colnames = F)
  dev.off()
}