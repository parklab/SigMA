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
  inds_pos <- which(df$pass_mva) 
  pos_tribase <- plot_tribase_dist(as.data.frame(t(df[inds_pos, 1:96])))
  neg_tribase <- plot_tribase_dist(as.data.frame(t(df[-inds_pos, 1:96])))
  

  if(detailed){
    # strip cosine similarity values
    inds <- grep('_c', colnames(df))
    inds_rm <- grep('_c_diff|_ml', colnames(df))
    df_cos <- df[, inds[-na.omit(match(inds_rm, inds))]]
    df_cos <- df_cos[, grep('Signature', colnames(df_cos))]
    print(colnames(df_cos))

    # strip likelihood values
    inds <- grep('_ml', colnames(df))  
    df_ml <- df[, inds]
    df_ml <- df_ml[, grep('Signature', colnames(df_ml))]
    print(colnames(df_ml))
  }
}