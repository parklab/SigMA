#' This function uses the trained GBM models in the package
#' to assign a probability to each genome that it has the 
#' signature of interest based on the likelihood and cosine 
#' similarity values
#'
#' @param input is a data frame that has likelihood cosine 
#' similarity and total snv values in it's columns
#' @param signame it's only used for naming the column of
#' output file at the moment in the future will be used
#' for identifying a specific signature

get_gbm_prediction <- function(input, signame, data, step = "mss"){
  predictions <- rep(0, dim(input)[[1]])
  if(step == "mss"){
    if(data == "seqcap") model <- gbm_exome_5 #gbm_ccle_sanger #gbm_ccle_v2 #gbm_exome_5
    else if(data == "msk") model <- gbm_msk_5 #gbm_msk_3
    else if(data == "found") model <- gbm_fo_5 #gbm_found_3
    else model <- gbm_wgs_5
  }
  else{
    if(data == "seqcap") model <- gbm_exome_msi
    else if(data == "msk") model <- gbm_msk_msi
    else if(data == "found") model <- gbm_fo_msi
    else model <- gbm_wgs_msi
  }
  p_this <- predict(object = model,
                    newdata = input,
                    n.trees = 1000,
                    type = "response")
  predictions <- p_this  
  output <- data.frame(prob = predictions)
  colnames(output)[[1]] <- sprintf('%s_gbm', signame)
  return(output)
}
