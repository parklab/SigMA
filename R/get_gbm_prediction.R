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

get_gbm_prediction <- function(input, signame, exome, panel){
  predictions <- rep(0, dim(input)[[1]])
  if(exome) model <- gbm_exome
  else if(panel) model <- gbm_msk
  else model <- gbm_wgs

  p_this <- predict(object = model,
                    newdata = input,
                    n.trees = 1000,
                    type = "response")
  predictions <- p_this  
  output <- data.frame(prob = predictions)
  colnames(output)[[1]] <- sprintf('%s_gbm', signame)
  return(output)
}
