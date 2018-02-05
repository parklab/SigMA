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

get_gbm_prediction <- function(input, signame, exome){
  input$total_snvs[input$total_snvs >= 500] <- 499
  boundary_snv <- c(4, 100, 200, 300, 500)
  predictions <- rep(0, dim(input)[[1]])
  for(ibin in 1:(length(boundary_snv) - 1)){
    indices <- which(input$total_snvs >= boundary_snv[[ibin]] &
                     input$total_snvs < boundary_snv[[ibin + 1]])
    input_this <- input[indices, ]
    if(exome) model <- gbm_exome
    else model <- list_gbm[[ibin]]
    p_this <- predict(object = model,
                          newdata = input_this,
                          n.trees = 1000,
                          type = "response")
    predictions[indices] <- p_this  
  }
  output <- data.frame(prob = predictions)
  colnames(output)[[1]] <- sprintf('%s_gbm', signame)

  return(output)
}
