#' This function uses the trained MVA, in particular gradient 
#' boosting  models, inside the package to assign a probability
#' for the existence of the signature of interest. 
#'
#' @param input is a data frame that has likelihood cosine 
#' similarity and total snv values in it's columns
#' @param signame name of the signature which is being identified
#' @param tumor_type tumor type tag see ?run
#' @param data determines the sequencing platform see run()
#'
#' @return a data.frame with a single column with the score
#' of MVA 

predict_mva <- function(input, signame, data, tumor_type, weight_cf, custom){  
  predictions <- rep(0, dim(input)[[1]])
  
  if(data == "msk" & weight_cf){
    model <- gbms_msk_cf[[tumor_type]] 
  }
  else if(data == "wgs_pancan"){
    input$tissue <- tumor_type
    if(tumor_type == "breast") model <- gbm_models[[data]][["breast"]]
    else  model <- gbm_models[[data]][["generic"]]
  } 
  else if(data %in% names(gbm_models)){
    model <- gbm_models[[data]][[tumor_type]] 
  }
  else if(custom){
    file_path <- system.file(paste0("extdata/gbm_models/", data, ".rda"),
                             package="SigMA")
    load(file_path)
    model <- gbm_models_custom[[tumor_type]]
  }
  else stop('the data option selected which indicates the sequencing platform is not valid')

  p_this <- predict(object = model,
                    newdata = input,
                    n.trees = model$n.trees,
                    type = "response")
  predictions <- p_this  
  output <- data.frame(prob = predictions)
  colnames(output)[[1]] <- sprintf('%s_mva', signame)
  return(output)
}
