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

predict_mva <- function(input, signame, data, tumor_type = "breast"){
  input$rat_sig3 <- input$exp_sig3/input$total_snvs
  predictions <- rep(0, dim(input)[[1]])
  
  if(data == "seqcap") model <- gbms_exome[[tumor_type]] #gbm_exome_5 
  else if(data == "msk"){
    model <- gbms_msk[[tumor_type]] 
  }

  else if(data == "wgs"){
    if(tumor_type == "breast") model <- gbms_wgs[["breast"]]
    else{
      inds <- na.omit(match(paste0('Signature_3_c', 1:10, '_ml'), colnames(input)))
      Signature_3_ml <- rowSums(input[,inds])
      model <- gbms_wgs[["generic"]]
    }
  }

  else stop('the data option selected which indicates the sequencing platform is not valid')

  p_this <- predict(object = model,
                    newdata = input,
                    n.trees = 1000,
                    type = "response")
  predictions <- p_this  
  output <- data.frame(prob = predictions)
  colnames(output)[[1]] <- sprintf('%s_mva', signame)
  return(output)
}
