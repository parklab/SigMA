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

  input$rat_sig3 <- input$exp_sig3/input$total_snvs
  
  if(length(grep(paste0(paste0('Signature_3_c', 1:10, '_ml'), collapse = '|'), colnames(input))) == 1)
    input$Signature_3_ml <- input$Signature_3_c1_ml
  else
    input$Signature_3_ml <- rowSums(input[, grep(paste0(paste0('Signature_3_c', 1:10, '_ml'), collapse = '|'), colnames(input))])

  predictions <- rep(0, dim(input)[[1]])
  
  if(data == "msk"){
    if(weight_cf)
      model <- gbms_msk_cf[[tumor_type]] 
    else 
      model <- gbms_msk[[tumor_type]] 
  }
  else if(data == "tcga_mc3"){
    model <- gbms_mc3[[tumor_type]]
  }
  else if(data == "seqcap"){
    model <- gbms_exome[[tumor_type]] 
  }
  else if(data == "seqcap_probe"){
    model <- gbms_seqcap_probe[[tumor_type]] 
  }
  else if(data == "wgs"){
    model <- gbms_wgs_tumor_type[[tumor_type]]
  }
  else if(data == "wgs_pancan"){
    input$tissue <- tumor_type
    if(tumor_type == "breast") model <- gbms_wgs[["breast"]]
    else{
      inds <- na.omit(match(paste0('Signature_3_c', 1:10, '_ml'), colnames(input)))
      if(length(inds) > 1) input$Signature_3_ml <- rowSums(input[,inds])
      if(length(inds) == 1) input$Signature_3_ml <- input[,inds]
      model <- gbms_wgs[["generic"]]
    }
  }
  else if(custom){
    file_path <- system.file(paste0("extdata/gbm_models/", data, ".rda"),
                             package="SigMA")
    load(file_path)
    model <- gbm_model
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
