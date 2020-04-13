#' The find_best_model() function determines the best value 
#' for the data parameter in run.R function based on the 
#' median count of SNVs in the dataset that is analyzed. 
#' 
#' @param input_file the path to the file that contains the
#' data table that will be used in run.R function
#' @param tumor_type see list_tumor_types() for more information
#' @param remove_msi_pole can be set to TRUE to remove 
#' microsattelite instable and POLE-exo domain mutated hypermutated
#' tumors because these bias the mutation count. If the data does
#' not contain hypermutated samples it can be set to FALSE to speed
#' up the calculation 

find_data_setting <- function(input_file, tumor_type, best_model = NULL, remove_msi_pole = T){
  if(!is.null(best_model)){
    output_file <- run(input_file, tumor_type = tumor_type, data = best_model, do_mva = F, do_assign = F, check_msi = T)
    df <- read.csv(output_file)
    df$Signature_msi_ml <- rowSums(df[,grepl('Signature_msi_c', colnames(df))])
    df <- df[!((df$Signature_pole_c1_ml_msi > 0.999 | df$Signature_msi_ml > 0.99) & df$total_snvs > median(df$total_snvs, na.rm = T)),]
  }
  else{
    df <- read.csv(input_file)
  }

  if(sum(grepl('total_snvs', colnames(df))) == 0){
    df$total_snvs <- rowSums(df[,1:96])
  }

  mean_this <- mean(df$total_snvs, na.rm = T)
  median_this <- median(df$total_snvs, na.rm = T)

  stat_this <- snv_counts[snv_counts$tumor_type == tumor_type,]

  mean_diff <- abs(mean_this - stat_this$mean_total_snvs)
  median_diff <- abs(median_this - stat_this$median_total_snvs)

  data1 <- stat_this$data[which(mean_diff == min(mean_diff))]
  data2 <- stat_this$data[which(median_diff == min(median_diff))]


  if(is.null(best_model) & remove_msi_pole){
    find_data_setting(input_file, tumor_type, best_model = data1)
  }
  else{
    if(data1 == data2) return(data1)
    else{
      warning(paste0(data1, ' is also feasible'))
      return(data2)
    }
    if(abs(min(median_diff) - median_this)/median_this > 0.25){
      warning('retuning suggested see quick_simulation() or contact author')
    }
  }
}

