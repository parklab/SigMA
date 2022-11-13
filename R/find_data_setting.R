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

find_data_setting <- function(input_file = NULL, tumor_type, best_model = NULL, remove_msi_pole = T, input_df = NULL, catalog_name = 'cosmic_v3p2_inhouse', with_classifier = F){

  if(!is.null(input_file) & is.null(input_df)){
    if(!is.null(best_model)){
      output_file <- run(input_file, tumor_type = tumor_type, data = best_model, do_mva = F, do_assign = F, check_msi = T, input_df = input_df, catalog_name = catalog_name)
      df <- read.csv(output_file)
      df$Signature_msi_ml <- rowSums(df[,grepl('Signature_msi_c', colnames(df))])
      df <- df[!((df$Signature_pole_c1_ml_msi > 0.999 | df$Signature_msi_ml > 0.99) & df$total_snvs > median(df$total_snvs, na.rm = T)),]
    } 
    else{
      df <- read.csv(input_file)
    }
  }
  else if(is.null(input_file) & !is.null(input_df)){
    df <- input_df
  }
  else if(!is.null(input_file) & !is.null(input_df)){
    stop('both input_file and input_df is provided only one is needed')
  }
  else{
    stop('both input_file and input_df is NULL')
  }


  if(sum(grepl('total_snvs', colnames(df))) == 0){
    df$total_snvs <- rowSums(df[,1:96])
  }

  mean_this <- mean(df$total_snvs, na.rm = T)
  median_this <- median(df$total_snvs, na.rm = T)

  stat_this <- snv_counts[snv_counts$tumor_type == tumor_type,]
  if(with_classifier){
    has_classifier <- unlist(apply(stat_this, 1, function(x, tumor_type){ has_model(data = x['data'], tumor_type = tumor_type)}, tumor_type = tumor_type))
    if(sum(has_classifier, na.rm = T) == 0) stop('no classifier for Sig3 for this tumor type, data setting can still be helpful for purposes so to find the best data setting without a classifier set parameter with_classifier = T')
    else
      stat_this <- stat_this[has_classifier,]
  }
  mean_diff <- abs(mean_this - stat_this$mean_total_snvs)
  median_diff <- abs(median_this - stat_this$median_total_snvs)

  data1 <- stat_this$data[which(mean_diff == min(mean_diff))]
  data2 <- stat_this$data[which(median_diff == min(median_diff))]
  if(length(data1) == 2){
    if(sum(!is.na(match(data1, c('msk', 'msk_with_zeros')))) == 2)
      data1 <- 'msk'
    if(sum(!is.na(match(data2, c('msk', 'msk_with_zeros')))) == 2)
      data2 <- 'msk'
  }
  
  if(data1 == "msk_with_zeros") data1 <- 'msk'
  if(data2 == "msk_with_zeros") data2 <- 'msk'
 
  if(length(data2) > 1){
    if(length(data1) == 1) data2 <- data1
    else data2 <- data2[[1]]
  }

  if(is.null(best_model) & remove_msi_pole){
    find_data_setting(input_file, tumor_type, best_model = data1)
  }
  else{
    if(abs(min(median_diff) - median_this)/median_this > 0.25){
      warning(paste0('retuning suggested because the difference in SNV counts between the existing data options and your data is large, see quick_simulation() or contact author. Median from tune:', median_diff, ', median in data:', median_this))
    }
    if(data1 == data2) return(data1)
    else{
      warning(paste0(data1, ' is also feasible'))
      return(data2)
    }
  }
}

adjust_cutoff <- function(input_df, data, tumor_type, below_cutoff){
  message('adjusting cutoff for MVA score')
  if(!is.null(below_cutoff)){
    message('there are samples below SNV cutoff which are taken into account')
    input_df <- rbind(input_df[,1:96], below_cutoff[,1:96])
    input_df$total_snvs <- rowSums(input_df[,1:96])
  }

  if(sum(names(dynamic_cutoff) == data) == 0) return(NULL)
  else if(sum(names(dynamic_cutoff[[data]]) == tumor_type) == 0) return(NULL)

  median_data <- median(input_df$total_snvs, na.rm = T)
  sd_data <- median(input_df$total_snvs, na.rm = T)
  median_data <- median(input_df$total_snvs[input_df$total_snvs < median_data + 5*sd_data], na.rm = T)
  
  this <- dynamic_cutoff[[data]][[tumor_type]]

  if(sum(this$median_total_snvs == median_data) == 0) return(NULL)

  cutoff_strict = this$cutoff[this$median_total_snvs == median_data & this$limit == 0.01]
  cutoff = this$cutoff[this$median_total_snvs == median_data & this$limit == 0.1]
  
  
  return(list(cutoffs = c(cutoff, cutoff_strict), limits = c(0.1, 0.01), cut_var = 'fpr'))
}