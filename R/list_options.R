#' list_data_options() and list_tumor_types() functions
#' to see information about available settings

list_data_options <- function(){
  for(i in 1:length(platform_names)){
    message(paste0('data option \'', names(platform_names)[[i]], '\' for ', platform_names[[i]]))
  }
}

list_tumor_types <- function(){
  matrices_data <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
  load(matrices_data)
  for (i in 1:length(tissue_names)){
    tumor_type <- tissue_names[[i]]
    short_name <- names(tissue_names)[[i]]
    message(paste0('tumor_type option \'', short_name, '\' for ', tumor_type))

    models_avail <- names(gbm_models)[sapply(gbm_models, function(x) short_name %in% names(x))]
    if (length(models_avail)>0){
        message(sprintf('models available for %s\n', paste0(models_avail, collapse = ', ')))
    } else {
        message('no models available\n')
    }
  }
}

info_stat <- function(data, tumor_type){
  stat_this <- stat[stat$tumor_type == tumor_type,]
  info <- character()
  for(i in 1:dim(stat_this)[[1]]){
    info <- c(info, paste0(' in ', stat_this$data[i], ' is ', round(stat_this$Sig3_frac[i], digit = 2)))
  }
  message(paste0('Sig3 fraction', paste0(info, collapse = ' and')))
  message(paste0('Median SNV count:', 
   snv_counts$median_total_snvs[snv_counts$data == data & snv_counts$tumor_type == tumor_type]))
}

has_model <- function(data, tumor_type){
  if(!is.na(match(data, names(platform_names)))){
    tumor_types <- names(gbm_models[[data]]) 
    if(tumor_type %in% tumor_types) return(TRUE)
    else return(FALSE)
  }
  else{
    file_path <- system.file("extdata/gbm_models/",
                             package="SigMA")
    filename <- paste0(file_path, data, '.rda')
    if(file.exists(filename)){
      load(filename)
      tumor_types <- names(gbm_model)
      if(!is.na(match(tumor_type, tumor_types))) return(TRUE)
      else return(FALSE)
    } 
    else{ 
      return(FALSE)
    }
  }
}