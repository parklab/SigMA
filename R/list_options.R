#' list_data_options() and list_tumor_types() functions
#' to see information about available settings

list_data_options <- function(){
  for(i in 1:length(platform_names)){
    message(paste0('data option \'', names(platform_names)[[i]], '\' for ', platform_names[[i]]))
  }
}

list_tumor_types <- function(){
  data_dir <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
  load(data_dir)
  for(i in 1:length(tissue_names)){
    message(paste0('tumor_type option \'', names(tissue_names)[[i]], '\' for ', tissue_names[[i]]))  
    if(!is.na(match(names(tissue_names)[[i]], names(matrices_96dim[['wgs']])))){
      printout <- 'matrices available for Sig3 tuning'
      names_wgs <- names(gbms_wgs_tumor_type)
      names_tcga <- names(gbms_mc3)
      names_seqcap <- names(gbms_exome)
      names_seqcap_probe <- names(gbms_seqcap_probe)
      names_msk <- names(gbms_msk)
      models_avail <- character()
      if(!is.na(match(names(tissue_names)[[i]], names_wgs))) 
        models_avail <- c(models_avail, 'wgs')
      if(!is.na(match(names(tissue_names)[[i]], names_tcga))) 
        models_avail <- c(models_avail, 'tcga_mc3')
      if(!is.na(match(names(tissue_names)[[i]], names_seqcap))) 
        models_avail <- c(models_avail, 'seqcap')
      if(!is.na(match(names(tissue_names)[[i]], names_seqcap_probe)))
        models_avail <- c(models_avail, 'seqcap_probe')
      if(!is.na(match(names(tissue_names)[[i]], names_msk)))
        models_avail <- c(models_avail, 'msk')
      message(paste0('models available for ', paste0(models_avail, collapse = ', '), '\n'))
    }
    else message('\n')
  }
}

info <- function(data, tumor_type){
  stat_this <- stat[stat$tumor_type == tumor_type,]
  info <- character()
  for(i in 1:dim(stat_this)[[1]]){
    info <- c(info, paste0(' in ', stat_this$data[i], ' is ', round(stat_this$Sig3_frac[i], digit = 2)))
  }
  message(paste0('Sig3 fraction', paste0(info, collapse = ' and')))
  message(paste0('Median SNV count:', 
   snv_counts$median_total_snvs[snv_counts$data == data & snv_counts$tumor_type == tumor_type]))
}