#' Using overlap_repeats() function calculates the indels 
#' overlapping with microsatellites
#'
#' @param input either a maf file path or path to directory containing
#' vcf files
#' @param output_file the output file path where the result will 
#' be saved
#' @param file_type 'maf' or 'vcf'

run_overlap_repeats <- function(input, output_file, file_type, return_df = NULL){

  file_path <- system.file("extdata/repeat_bed_files/",
                           package="SigMA")

  bed_files <- list.files(file_path)
  bed_files <- bed_files[!grepl('tar.gz', bed_files)]  
  for(bed_file in bed_files){

    if(file_type == "maf")
      intersept_this <- overlap_repeats(maf_file = input, repeat_bed_file=paste0(file_path, '/', bed_file))
    else if(file_type == "vcf")
      intersept_this <- overlap_repeats(vcf_dir = input, repeat_bed_file=paste0(file_path, '/', bed_file))
    else
      stop('allowed file_type options: maf, vcf')

    if(exists('intersept'))
      intersept <- rbind(intersept, intersept_this)
    else
      intersept <- intersept_this
    
  }
  intersept$tumor <- as.factor(intersept$tumor)
  nmsi_ins_counts <- by(intersept, intersept$tumor, function(x){sum(x[,'nmsi_ins'])})
  nmsi_del_counts <- by(intersept, intersept$tumor, function(x){sum(x[,'nmsi_del'])})

  nins_counts <- by(intersept_this, intersept_this$tumor, function(x){sum(x[,'nins'])})
  ndel_counts <- by(intersept_this, intersept_this$tumor, function(x){sum(x[,'ndel'])})
  df_nmsi <- data.frame(tumor = unlist(names(nmsi_ins_counts)),
                        nmsi_ins = c(unlist(nmsi_ins_counts)),
                        nmsi_del = c(unlist(nmsi_del_counts)),
                        nins = c(unlist(nins_counts)), 
                        ndel = c(unlist(ndel_counts)))

  if(!is.null(return_df)) return(df_nmsi)
  else
    write.table(df_nmsi, file = output_file, row.names = F, sep = ',', quote = F)
}

