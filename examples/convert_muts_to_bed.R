args = commandArgs(trailingOnly=TRUE)

file_type <- args[[1]]
input <- args[[2]]
output <- args[[3]]

if(file_type == "maf"){
  df_muts <- read.delim(input)
  if(sum(colnames(df_muts) == "Start_Position", na.rm = T) > 0){
    colnames(df_muts)[colnames(df_muts) == "Start_Position"] <- 'start'
  }
  else if(sum(colnames(df_muts) == "Start_position", na.rm = T) > 0){
    colnames(df_muts)[colnames(df_muts) == "Start_position"] <- 'start'
  }
  else{
    stop('input MAF file does not contain Start_Position or Start_position columns')
  }
  if(sum(colnames(df_muts) == "End_Position", na.rm = T) > 0){
    colnames(df_muts)[colnames(df_muts) == "End_Position"] <- 'end'
  }
  else if(sum(colnames(df_muts) == "End_position", na.rm = T) > 0){
    colnames(df_muts)[colnames(df_muts) == "End_position"] <- 'end'
  }
  else{
    stop('input MAF file does not contain End_Position or End_position columns')
  }
  if(sum(colnames(df_muts) %in% c('Variant_Type', 'VARIANT_TYPE', 'VariantType', 'variant_type')) > 0){
    colnames(df_muts)[colnames(df_muts) %in% c('Variant_Type', 'VARIANT_TYPE', 'VariantType', 'variant_type')] <- 'variant_type'
  }
  else{
    stop('input MAF file does not contain variant type in either of the following formats: Variant_Type, VARIANT_TYPE, VariantType, variant_type')
  }
}

if(file_type == "vcf"){
  files <- list.files(input)
  df_muts <- data.frame(chromosome = character(),
             start = integer(),
             end = integer(),
             tumor = character(), 
             variant_type = character())

  for(file in files){
    vcf <- VariantAnnotation::readVcf(paste0(input, '/', file))
    ref_vector <- as.character(VariantAnnotation::ref(vcf))
    alt_vector <- as.character(unlist(VariantAnnotation::alt(vcf)))

    if(length(ref_vector) == 0) next

    gr <- GenomicRanges::granges(vcf)

    seq_start <- GenomicRanges::start(gr)
    seq_end <- GenomicRanges::end(gr)

    chrom_nums <- paste0('chr', as.vector(GenomicRanges::seqnames(gr)))

    variant_type <- rep('SNP', length(seq_start))

    variant_type[which(nchar(ref_vector) < nchar(alt_vector) | ref_vector == "-")] <- 'INS'
    variant_type[which(nchar(ref_vector) > nchar(alt_vector) | alt_vector == "-")] <- 'DEL'
    
    df_muts <- rbind(df_muts,
                    data.frame(Chromosome = chrom_nums,
                              start = seq_start,
                              end = seq_end, 
                              Tumor_Sample_Barcode = file, 
                              variant_type = variant_type))

    
  }
}

df_muts <- df_muts[df_muts$variant_type %in% c('INS','DEL'),]
df_muts$mut_id <- paste0(df_muts$Tumor_Sample_Barcode ,'_', df_muts$Chromosome, '_', df_muts$start , '_', df_muts$end)
write.table(df_muts[, c('Chromosome', 'start', 'end', 'Tumor_Sample_Barcode', 'variant_type', 'mut_id')],
            file = output, row.names = F, sep = '\t', quote = F, col.names = F)
