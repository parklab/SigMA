#' Calculates the indels overlapping with repeat regions.
#' Takes the bed file defining repeat regions and the mutation
#' calls as input and returns a data table with number of indels.
#' Either vcf_dir or maf_file should be provided
#'
#' @param repeat_bed_file bed file containing repeat regions
#' @param output_file output file name in csv format
#' @param vcf_dir directory containing vcf files
#' @param maf_file file path to mutation file in MAF format

overlap_repeats <- function(repeat_bed_file,
                           output_file = NULL,
                           vcf_dir = NULL,
                           maf_file = NULL){

  bed <- read.delim(repeat_bed_file, header = F)
  bed_list <- list()

  for(chr in unique(bed$V1)){
    bed_list[[chr]] <- bed[bed$V1 == chr,]
  }

  tmp <- unlist(strsplit(repeat_bed_file, split = '\\/'))
  bed_name <- tmp[length(tmp)]

  if(!is.null(vcf_dir)){
    tumors <- list.files(vcf_dir)
  }
  else if(!is.null(maf_file)){
    muts <- read.delim(maf_file)
    tumors <- unique(muts$Tumor_Sample_Barcode)
  }
  else{
    error('both MAF and and VCF files are empty')
  }
 
  df <- data.frame(tumor = tumors, count = 0)

  nfile <- length(tumors)
  ifile <- 0

  df_count <- data.frame(nmsi_ins = integer(),
                         nmsi_del = integer(),
                         nins = integer(),
                         ndel = integer(), 
                         tumor = character())
  for(tumor in tumors){
    ifile <- ifile + 1
    if(ifile/nfile %% 100 == 0) print(round(100*ifile/nfile, digit = 0))
    if(!is.null(vcf_dir)){
      vcf <- VariantAnnotation::readVcf(paste0(vcf_dir,'/',tumor))
      gr <- GenomicRanges::granges(vcf)

      chrom_nums <- paste0('chr', as.vector(GenomicRanges::seqnames(gr)))
      seq_start <- GenomicRanges::start(gr)
      seq_end <- GenomicRanges::end(gr)

      ref_vector <- as.character(VariantAnnotation::ref(vcf))
      alt_vector <- as.character(unlist(VariantAnnotation::alt(vcf)))
   
    }
    else if(!is.null(maf_file)){
      this <- muts[muts$Tumor_Sample_Barcode == tumor,]    

      ref_vector <- this$Reference_Allele
      alt_vector <- this$Tumor_Seq_Allele2

      chrom_nums <- paste0('chr', this$Chromosome)
      seq_start <- this$Start_Position
      seq_end <- this$End_Position
    }

    inds_ins <- which(ref_vector == '-' | nchar(ref_vector) < nchar(alt_vector))
    inds_del <- which(alt_vector == '-' | nchar(alt_vector) < nchar(ref_vector))

    if(length(inds_ins) == 0) ins_count <- 0
    else{
      start_ins <- seq_start[inds_ins]
      end_ins <- seq_end[inds_ins]
      chr_ins <- gsub(chrom_nums[inds_ins], pattern = 'chr', replace = '')
    }
    if(length(inds_del) == 0) del_count <- 0
    else{
      start_del <- seq_start[inds_del]
      end_del <- seq_end[inds_del]
      chr_del <- gsub(chrom_nums[inds_del], pattern = 'chr', replace = '')
    }
  
    del_count <- 0
    if(length(inds_del) != 0){
      for(chr in unique(bed$V1)){
        bed_this <- bed_list[[chr]]
        inds <- which(chr_del == chr)
        chr_del_this <- chr_del[inds]
        start_del_this <- start_del[inds]
        end_del_this <- end_del[inds]
        if(length(start_del_this) == 0) next
        for(i in 1:length(start_del_this)){
          bed_this2 <- bed_this[end_del_this[[i]] + 1 >= bed_this$V2 & start_del_this[[i]] - 1 <= bed_this$V3,]
          del <- as.character(bed_this2$V1) == as.character(chr_del_this[[i]]) & 
                ((bed_this2$V2 <= start_del_this[[i]] & bed_this2$V3 >= start_del_this[[i]]) | 
                 (bed_this2$V2 <= end_del_this[[i]]  & bed_this2$V3 >= end_del_this[[i]]) | 
                  bed_this2$V3 == start_del_this[[i]] |
                  bed_this2$V2 == end_del_this[[i]] | 
                  bed_this2$V3 == start_del_this[[i]] - 1 | 
                  bed_this2$V2 == end_del_this[[i]] +1)
          if(sum(del) > 0) del_count <- del_count + 1
        }
      }
    }
    ins_count <- 0
    if(length(inds_ins) != 0){
      for(chr in unique(bed$V1)){
        bed_this <- bed_list[[chr]]
        inds <- which(chr_ins == chr)
        chr_ins_this <- chr_ins[inds]
        start_ins_this <- start_ins[inds]
        if(length(start_ins_this) == 0) next
        end_ins_this <- end_ins[inds]
        for(i in 1:length(start_ins_this)){
          bed_this2 <- bed_this[end_ins_this[[i]] + 1 >= bed_this$V2 & start_ins_this[[i]] - 1 <= bed_this$V3,]

          ins <- as.character(bed_this2$V1) == as.character(chr_ins_this[[i]]) & 
                ((bed_this2$V2 <= start_ins_this[[i]] & bed_this2$V3 >= start_ins_this[[i]]) | 
                 (bed_this2$V2 <= end_ins_this[[i]]  & bed_this2$V3 >= end_ins_this[[i]]) | 
                  bed_this2$V3 == start_ins_this[[i]] |
                  bed_this2$V2 == end_ins_this[[i]] | 
                  bed_this2$V3 == start_ins_this[[i]] - 1 | 
                  bed_this2$V2 == end_ins_this[[i]] +1)
          if(sum(ins) > 0) ins_count <- ins_count + 1
        }
      }
    }
    df_count <- rbind(df_count, data.frame(nmsi_del = del_count, 
                                           nmsi_ins = ins_count, 
                                           ndel = length(inds_del), 
                                           nins = length(inds_ins),
                                           tumor = tumor))
  }
  if(!is.null(output_file)){ 
    write.table(df_count, file = output_file, row.names = F, sep = ',', quote = F)
  }
  return(df_count)
}
