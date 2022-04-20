mmej_nhej_indels_from_vcfs <- function(vcf_dir,
                             output_file = 'test.csv',
			     ref_genome_name = 'hg19',
			     min_size_mh = 2,
                             min_size_del = 5,
			     save_detailed_stats = FALSE,
			     snv_matrix_file = NULL){
  
  if(ref_genome_name == "hg19")
    ref_genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  else if(ref_genome_name == "hg38")
    ref_genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  else
    stop('ref_genome_name invalid')
  
  files <- list.files(vcf_dir)
  for(file in files){
    df_this <- mh_from_vcf(vcf_file_path = paste0(vcf_dir, '/', file),
                           ref_genome = ref_genome,
			   min_size_del = min_size_del)
    if(is.null(df_this)) next
    if(exists('df_comb')) df_comb  <- rbind(df_comb, df_this)
    else df_comb <- df_this
  }

  if(save_detailed_stats){
    write.table(df_comb, file = gsub(output_file, pattern = '.csv', replace = '_annotated_muts.csv'), row.names = F, sep = ',', quote = F)
  }
  t_mh <- table(df_comb$file[df_comb$nmh >= min_size_mh & df_comb$nmh != nchar(df_comb$ref)])
  t_nmh <- table(df_comb$file[df_comb$nmh < min_size_mh & nchar(df_comb$ref) - nchar(df_comb$alt) + as.integer(df_comb$alt == "-") >= min_size_del])
    
  df_mh  <- data.frame(file = names(unlist(t_mh)),
                       nmh = c(unlist(t_mh)))
	      
  df_mh$tumor <- unlist(lapply(strsplit(as.character(df_mh$file), split = '\\/'), function(x){x[[length(x)]]}))
  df_nmh <- data.frame(file = names(unlist(t_nmh)),
                       nhej = c(unlist(t_nmh)))

  df_nmh$tumor <- unlist(lapply(strsplit(as.character(df_nmh$file), split = '\\/'), function(x){x[[length(x)]]}))

  df_summary_comb <- data.frame(file = files)
  df_summary_comb$tumor <- unlist(lapply(strsplit(as.character(df_summary_comb$file), split = '\\/'), function(x){x[[length(x)]]}))

  df_summary_comb$mmej <- df_mh$nmh[match(df_summary_comb$tumor, df_mh$tumor)]
  df_summary_comb$mmej[is.na(df_summary_comb$mmej)] <- 0
  df_summary_comb$nhej <- df_nmh$nhej[match(df_summary_comb$tumor, df_nmh$tumor)]
  df_summary_comb$nhej[is.na(df_summary_comb$nhej)] <- 0

  if(!is.null(snv_matrix_file)){
    df_snv <- read.csv(snv_matrix_file)
    df_snv$tumor <- unlist(lapply(strsplit(as.character(df_snv$tumor), split = '\\/'), function(x){ x[[length(x)]]}))
    df_snv$tumor <- gsub(df_snv$tumor, pattern = '.vcf', replace = '')

    df_summary_comb$tumor <- gsub(df_summary_comb$tumor, pattern = '.vcf', replace = '')

    df_snv$mmej <- df_summary_comb$mmej[match(df_snv$tumor, df_summary_comb$tumor)]
    df_snv$nhej <- df_summary_comb$nhej[match(df_snv$tumor, df_summary_comb$tumor)]
    write.table(df_snv, file = output_file, row.names = F, sep = ',', quote = F)
  }
  else{
    write.table(df_summary_comb,  file = output_file, row.names = F, sep = ',', quote = F)
  }
}


mh_from_vcf <- function(vcf_file_path = 'test_vcf/BL-18-F43588_L_2-96252_27-6131_George_PDX_HC_2-96965_676_nalt1.vcf', 
                        output_file = NULL,
                        ref_genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                        min_size_del = 5){
  vcf <- VariantAnnotation::readVcf(vcf_file_path)
  
  chrom <- GenomicRanges::seqnames(GenomicRanges::granges(vcf))
  chrom <- gsub(chrom, pattern ='chr', replace = '')
  inds <- which(chrom %in% c(1:22, 'X', 'Y'))
  vcf <- vcf[inds]
  inds <- which(unlist(lapply(VariantAnnotation::alt(vcf), function(x){ length(x)})) == 1)
  vcf <- vcf[inds]
  if(dim(vcf) == 0) return(NULL) 
  chrom <- GenomicRanges::seqnames(GenomicRanges::granges(vcf))
  
  refs <- as.character(VariantAnnotation::ref(vcf))
  alts <- as.character(unlist(VariantAnnotation::alt(vcf)))
  starts <- GenomicRanges::start(GenomicRanges::ranges(GenomicRanges::granges(vcf)))
  ends <- GenomicRanges::end(GenomicRanges::ranges(GenomicRanges::granges(vcf)))

  df_out <- data.frame()
  if(sum(alts == '-' & nchar(refs) >= min_size_del, na.rm = T) > 0){
    inds <-  which(alts == '-' & nchar(refs) >= min_size_del)
    vcf_this <- vcf[inds,]
    df_this <- calculate_mh(start = starts[inds], 
                            end = ends[inds], 
                            chr = chrom[inds],
                            ref = refs[inds], 
                            alt = alts[inds],
                            ref_genome = ref_genome)
    df_this$file <- vcf_file_path
    df_out <- rbind(df_out, df_this)
  }

  if(sum(alts != '-' & nchar(refs) - nchar(alts) >= min_size_del, na.rm = T) > 0){
    inds <- which(alts != '-' & nchar(refs) - nchar(alts) >= min_size_del & nchar(alts) == 1)
    if(length(inds) == 0) next
    ref_first <- unlist(lapply(VariantAnnotation::ref(vcf[inds,]), function(x){substr(as.character(x), start =1 , stop =1)}))
    ref_rest <- unlist(lapply(VariantAnnotation::ref(vcf[inds,]), function(x){substr(as.character(x), start =2 , stop =length(x))}))

    inds2 <- which(ref_first == alts[inds])
    inds_tmp <- inds[inds2]
    ref_tmp <- ref_rest[inds2]
    chrom_tmp <- chrom[inds_tmp]
    start_tmp <- starts[inds_tmp] + 1
    end_tmp <- ends[inds_tmp] 
    alt_tmp <- '-'

    if(length(inds2) != length(inds)){

      ref_last <- unlist(lapply(ref(vcf[inds,]), function(x){substr(as.character(x), start =length(x) , stop =length(x))}))
      ref_rest2 <- unlist(lapply(ref(vcf[inds,]), function(x){substr(as.character(x), start =1 , stop =length(x) - 1)}))
      inds2 <- which(ref_last == alts[inds])
      if(length(inds2) != length(inds)){
        inds_tmp2 <- inds[inds2]
      }
      inds_tmp2 <- inds_tmp2[is.na(match(inds_tmp2, inds_tmp))]

      if(length(inds_tmp2) > 0){
        ref_tmp2 <- ref_rest2[ind2]
        alt_tmp2 <- '-'
        chrom_tmp2 <- chrom[inds_tmp2]
        start_tmp2 <- starts[inds_tmp2] 
        end_tmp2 <- ends[inds_tmp2] - 1
   
        ref_tmp <- c(ref_tmp, ref_tmp2)
        alt_tmp <- c(alt_tmp, alt_tmp2)
        chrom_tmp <- c(chrom_tmp, chrom_tmp2)
        start_tmp <- c(start_tmp, start_tmp2)
        end_tmp <- c(end_tmp, end_tmp2)
      }
    }

    df_this <- calculate_mh(start = start_tmp, 
                            end = end_tmp,
                            chr = chrom_tmp,
                            ref = ref_tmp,
                            alt = alt_tmp, 
                            ref_genome = ref_genome)
   
    df_this$file <- vcf_file_path
    
    df_out <- rbind(df_out, df_this)
  }
  if(!is.null(output_file))
    write.table(df_out, output_file, row.names = F, sep = ',', quote = F)
  return(df_out)
}
  
calculate_mh <- function(start, end, chr, ref, alt, ref_genome){

  nmh_vec <- rep(0, length(start))
  if(length(start) == 0){
    df_nmh <- data.frame(chrom = character(), start =integer(), end = integer(),
                         ref = character(), alt = character(), nmh = integer()) # add the tumor and alt
  }
  else{
    for(irow in 1:length(start)){ 
      start_this <- start[irow]
      end_this <- end[irow]
      chr_this <- chr[irow]
      ref_this <- ref[irow]

      l_del <- end_this - start_this + 1

      for(il in 1:l_del){
        ncontext = il
        gr_lhs_start <- GenomicRanges::GRanges(chr_this, IRanges::IRanges(start_this - 1, start_this - 1))
        gr_rhs_start <- GenomicRanges::GRanges(chr_this, IRanges::IRanges(start_this , start_this))
        gr_lhs_end <- GenomicRanges::GRanges(chr_this, IRanges::IRanges(end_this , end_this))
        gr_rhs_end <- GenomicRanges::GRanges(chr_this, IRanges::IRanges(end_this + 1, end_this + 1))

        context_lhs_start <- GenomicRanges::resize(gr_lhs_start, ncontext, fix = 'end') 
        context_rhs_start <- GenomicRanges::resize(gr_rhs_start, ncontext, fix = 'start')
        context_lhs_end <- GenomicRanges::resize(gr_lhs_end, ncontext, fix = 'end')
        context_rhs_end <- GenomicRanges::resize(gr_rhs_end, ncontext, fix = 'start')
            
        first_lhs_start <- GenomicRanges::start(context_lhs_start)
        last_lhs_start <- GenomicRanges::end(context_lhs_start)
        first_lhs_end <- GenomicRanges::start(context_rhs_end)
        last_lhs_end <- GenomicRanges::end(context_rhs_end)


        chrom_nums_start <- paste0('chr', as.vector(GenomicRanges::seqnames(context_lhs_start)))
        chrom_nums_end <- paste0('chr', as.vector(GenomicRanges::seqnames(context_lhs_end)))


        context_lhs_start <- VariantAnnotation::getSeq(ref_genome,
                                                 names = chrom_nums_start,
                                                 start = first_lhs_start,
                                                 end = last_lhs_start,
                                                 as.character = TRUE)

        context_lhs_end <- VariantAnnotation::getSeq(ref_genome,
                                                 names = chrom_nums_end,
                                                 start = first_lhs_end,
                                                 end = last_lhs_end,
                                                 as.character = TRUE)


        first_rhs_start <- GenomicRanges::start(context_rhs_start)
        last_rhs_start <- GenomicRanges::end(context_rhs_start)
        first_rhs_end <- GenomicRanges::start(context_rhs_end)
        last_rhs_end <- GenomicRanges::end(context_rhs_end)
        chrom_nums_start <- paste0('chr', as.vector(GenomicRanges::seqnames(context_rhs_start)))
        chrom_nums_end <- paste0('chr', as.vector(GenomicRanges::seqnames(context_rhs_end)))

        context_rhs_start <- VariantAnnotation::getSeq(ref_genome,
                                                 names = chrom_nums_start,
                                                 start = first_rhs_start,
                                                 end = last_rhs_start,
                                                 as.character = TRUE)

      context_rhs_end <- VariantAnnotation::getSeq(ref_genome,
                                               names = chrom_nums_end,
                                               start = first_rhs_end,
                                               end = last_rhs_end,
                                               as.character = TRUE)

 

        if(context_lhs_start != context_lhs_end & context_rhs_start != context_rhs_end){
          nmh <- il - 1
          break
        }else{
          nmh <- il 
        }
      }
      nmh_vec[[irow]] <- nmh
    }
    df_nmh <- data.frame(chrom = chr, start =start, end = end, ref = ref, alt = alt, nmh = nmh_vec) # add the tumor and alt
  }

  return(df_nmh)
}

