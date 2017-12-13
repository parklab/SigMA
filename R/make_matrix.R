#' Converts an input file list with vcf or maf file paths 
#' to a matrix, it works for general number of context
#' and for 1 or 2 strands
#' 
#' @param directory pointer to the directory where input vcf 
#' maf files reside
#' @param file_type 'maf' or 'vcf'
#' @param ref_genome name of the BSgenome currently set by default to
#' BSgenome.Hsapiens.UCSC.hg19
#' @param ncontext number of bases in the nucleotide sequence which 
#' makes up the spectrum
#' @param nstrand number of strands to be considered, 1 contracts to a single 
#' strand which for ncontext = 3 gives the commonly used 96 dimensions

#' @examples
#' by default runs on vcf input and produces 96 dimensional spectra 
#' make_matrix(directory = 'input')
#' make_matrix(directory = 'input', 
#'             file_type = 'vcf',
#'             ref_genome = BSgenome.Hsapiens.UCSC.hg19,
#'             ncontext = 5,
#'             nstrand = 2)

make_matrix <- function(directory, file_type = 'vcf', ref_genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, ncontext = 3, nstrand = 1){
  file_list <- list.files(directory)
  file_list <- paste0(directory, '/', file_list)
 
  if(file_type == "vcf") .make_matrix_from_vcf_list(file_list, ref_genome, ncontext, nstrand)
  else if(file_type == "maf") .make_matrix_from_maf_list(file_list, ref_genome, ncontext, nstrand)
  else stop('file_type should either be maf or vcf')
}

conv_snv_matrix_to_df <- function(genomes_matrix){
  genomes <- as.data.frame(t(genomes_matrix))
  colnames(genomes) <- rownames(genomes_matrix)
  genomes$tumor <- colnames(genomes_matrix)
  return(genomes)
}


##flips a single base
.flip_base <- function(base){
  switch(base, 
  a = {
    return('t')
  },
  c = {
    return('g')    
  },
  g = {
    return('c')
  },
  t = {
    return('a')
  },{
    stop(sprintf('invalid base input %s', base))
  })
}

##flips the bases in a combination to the opposite strand
.flip_strand <- function(comb_in){
  comb_out <- rep('', length(comb_in))
  for(isnv in 1:length(comb_in)){
    snv <- comb_in[[isnv]]
    ref <- substr(snv, 1, 1)
    if( ref == 'c' || ref == 't'){
      comb_out[[isnv]] <- snv
      next
    }else{
      ref <- .flip_base(ref)
      alt <- .flip_base(substr(snv, 2, 2))
      prime5 <- .flip_base(substr(snv, 3, 3))
      prime3 <- .flip_base(substr(snv, 4, 4))
      comb_out[[isnv]] <- paste0(ref, alt, prime3, prime5)
    }
  }
  return(comb_out)
}

##given nmer context, ref, alt information calculates 
##counts for each type of snv
.convert_seq_to_vector <- function(context, ref_vector, alt_vector, types, nstrand = 1){
  context <- tolower(context)
  ref_vector <- tolower(ref_vector)
  alt_vector <- tolower(alt_vector)
  nsnv <- length(context)

  if(nsnv != length(ref_vector) || nsnv != length(alt_vector))
    stop('dimensions of context ref alt are different')
  if(nsnv == 0) 
    stop('empty snv array')
  
 
  range <- nchar(context[[1]])
  if(range %% 2 == 0)
    stop('context should be an odd number')
  

  combined <- paste0(ref_vector, alt_vector) 
  combined <- paste0(combined, substr(context, 1, floor(range/2)))
  combined <- paste0(combined, substr(context, ceiling(range/2)+1, range))

  if(nstrand == 1) combined <- .flip_strand(combined)
  

  count_vector <- rep(0, length(types))
 
  #binary search for counting the occurances of each type

  for(snv in combined){
    L <- 1
    R <- length(types) 
    m <- 1
    while(L <= R){
      m <- floor((L + R)/2)
      if(types[[m]] < snv) L <- m+1
      else if (types[[m]] > snv) R <- m-1
      else{
        count_vector[[m]] <- count_vector[[m]] + 1
        break
      }
    }
  }
  return(count_vector)
}

##given nmer size and strand choice returns 
## an array of types of snvs
.make_type <- function(ncontext = 3, 
                       nstrand = 1){
  components <- c('a', 'c', 'g', 't')
  base_in <- ''
  if(nstrand == 1) base_in <- c('c', 't')
  else base_in <- components  
  types <- rep('', 6*4^(ncontext -1)*nstrand )
  index <- 1
  for(base in base_in){
    for(alt in components[!grepl(base,components)]){
      for(prime5 in components){
        for(prime3 in components){     
          types[[index]] <- paste0(base, alt, prime5, prime3)
           index <- index + 1
        }
      }
    }
  }
  types <- sort(types)
}

.make_vector_from_gr <- function(gr, 
                                 ref_vector, 
                                 alt_vector, 
                                 ref_genome, 
                                 types, 
                                 ncontext = 3, 
                                 nstrand = 1){
  #get a range of ncontext around snv
  gr_context <- GenomicRanges::resize(gr, ncontext, fix = 'center')

  #read the context around snv from the reference
  seq_start <- GenomicRanges::start(gr_context)
  seq_end <- GenomicRanges::end(gr_context)
  
  chrom_nums <- paste0('chr',as.vector(GenomicRanges::seqnames(gr_context)))
  context_seq <- VariantAnnotation::getSeq(ref_genome, 
                        names = chrom_nums, 
                        start = seq_start, 
                        end = seq_end, 
                        as.character = TRUE)

  #convert snv arrays into counts specified by types
  count_vector <- .convert_seq_to_vector(context_seq, 
                                         ref_vector, 
                                         alt_vector, 
                                         types, 
                                         nstrand)
  return(count_vector)
}

##converts a single vcf to a vector of counts of types
.make_vector_from_vcf <- function(vcf_file, 
                                  ref_genome,
                                  types, 
                                  ncontext = 3, 
                                  nstrand = 1){
  #get vcf obj using genomic ranges
  vcf <- VariantAnnotation::readVcf(vcf_file)
  if(dim(vcf)[[1]] == 0)
    return(rep(0, length(types)))

  gr <- GenomicRanges::granges(vcf)  

  #get ref and alt
  ref_vector <- as.character(VariantAnnotation::ref(vcf))
  alt_vector <- as.character(unlist(VariantAnnotation::alt(vcf)))
 
  count_vector <- .make_vector_from_gr(gr, 
                                       ref_vector, 
                                       alt_vector, 
                                       ref_genome, 
                                       types, 
                                       ncontext, 
                                       nstrand)
  return(count_vector)
} 

.make_vector_from_maf <- function(maf_file, 
                                  ref_genome, 
                                  types, 
                                  ncontext = 3, 
                                  nstrand = 1){
  maf <- read.delim(maf_file, 
                    comment.char = "#", 
                    sep = "\t", 
                    header = TRUE, 
                    stringsAsFactors = TRUE)

  if(dim(maf)[[1]] == 0) return(rep(0, 96))

  maf <- maf[maf$Variant_Type == "SNP",]
  maf$Tumor_Seq_Allele2[maf$Tumor_Seq_Allele2 == "TRUE"] <- "T"
  maf <-  maf[maf$Chromosome != "MT",]

  gr <- with(maf, GenomicRanges::GRanges(Chromosome, 
                                         IRanges::IRanges(Start_Position, End_Position)))
  
  ref_vector <- maf$Reference_Allele
  alt_vector <- maf$Tumor_Seq_Allele2

  count_vector <- .make_vector_from_gr(gr, 
                                       ref_vector, 
                                       alt_vector, 
                                       ref_genome, 
                                       types, 
                                       ncontext, 
                                       nstrand)
  return(count_vector)
}

##given a list of vcf files returns a matrix
.make_matrix_from_vcf_list <- function(vcf_files, 
                                       ref_genome, 
                                       ncontext = 3, 
                                       nstrand = 1){
  matrix_snvs <- matrix(0, 6*4^(ncontext - 1)*nstrand, length(vcf_files))
  types <- .make_type(ncontext, nstrand)
  for(ifile in 1:length(vcf_files)){
    print(vcf_files[[ifile]])
    count_vector <- .make_vector_from_vcf(vcf_files[[ifile]], 
                                          ref_genome, 
                                          types, 
                                          ncontext, 
                                          nstrand)
    matrix_snvs[, ifile] <- count_vector
  }  
  rownames(matrix_snvs) <- types
  colnames(matrix_snvs) <- vcf_files
  return(matrix_snvs)
}

.make_matrix_from_maf_list <- function(maf_files, 
                                       ref_genome, 
                                       ncontext = 3, 
                                       nstrand = 1){
  matrix_snvs <- matrix(0, 6*4^(ncontext - 1)*nstrand, length(maf_files))

  types <- .make_type(ncontext, nstrand)
  for(ifile in 1:length(maf_files)){
    count_vector <- .make_vector_from_maf(maf_files[[ifile]], ref_genome, types, ncontext, nstrand)   
    matrix_snvs[, ifile] <- count_vector
  }  
  rownames(matrix_snvs) <- types
  colnames(matrix_snvs) <- maf_files

  return(matrix_snvs)
}

