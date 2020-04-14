#' The get_trinuc_norm() can be used to determine the
#' normalization to be used to match the frequency in the 
#' sequencing platform to the frequency in the whole genomes
#'
#' @param bed_file the path to the bed file that defines the 
#' library of the sequencing platform
#' @param do_MC_sampling if set to TRUE, to speed up 
#' the calculation the regions are sampled randomly rather 
# than providing a full count

get_trinuc_norm <- function(bed_file, do_MC_sampling = F,
  ref_genome = BSgenome.Hsapiens.UCSC.hg19:: BSgenome.Hsapiens.UCSC.hg19){
  counts <- get_trinuc_counts(bed_file, do_MC_sampling, ref_genome)

  counts_expanded <- c(rep(counts[1:16], 3), rep(counts[17:32], 3))
  norm <- counts_expanded/counts_trinuc_genome
  return(norm=100*norm/sum(norm))
}

get_trinuc_counts <- function(bed_file, do_MC_sampling = F,
  ref_genome = BSgenome.Hsapiens.UCSC.hg19:: BSgenome.Hsapiens.UCSC.hg19){

  counts_all <- rep(0, 64)
  names(counts_all) <- c('ACA', 'ACC', 'ACG', 'ACT',
                         'CCA', 'CCC', 'CCG', 'CCT',
                         'GCA', 'GCC', 'GCG', 'GCT',
                         'TCA', 'TCC', 'TCG', 'TCT',
                         'TGT', 'GGT', 'CGT', 'AGT',
                         'TGG', 'GGG', 'CGG', 'AGG',
                         'TGC', 'GGC', 'CGC', 'AGC',
                         'TGA', 'GGA', 'CGA', 'AGA',
                         'ATA', 'ATC', 'ATG', 'ATT',
                         'CTA', 'CTC', 'CTG', 'CTT',
                         'GTA', 'GTC', 'GTG', 'GTT',
                         'TTA', 'TTC', 'TTG', 'TTT',
                         'TAT', 'GAT', 'CAT', 'AAT',
                         'TAG', 'GAG', 'CAG', 'AAG',
                         'TAC', 'GAC', 'CAC', 'AAC',
                         'TAA', 'GAA', 'CAA', 'AAA')
                         
  bed <- read.table(bed_file,header=F)
  
  colnames(bed) <- c('chr','start','end','id')[1:dim(bed)[[2]]]
  context_seq_vec <- character()

  if(!grepl('chr', bed$chr[[1]])) bed$chr <- paste0('chr', bed$chr)

  if(do_MC_sampling){
    library_size <- sum(bed$end - bed$start)
    scale_MC_sampling <- (library_size/4000000)
    if(scale_MC_sampling < 1){
      warning('library smaller than 4Mb quick exact calculation feasible')
      do_MC_sampling <- F
    }
  }

  for(i in 1:dim(bed)[[1]]){
    starts <- (bed$start[[i]] - 1):(bed$end[[i]] - 1)
    ends <- (bed$start[[i]] + 1):(bed$end[[i]] + 1)
  
    if(do_MC_sampling){
      inds <- 1:length(starts) 
      n <- length(inds)/scale_MC_sampling
      inds <- unique(round(runif(n, 0.5, inds + 0.5), digit = 0))
      starts <- starts[inds]
      ends <- ends[inds]
    }

    context_seq <- VariantAnnotation::getSeq(ref_genome,
                                           names = bed$chr[[i]],
                                           start = starts,
                                           end = ends,
                                           as.character = TRUE)

    context_seq_vec <- c(context_seq_vec, context_seq)
    counts_context <- table(context_seq)
    inds <- match(names(counts_context), names(counts_all))
    counts_all[na.omit(inds)] <- counts_all[na.omit(inds)] + as.numeric(unlist(counts_context))
  }
  counts <- c(counts_all[1:16] + counts_all[17:32], 
                  counts_all[33:48] + counts_all[49:64]) 
    
  return(counts = counts)
}

