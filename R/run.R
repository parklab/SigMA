#' Loads signatures and genome snv matrices and runs the matching
#'
#' @param genome_file a file with snv spectra info can be created
#' from vcf file using @make_genome_matrix() function:
#' see ?make_genome_matrix
#' @param output_file the output file name
#' @param sig_catalog an array of 'cosmic' 'pcawg' 'default
#' 'custom', 'default option uses the intersection of cosmic
#' and pcawg catalogs adding the subsignatures from pcawg of 
#' existing cosmic signatures, e.g. Signature 7a, 7b, etc, but
#' the new signatures introduced by PCAWG catalog are not used
#' @param custom_sig_file a file with signature distributions
#' @param rm_sigs an array with signature names that the user
#' would like to remove in the decomposition
#'
#' @examples
#' run(genome_file = 'input_genomes.csv', 
#'     sig_catalog = c('cosmic', 'pcawg'))
#' run(genome_file = 'input_genomes.csv', 
#'     sig_catalog = 'custom', 
#'     custom_sig_file = 'signatures_me.csv')
#' run(genome_file = 'input_genomes.csv',
#'     sig_catalog = c('pcawg')
#'     rm_sigs = 'Signature_40')

run <- function(genome_file, 
                output_file = 'output.csv',
                sig_catalog = 'default', 
                custom_sig_file = NULL,
                rm_sigs = NULL){

  genomes <- read.csv(genome_file)
  signatures_cosmic_file <- system.file("extdata", 
                                        "sorted_cosmic_signatures.csv", 
                                        package = "LowRStat")
  signatures_pcawg_file <- system.file("extdata", 
                                       "sorted_PCAWG_signature_patterns_beta2.csv", 
                                       package = "LowRStat")


  # set signature catalog to 'cosmic' or 'pcawg' if only one is provided
  # if both 'cosmic' and 'pcawg' is set then take the union of the two
  # catalogs use 'pcawg' when available
  if(length(grep("pcawg", sig_catalog)) > 0){
    signatures_pcawg <- read.csv(signatures_pcawg_file)
    signatures_pcawg <- signatures_pcawg[, grep('Signature', colnames(signatures_pcawg))]
    signatures <- signatures_pcawg
  }
  if(length(grep("cosmic", sig_catalog)) > 0){
    signatures_cosmic <- read.csv(signatures_cosmic_file)
    signatures_cosmic <- signatures_cosmic[, grep('Signature', colnames(signatures_cosmic))]
    if(!exists('signatures')) signatures <- signatures_cosmic
    else{
      matched_inds <- match(colnames(signatures_pcawg), colnames(signatures_cosmic))
      matched_inds <- matched_inds[!is.na(matched_inds)]
      signatures_cosmic <- signatures_cosmic[, -matched_inds]
      signatures <- cbind(signatures_pcawg, signatures_cosmic)
    }
  }

  # if sig_catalog has 'default' it overwrites other options the intersection of
  # two catalogs rather than the union is used
  if(length(grep("default", sig_catalog)) > 0){
    if(!exists('signatures_pcawg')){
      signatures_pcawg <- read.csv(signatures_pcawg_file)
      signatures_pcawg <- signatures_pcawg[, grep('Signature', colnames(signatures_pcawg))]
    }
    if(!exists('signtures_cosmic')){
      signatures_cosmic <- read.csv(signatures_cosmic_file)
      signatures_cosmic <- signatures_cosmic[, grep('Signature', colnames(signatures_cosmic))]
    }
 
    matched_inds <- match(colnames(signatures_pcawg), colnames(signatures_cosmic))
    matched_inds <- matched_inds[!is.na(matched_inds)]
    signatures_cosmic <- signatures_cosmic[, matched_inds]
    sub_signatures <- c('Signature_7a', 'Signature_7b', 'Signature_17a', 'Signature_17b',
                        'Signature_10a', 'Signature_10b')
    signatures <- cbind(signatures_cosmic, signatures_pcawg[, sub_signatures])
    # signatures_arti <- signatures_pcawg[, grep('Signature_R', colnames(signatures_pcawg))]
    # signatures <- cbind(signatures, signatures_arti)
  }

  # custom overwrites all the options above and uses the user defined input file 
  if(length(grep("custom", sig_catalog) > 0)){
    signatures <- read.csv(custom_sig_file)
  }
 
  genomes <- read.csv(genome_file)
  output <- match_to_catalog(genomes, signatures)
  

  write.table(output, output_file, sep = ',', row.names = F, col.names = T, quote = F)
  return(output)
}