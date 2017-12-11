#' This function returns a data table with mutational spectra
#' of MC sample generated based on a single signatures probability
#' distribution function
#' 
#' @param signame signature name to use to be used to generate the 
#' distributions
#' @param nsnv_low lowest mutation count 
#' @param nsnv_high highest mutation count
#' @param distribution of mutation counts set to 'flat by default
#' other user specific distributions will be implemented
#' @param file_name output file name
#' @param sig_catalog whether to use 'pcawg' or 'cosmic' catalogs

generate_single_sig_mc <- function(signame,
                                   nsnv_low = 5,
                                   nsnv_high = 104,
                                   distribution = 'flat',
                                   ngenomes = 1000,
                                   file_name = sprintf('mc_%s_snvs.csv',
                                               signame),
                                   sig_catalog = 'cosmic'){

  signatures_cosmic_file <- system.file("extdata",
                                        "sorted_cosmic_signatures.csv",
                                        package = "LowRStat")

  signatures_pcawg_file <- system.file("extdata",
                                       "sorted_PCAWG_signature_patterns_beta2.csv",
                                       package = "LowRStat")

  if(sig_catalog == "cosmic") signatures <- read.csv(signatures_cosmic_file)
  if(sig_catalog == "pcawg") signatures <- read.csv(signatures_pcawg_file)
 
  dist <- signatures[,signame]
  ntype <- length(dist)
  genomes <- matrix(0, ngenomes, ntype)
  total_snvs <- rep(0, ngenomes)

  if(distribution == 'flat'){
    per_value <- ngenomes/(nsnv_high - nsnv_low + 1)
    per_value <- rep(as.integer(per_value), nsnv_high - nsnv_low + 1)
  }else{
    # write here a function that generates number of snv distribution 
    # based on the distribution parameter
  }

  index <- 1
  for(isnv in nsnv_low:nsnv_high){
    local_nsnv <- per_value[[isnv - nsnv_low + 1]]
    for(ite in 1:local_nsnv){
      total_snvs[index] <- isnv
      while(sum(genomes[index, ]) < isnv){
        location <- as.integer(round(runif(1, 0.5, ntype + 0.5), digit = 0))
        rand_prob <- runif(1, 0, 1)
        if(rand_prob < dist[[location]])
          genomes[index, location] <- genomes[index, location] + 1       
      }
      index <- index + 1       
    }
  } 

  colnames_snv <- paste0(signatures$Somatic_Mutation_Type)

  df <- data.frame(genomes, total_snvs = total_snvs)
  colnames(df)[1:96] <- colnames_snv
  write.table(df, file_name, sep = ',', row.names = F, col.names = T, quote = F)
}