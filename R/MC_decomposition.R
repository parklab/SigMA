#' Associates each SNV to a specific signature based on 
#' probabilities. The signatures and their exposures must
#' be discovered by another algorithm and used as input for
#' this function
#' 
#' @param genomes SNV spectra of the genome.
#' @param exposures exposures of the signatures in the .
#' @param signatures signatures found in the genome
#' @return A matrix of size ntype x ngenome x nsig
#' where ntype: SNV spectra dimensions
#' ngenome: number of genomes
#' nsig: number of signatures 
#' @examples
#' run_decomposition('input_snvs_genome.csv', 
#'                   'signatures.csv',
#'                   'exposures.csv')
#' MC_decomposition(genomes, signatures, exposures)

MC_decomposition <- function(genomes, signatures, exposures){
  ntype <- dim(genomes)[[1]]
  ngenome <- dim(genomes)[[2]]
  nsig <- dim(signatures)[[2]]
  reco <- signatures %*% exposures
 
  MAX <- max(max(genomes), max(reco))
  genome_sig_assoc = array(0, c(ntype, ngenome, nsig))
 
  for(i in 1:ntype){
    for(j in 1:1){
      this_snv_num <- genomes[i,j]
      if(this_snv_num == 0) next
      
      total = 0
      for(isnv in 1:this_snv_num){
        trial = 0
        while(TRUE){
          dice_sig <- floor(nsig*runif(1, 0, 1)) + 1
          number_prob <- signatures[i, dice_sig] * exposures[dice_sig, j]
          dice_num <- runif(1, 0, MAX)
          if(number_prob >= dice_num){
            genome_sig_assoc[i, j, dice_sig] = genome_sig_assoc[i, j, dice_sig] + 1
            total = total + 1
            break
          }
          trial = trial + 1
        }
      }
    }
  }
  save(genome_sig_assoc, file = 'out_genome_sig_asso_genome_by_genome.Rda')
  test <- rep(0,96)
  for(i in 1:nsig){
    test = test + genome_sig_assoc[, 1, i]
  }
}

#run_decomposition<-function(genome_file, sig_file, exp_file){
#  exposures <- read.csv(exp_file)
#  signatures <- read.csv(sig_file)
#  snvs <- read.csv(genome_file)
#  if(dim(exposures)[[1]] != dim(signatures)[[2]]) 
#    stop("dimensions of signature and exposure matrices do not match")
#  if(dim(exposures)[[2]] != dim(genomes)[[2]]) 
#    stop("dimensions of genomes and exposure matrices do not match")
  
#  MC_decomposition(snvs, exposures, signatures)
#}