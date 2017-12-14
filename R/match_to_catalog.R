#' Calculates the compatibility of a list of genomes
#' to an input catalog based on likelihood and cosine 
#' similarity
#' 
#' @param genomes a data table or matrix with snv  
#' spectra in the first ntype columns and genomes in each row 
#' @param signatures the input catalog, a data table
#' with signature spectra in each column 
#' @return A data frame that contains the input genomes
#' and in addition columns associated to each signature in
#' in the catalog with likelihood and cosine simil values

match_to_catalog <- function(genomes, signatures, ntype = 96, use_weight = F){

  # if the dataset is bc this part is using weights for each signature to
  # pick the most distinct context for calculating the probabilities
  if(use_weight){
    ind_bc <- grep('Signature_', colnames(weights_560_bc_cooccur_PCAWG_sig))
    signatures <- signatures[, colnames(weights_560_bc_cooccur_PCAWG_sig)[ind_bc]]
    weights <- weights_560_bc_cooccur_PCAWG_sig[, ind_bc]
    signatures <- signatures * weights 
  }else{
    weights <- matrix(1, dim(weights_560_bc_cooccur_PCAWG_sig)[[1]],
                      dim(weights_560_bc_cooccur_PCAWG_sig)[[2]])
    colnames(weights) <- colnames(weights_560_bc_cooccur_PCAWG_sig)
  }  

  calc_prob <- function(this_genome, signatures){
    eps <- 0.00001

    if(sum(this_genome) == 0){
      return(list(probs = rep(0, dim(signatures)[[2]]), 
                  ind_max = NA, 
                  sig_max = NA, 
                  max_val = NA))
    }

    probs <- apply(signatures, 2, 
	           function(x, pow){ 
                     prob <- 0
                     for(i in 1:length(x)){
                       if(x[[i]] == 0) x[[i]] <- eps
	               if(pow[[i]] > 0){
                         prob <- prob + pow[[i]] * log(x[[i]]) 
                       }
	             }
	             return(prob)
                   },
                   pow = this_genome)


    mean_probs <- mean(probs)
    q3_probs <- mean(probs[probs > mean(probs)])
    max_probs <- max(probs)

    inds_keep <- which(probs >= q3_probs)
    inds_rm <- which(probs < q3_probs)
    
    probs[inds_keep] <- exp(probs[inds_keep] - max_probs)
    probs[inds_rm] <- 0

    ind_max <- which(max(probs) == probs)[[1]]

    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- probs[[ind_max]]
    return(list(probs = probs, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
  }

  cosine <- function(x, y){
    return(x%*%y / sqrt(x%*%x * y%*%y))
  }

  calc_cos <- function(this_genome, signatures){
    if(sum(this_genome) == 0){
      return(list(simils = rep(0, dim(signatures)[[2]]), 
                  ind_max = NA, 
                  sig_max = NA, 
                  max_val = NA))
    }

    simils <- apply(signatures, 2,
                    function(x){ cosine(this_genome, x) })
    ind_max <- which(max(simils) == simils)
    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- simils[[ind_max]]
    return(list(simils = simils, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
  }

  signature_names <- colnames(signatures)
  genome_matrix <- genomes[, 1:dim(signatures)[[1]]]
  
  probs_all    <- apply(genome_matrix, 1, 
                        function(x, y){calc_prob(x, y)$probs}, y = signatures)
  max_sigs_all <- apply(genome_matrix, 1, 
                        function(x, y){calc_prob(x, y)$sig_max}, y = signatures)
  max_val_all  <- apply(genome_matrix, 1, 
                        function(x, y){calc_prob(x, y)$max_val}, y = signatures)

  simils_all    <- apply(genome_matrix, 1, 
                        function(x, y){calc_cos(x, y)$simils}, y = signatures)
  max_sigs_cos_all <- apply(genome_matrix, 1, 
                        function(x, y){calc_cos(x, y)$sig_max}, y = signatures)
  max_val_cos_all  <- apply(genome_matrix, 1, 
                        function(x, y){calc_cos(x, y)$max_val}, y = signatures)

  colnames_before <- colnames(genomes)
  genomes <- as.data.frame(genomes)

  genomes <- cbind(genomes, t(probs_all), max_sigs_all, max_val_all, t(simils_all), max_sigs_cos_all, max_val_cos_all)
  colnames(genomes) <- c(colnames_before, paste0(colnames(signatures), '_l'), 'sig_max_l', 'max_l', 
                      paste0(colnames(signatures), '_c'), 'sig_max_c', 'max_c')

  return(genomes)
}