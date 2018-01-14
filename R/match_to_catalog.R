#' Calculates the compatibility of a list of genomes
#' to an input catalog based on likelihood and cosine 
#' similarity
#' 
#' @param genomes a data table or matrix with snv  
#' spectra in the first ntype columns and genomes in each row 
#' @param signatures the input catalog, a data table
#' with signature spectra in each column 
#' @param method can be 'median_catalog', 'weighted_catalog'
#' 'cosine_simil'. 'median_catalog' uses a custom signature
#' catalog formed by clustering genome SNV spectra and 
#' using it as a probability distribution. The 'median_catalog'
#' method can be used with any custom signatures data frame 
#' if the user intends to provide their own signature table.
#'
#' @return A data frame that contains the input genomes
#' and in addition columns associated to each signature in
#' in the catalog with likelihood and cosine simil values

match_to_catalog <- function(genomes, signatures, method = 'median_catalog'){

  if(method == 'weighted_catalog'){
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
    q1_probs <- mean(probs[probs < mean(probs)])
    max_probs <- max(probs)

    inds_keep <- which(probs >= q1_probs)
    inds_rm <- which(probs < q1_probs)
    
    probs[inds_keep] <- exp(probs[inds_keep] - max_probs)
    probs[inds_keep] <- probs[inds_keep]/sum(probs[inds_keep])

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
    simils <- apply(signatures, 2,
                    function(x){ cosine(this_genome, x) })
    ind_max <- which(max(simils) == simils)
    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- simils[[ind_max]]
    return(list(simils = simils, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
  }

  signature_names <- colnames(signatures)
  genome_matrix <- genomes[, 1:dim(signatures)[[1]]]
  
  if(method == 'weighted_catalog' | method == 'median_catalog'){
    probs_all    <- apply(genome_matrix, 1, 
                          function(x, y){calc_prob(x, y)$probs}, y = signatures)
    max_sigs_all <- apply(genome_matrix, 1, 
                          function(x, y){calc_prob(x, y)$sig_max}, y = signatures)
    max_val_all  <- apply(genome_matrix, 1, 
                          function(x, y){calc_prob(x, y)$max_val}, y = signatures)
  }
  if(method == 'cosine_simil'){
    simils_all    <- apply(genome_matrix, 1, 
                          function(x, y){calc_cos(x, y)$simils}, y = signatures)
    max_sigs_cos_all <- apply(genome_matrix, 1, 
                          function(x, y){calc_cos(x, y)$sig_max}, y = signatures)
    max_val_cos_all  <- apply(genome_matrix, 1, 
                          function(x, y){calc_cos(x, y)$max_val}, y = signatures) 
  }
  

  if(method == 'weighted_catalog' | method == 'median_catalog'){
    output <- cbind(t(probs_all), 
                    max_sigs_all, 
                    max_val_all)
    if(method == 'weighted_catalog') mname = 'wl'
    if(method == 'median_catalog') mname = 'ml'
    colnames(output) <- c(paste0(colnames(signatures), '_', mname), 
                          paste0('sig_max_', mname),
                          paste0('max_', mname))
  }

  if(method == 'cosine_simil'){ 
    output <- cbind(t(simils_all),
                    max_sigs_cos_all,
                    max_val_cos_all)
    
    colnames(output) <- c(paste0(colnames(signatures), '_c'), 
                          'sig_max_c',
                          'max_c')
  }

  return(output)
}