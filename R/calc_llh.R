#' Calculates likelihood of the genome with respect to the available 
#' signature probability distributions
#' 
#' @param spectrum is the mutational spectrum
#' @param signatures is the reference signature catalog with the 
#' probability distributions
#' @param counts is the number of cases in each cluster that is
#' represented in the catalog. They are used as weights for each
#' signature in the catalog
#' @param normalize is true by default only for when it is used
#' together with NNLS in the match_to_catalog function it is not 
#' normalized here but outside of the function

calc_llh <- function(spectrum, 
                     signatures, 
                     counts = NULL, 
                     normalize = T){

  eps <- 0.00001
  probs <- apply(signatures, 2,
                 function(x, pow){
                   if(sum(is.na(x)) > 0) return(0)
                   prob <- 0
                   for(i in 1:length(x)){
                     if(x[[i]] == 0) x[[i]] <- eps
                     if(pow[[i]] > 0){
                       prob <- prob + pow[[i]] * log(x[[i]])
                     }
                   }
                   return(prob)
                 },
                 pow = spectrum)

  if(normalize){ # normalizes the total prob for all clusters to be 1
    if(is.null(counts) | length(counts) != dim(signatures)[[2]]){
      mean_probs <- mean(probs)

      q1_probs <- mean(probs[probs < mean(probs)])
      max_probs <- max(probs)

      # often having too small values leads to infinities so
      # only the probabilities in the first quartile are kept
      # to be nonzero
      inds_keep <- which(probs >= q1_probs)
      inds_rm <- which(probs < q1_probs)
      probs[inds_rm] <- 0

      probs[inds_keep] <- exp(probs[inds_keep] - max_probs)
      probs[inds_keep] <- probs[inds_keep]/sum(probs[inds_keep])
    }
    else{
      mean_probs <- mean(probs)

      q1_probs <- mean(probs[probs < mean(probs)])
      max_probs <- max(probs)

      # often having too small values leads to infinities so
      # only the probabilities in the first quartile are kept
      # to be nonzero
      inds_keep <- which(probs >= q1_probs)
      inds_rm <- which(probs < q1_probs)
      probs[inds_rm] <- 0
      probs[inds_keep] <- exp(probs[inds_keep] - max_probs)
      probs <- as.numeric(unlist(probs))*as.numeric(unlist(counts))
      probs[inds_keep] <- probs[inds_keep]/sum(probs[inds_keep])
    }
  }
 
  ind_max <- which(max(probs) == probs)
 
  # signature with highest likelihood
  if(length(ind_max) > 0) ind_max <- ind_max[[1]]
  sig_max <- colnames(signatures)[[ind_max]]
  max_val <- probs[[ind_max]]

  names(probs) <- colnames(signatures)

  return(list(probs = probs, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
}

