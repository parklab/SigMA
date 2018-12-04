#' calculates cosine similarity between the spectrum and 
#' a set of signatures
#'
#' @param spectrum is the mutational spectrum
#' @param signatures is the reference signature catalog with
#' the probability distribution

# cosine similarity
cosine <- function(x, y){
  return(x%*%y / sqrt(x%*%x * y%*%y))
}

# calculates cosine similarity for all signature options
calc_cos <- function(spectrum, signatures){
  simils <- apply(signatures, 2,
                  function(x){ cosine(spectrum, x) })

  ind_max <- which(max(simils) == simils)

  # signature with maximum cosine similarity and its value
  if(length(ind_max) > 0) ind_max = ind_max[[1]]
  sig_max <- colnames(signatures)[[ind_max]]
  max_val <- simils[[ind_max]]

  return(list(simils = simils,
              ind_max = ind_max,
              sig_max = sig_max,
              max_val = max_val))
}
