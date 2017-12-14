#' calculates the weights for each 96 dimensional mutation types according
#' to how distinctive this context is for that signature
#' 
#' @param signatures table of signatures
#' @param cooccur_freq_matrix symmetric matrix telling the cooccurence
#' frequencies of two signatures
#' @param cooccur_median_matrix symmetrix matrix telling the ratio of 
#' a signature pair when they cooccur
#' @return the weight that signifies the importance of that particular
#' type for each column

choose_optimal_context <- function(signatures, 
                                   cooccur_freq_matrix, 
                                   cooccur_median_matrix){

  print(colnames(cooccur_freq_matrix))
  print(colnames(cooccur_median_matrix))
  print(colnames(signatures))
  print(sum(colnames(cooccur_freq_matrix) == colnames(signatures)))
  print(sum(colnames(cooccur_median_matrix) == colnames(signatures)))
  print(dim(signatures))
  if(sum(colnames(cooccur_freq_matrix) == colnames(signatures)) != dim(signatures)[[2]] |  
    sum(colnames(cooccur_median_matrix) == colnames(signatures)) != dim(signatures)[[2]]) 
    stop('signatures in cooccurrence matrices and the signature table are incompatible') 
 
 
  signames <- colnames(signatures)
  most_distinct <- matrix(0, dim(signatures)[[1]], dim(signatures)[[2]])
  for(isig in 1:length(signames)){
    others <- signatures[, -isig]
    others_contrib <- apply(others, 1, 
                            function(x, cooccur_freq, cooccur_median, isig){
                              coeff <- cooccur_freq[isig, -isig] * cooccur_median[isig, -isig]
                              return(coeff %*% x)
                            }, 
                            cooccur_freq = cooccur_freq_matrix, 
                            cooccur_median = cooccur_median_matrix,
                            isig = isig)


    others_contrib <- c(unlist(others_contrib))

    print(length(1/others_contrib))
    print(length(signatures[,isig]))
    most_distinct[, isig] <- signatures[, isig] * (1/others_contrib)
  }
  colnames(most_distinct) <- signames
  return(most_distinct)
}