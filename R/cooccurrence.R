#' Calculate cooccurence frequencies and relative ratios of signatures 
#' in a given training data
#' 
#' @param df an input dataset with signature assignments
#' @param signames vector of signature names that are of interest
#' @return list of two matrices cooccur_amp_matrix, mean ratio between
#' amplitudes of two signatures when they are nonzero,
#' cooccur_freq_matrix, the fraction of samples that have both signature

cooccurrence <- function(df, signames){
  df <- df[,match(signames, colnames(df))]

  cooccur_freq_matrix <- matrix(0, length(signames), length(signames))
  cooccur_amp_matrix <- matrix(0, length(signames), length(signames))

  for(isig in 1:length(signames)){
    coexists <- apply(df, 1, 
                      function(x, isig){
                        return(x[[isig]] > 0 & x > 0)
                      }, isig = isig)
    coexists <- matrix(unlist(coexists), byrow = T, nrow = dim(df)[[1]])


    #add if else for cases where only one sample has nonzero value for the 
    #signature with index isig
    if(sum(coexists[, isig]) > 1){
      coexists_sum <- colSums(coexists[coexists[, isig], ])
      count_nonzero <- sum(coexists[, isig] > 0)
      cooccur_freq_matrix[isig, ] <- coexists_sum / count_nonzero
    }else{
      cooccur_freq_matrix[isig, ] <- as.integer(coexists[coexists[, isig], ])
    }
    
    cooccur_rat <- apply(df[coexists[, isig], ], 1,
                         function(x, isig){
                           x/x[[isig]]
                         }, isig = isig)
    rats <- unlist(apply(cooccur_rat, 1,
                         function(x){
                           median(x[x > 0])  
                         }))
    rats[is.na(rats)] <- 0
    cooccur_amp_matrix[isig, ] <- rats
  }
  colnames(cooccur_amp_matrix) <- signames
  rownames(cooccur_amp_matrix) <- signames

  colnames(cooccur_freq_matrix) <- signames
  rownames(cooccur_freq_matrix) <- signames

  return(list(cooccur_amp_matrix = cooccur_amp_matrix,
              cooccur_freq_matrix = cooccur_freq_matrix))
}