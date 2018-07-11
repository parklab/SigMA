# decomposes the mutational spectrum of a genome using the signatures
# provided in the table, different number of signatures are tried
decompose <- function(this_genome, signatures, data){

  # calculates frobenius error
  error <- function(this_genome, signatures, exposures){
    reco <- (as.matrix(signatures) %*% exposures)
    error_frac <- sqrt(c(this_genome - reco) %*% c(this_genome - reco))/sqrt(this_genome %*% this_genome)
    return(error_frac)
  }

  dim2 <- dim(signatures)[[2]]
  if(dim2 == 2){
    exps <- coef(nnls::nnls(as.matrix(signatures), this_genome))
    inds <- which(exps != 0)
    exps <- exps[inds]
    sigs <- colnames(signatures)[inds]

    error <- error(this_genome, signatures[, sigs], exps)
    return(list(signatures = sigs,
                exposures = exps,
                error = error))
  }

  min_error <- 1
  min_indices <- NULL
  min_exposures <- NULL

  nloop <- min(dim(signatures)[[2]], 6)
  # increasing nloop to too large values causes overfitting
  # maximum 5 signatures are considered for panels and WES
  if(data == "msk" | data == "found" | data == "seqcap")
    nloop <- min(nloop, 5)

  if(nloop < 2)
    stop('linear decomposition requires at least 2 signatures')
  for(i in 1:(dim2 - nloop + 1)){
    for(j in (i+1):(dim2 - nloop + 2)){
      if(nloop >= 3){
        for(k in (j+1):(dim2 - nloop + 3)){
          if(nloop >= 4){
            for(l in (k+1):(dim2 - nloop + 4)){
              if(nloop >= 5){
                for(m in (l+1):(dim2 - nloop + 5)){
                  if(nloop >= 6){
                    for(n in (m+1):(dim2 - nloop + 6)){
                      exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m, n)]), this_genome))
                      error_this <- error(this_genome, as.matrix(signatures[, c(i, j, k, l, m, n)]), exposures)
                      if(error_this < min_error){
                        min_error <- error_this
                        min_indices <- c(i, j, k, l, m, n)
                        min_exposures <- exposures
                      }
                    }
                  }else{
                    exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m)]), this_genome))
                    error_this <- error(this_genome, as.matrix(signatures[, c(i, j, k, l, m)]), exposures)
                    if(error_this < min_error){
                      min_error <- error_this
                      min_indices <- c(i, j, k, l, m)
                      min_exposures <- exposures
                    }
                  }
                }
              }else{
                exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l)]), this_genome))
                error_this <- error(this_genome, as.matrix(signatures[, c(i, j, k, l)]), exposures)
                if(error_this < min_error){
                  min_error <- error_this
                  min_indices <- c(i, j, k, l)
                  min_exposures <- exposures
                }
              }
            }
          }else{
            exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k)]), this_genome))
            error_this <- error(this_genome, as.matrix(signatures[, c(i, j, k)]), exposures)
            if(error_this < min_error){
              min_error <- error_this
              min_indices <- c(i, j, k)
              min_exposures <- exposures
            }
          }
        }
      }else{
        exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j)]), this_genome))
        error_this <- error(this_genome, as.matrix(signatures[, c(i, j)]), exposures)
        if(error_this < min_error){
          min_error <- error_this
          min_indices <- c(i, j)
          min_exposures <- exposures
        }
      }
    }
  }
  return(list(signatures = colnames(signatures)[min_indices],
              exposures = min_exposures,
              error = min_error))
}
