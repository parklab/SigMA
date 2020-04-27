#' Decomposes the mutational spectrum of a genome in terms of
#' tumor type specific signatures that were calculated through analysis
#' of public WGS samples from ICGC and TCGA consortia, and contained
#' as a list in the package. Non-negative-least squares algorithm
#' is used and the number of signatures to be considered in the 
#' decomposition is increased gradually, first all pairs from among the
#' available signatures are considered and minimal error pair is kept. 
#' Then all 3-signature combinations, 4-signature combinations and so on
#' are considered. The result is updated if the error is smaller with
#' larger number of signatures 
#'
#' @param spect composite spectrum that is being decomposed
#' @param signatures a data.frame that contains the signatures in its
#' columns
#' @param data sequencing platform that as in run(), used for setting
#' the maximum number of signatures that is allowed in the decomposition

decompose <- function(spect, signatures, data, nloop_user = NULL, delta=0){
  # calculates frobenius error
  error <- function(spect, signatures, exposures){
    reco <- (as.matrix(signatures) %*% exposures)
    error_frac <- sqrt(sum(diag((spect - reco) %*% t(spect - reco))))/sqrt(sum(diag(spect %*% t(spect))))
    return(error_frac)
  }

  dim2 <- dim(signatures)[[2]]
  if(dim2 == 2){
    exps <- coef(nnls::nnls(as.matrix(signatures), spect))
    inds <- which(exps != 0)
    exps <- exps[inds]
    sigs <- colnames(signatures)[inds]

    error <- error(spect, signatures[, sigs], exps)
    return(list(signatures = sigs,
                exposures = exps,
                error = error))
  }

  min_error <- 10000000
  min_indices <- NULL
  min_exposures <- NULL

  if(!is.null(nloop_user)) nloop <- min(dim(signatures)[[2]], nloop_user)
  else{
    nloop <- min(dim(signatures)[[2]], 8) 
    # increasing nloop to too large values causes overfitting
    # maximum 5 signatures are considered for panels and WES 
    if(data == "wgs_pancan"){ 
       nloop <- 6
    }
    else if(data == "tcga_mc3"){
      if(dim(signatures)[[2]] > 10){ 
        nloop <- 8
      }
      else{
        nloop <- min(dim(signatures)[[2]], 7)
      }
    }
    else if(data == "seqcap" | data == "seqcap_probe"){
      nloop <- min(dim(signatures)[[2]], 5)
    }
    else if(data == "msk" | data == "fo" | data == "op")
      nloop <- min(nloop, 5)
  }

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
                      if(nloop >= 7){
                        for(o in (n+1):(dim2 - nloop + 7)){
                          if(nloop >= 8){
                            for(p in (o+1):(dim2 - nloop + 8)){
                              if(nloop >= 9){
                                for(r in (p+1):(dim2 - nloop + 9)){
                                  exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m, n, o, p, r)]), spect))
                                  error_this <- error(spect, as.matrix(signatures[, c(i, j, k, l, m, n, o, p, r)]), exposures)
                                  if(error_this < min_error - delta){
                                    min_error <- error_this
                                    min_indices <- c(i, j, k, l, m, n, o, p, r)
                                    min_exposures <- exposures
                                  }
                                }
                              }
                              else{
                              exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m, n, o, p)]), spect))
                              error_this <- error(spect, as.matrix(signatures[, c(i, j, k, l, m, n, o, p)]), exposures)
                              if(error_this < min_error - delta){
                                min_error <- error_this
                                min_indices <- c(i, j, k, l, m, n, o, p)
                                min_exposures <- exposures
                              }
                              }
                            }
                          }
                          else{
                            exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m, n, o)]), spect))
                            error_this <- error(spect, as.matrix(signatures[, c(i, j, k, l, m, n, o)]), exposures)
                            if(error_this < min_error - delta){
                              min_error <- error_this
                              min_indices <- c(i, j, k, l, m, n, o)
                              min_exposures <- exposures
                            }
                          }
                        }
                      }else{
                        exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m, n)]), spect))
                        error_this <- error(spect, as.matrix(signatures[, c(i, j, k, l, m, n)]), exposures)
                        if(error_this < min_error - delta){
                          min_error <- error_this
                          min_indices <- c(i, j, k, l, m, n)
                          min_exposures <- exposures
                        }
                      }
                    }
                  }else{
                    exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l, m)]), spect))
                    error_this <- error(spect, as.matrix(signatures[, c(i, j, k, l, m)]), exposures)
                    if(error_this < min_error - delta){
                      min_error <- error_this
                      min_indices <- c(i, j, k, l, m)
                      min_exposures <- exposures
                    }
                  }
                }
              }else{
                exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k, l)]), spect))
                error_this <- error(spect, as.matrix(signatures[, c(i, j, k, l)]), exposures)
                if(error_this < min_error - delta){
                  min_error <- error_this
                  min_indices <- c(i, j, k, l)
                  min_exposures <- exposures
                }
              }
            }
          }else{
            exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k)]), spect))
            error_this <- error(spect, as.matrix(signatures[, c(i, j, k)]), exposures)
            if(error_this < min_error - delta){
              min_error <- error_this
              min_indices <- c(i, j, k)
              min_exposures <- exposures
            }
          }
        }
      }else{
        exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j)]), spect))
        error_this <- error(spect, as.matrix(signatures[, c(i, j)]), exposures)
        if(error_this < min_error - delta){
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
